"""Kör centroidberäkningen för valda år och nivåer.

Exempel (från projektmappen ~/wpc):
  python -m wpc.run --years all --levels 0 1 2

Varje kombination av gränsuppsättning, ljuskälla, nivå och år sparas som en
CSV-fil. Färdiga filer hoppas över, så en avbruten körning kan startas om med
samma kommando.

Utdata:  <out>/<gränser>/<källa>/adm<N>/<år>.csv  samt  <år>.json (körningens metadata)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely

from . import __version__
from .centroids import UnitResult, compute_unit
from .lights import SOURCES, available_years, make_source
from .progress import Progress, fmt_seconds, log


# Gränsuppsättningar: lager och ID-kolumner per nivå.
BOUNDARIES = {
    "gadm41": {"layer": "ADM_{level}", "id": "GID_{level}",
               "keep": {0: ["GID_0", "COUNTRY"],
                        1: ["GID_0", "COUNTRY", "GID_1", "NAME_1"],
                        2: ["GID_0", "COUNTRY", "GID_1", "NAME_1", "GID_2", "NAME_2"]}},
    "gadm28": {"layer": "adm{level}", "id": "OBJECTID",
               "keep": {0: ["OBJECTID", "ID_0", "ISO", "NAME_ENGLI"],
                        1: ["OBJECTID", "ID_0", "ISO", "NAME_0", "ID_1", "NAME_1"],
                        2: ["OBJECTID", "ID_0", "ISO", "NAME_0", "ID_1", "ID_2", "NAME_1", "NAME_2"]}},
}


def boundary_path(path: Path) -> str:
    """En zip-fil läses direkt; GDAL behöver då /vsizip/ och gpkg-filens namn."""
    if path.suffix == ".zip":
        import zipfile
        inner = [n for n in zipfile.ZipFile(path).namelist() if n.endswith(".gpkg")]
        if len(inner) != 1:
            raise ValueError(f"{path}: väntade en .gpkg i zip-filen, hittade {inner}")
        return f"/vsizip/{path}/{inner[0]}"
    return str(path)


def cache_path(cache_dir: Path, bset: str, level: int) -> Path:
    return cache_dir / f"{bset}_adm{level}.gpkg"


def build_cache(path: Path, bset: str, level: int, cache_dir: Path) -> None:
    """Läs ett lager ur gränsfilen (t.ex. en zip) och spara det okomprimerat.

    Att läsa en GeoPackage inuti en zip-fil är mycket långsamt, eftersom
    GDAL måste packa upp från början vid varje hopp i filen. Lagret sparas
    därför en gång som en egen GeoPackage med bara de kolumner som behövs.
    """
    spec = BOUNDARIES[bset]
    gdf = gpd.read_file(boundary_path(path), layer=spec["layer"].format(level=level),
                        columns=spec["keep"][level], engine="pyogrio")
    out = cache_path(cache_dir, bset, level)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part.gpkg")
    tmp.unlink(missing_ok=True)
    gdf.to_file(tmp, driver="GPKG", engine="pyogrio")
    tmp.replace(out)


def load_units(path: Path, bset: str, level: int, cache_dir: Path | None = None) -> gpd.GeoDataFrame:
    """Läs enheterna på en nivå, via cachen om den finns.

    Byggs cachen körs läsningen i en egen process, så att förloppsvisningen i
    huvudprocessen fortsätter att uppdateras (GDAL-läsningen blockerar annars
    Python och klockan står still).
    """
    spec = BOUNDARIES[bset]
    if cache_dir is not None:
        cached = cache_path(cache_dir, bset, level)
        if not cached.exists():
            with ProcessPoolExecutor(max_workers=1) as ex:
                ex.submit(build_cache, path, bset, level, cache_dir).result()
        gdf = gpd.read_file(cached, engine="pyogrio")
    else:
        gdf = gpd.read_file(boundary_path(path), layer=spec["layer"].format(level=level),
                            columns=spec["keep"][level], engine="pyogrio")
    gdf = gdf.rename(columns={spec["id"].format(level=level): "unit_id"}) if spec["id"] != "OBJECTID" \
        else gdf.assign(unit_id=gdf["OBJECTID"])
    if gdf["unit_id"].duplicated().any():
        raise ValueError(f"Dubbletter av enhets-ID på nivå {level}")
    return gdf


# --------------------------------------------------------------------------
# Arbetsprocesser
# --------------------------------------------------------------------------

_SOURCE = None
_SOURCE_KEY = None


def _get_source(key: tuple):
    global _SOURCE, _SOURCE_KEY
    if _SOURCE_KEY != key:
        if _SOURCE is not None:
            _SOURCE.close()
        name, year, raw_dir, calibration, flares = key
        _SOURCE = make_source(name, year, Path(raw_dir), Path(calibration),
                              Path(flares) if flares else None)
        _SOURCE_KEY = key
    return _SOURCE


def _work(args) -> tuple[int, dict]:
    idx, wkb, key = args
    geom = shapely.from_wkb(wkb) if wkb is not None else None
    return idx, asdict(compute_unit(geom, _get_source(key)))


# --------------------------------------------------------------------------
# Flaggor
# --------------------------------------------------------------------------

def flags(units: gpd.GeoDataFrame, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """0 = i egen enhet, 1 = utanför alla enheter, 2 = i en annan enhet, -1 = ingen centroid.

    Samma kodning som OUT_FLAG i v1.
    """
    out = np.full(len(units), -1, dtype="int8")
    have = ~np.isnan(lon)
    pts = shapely.points(lon[have], lat[have])
    own = units.geometry.values[have]
    inside = shapely.covers(own, pts)
    idx_have = np.flatnonzero(have)
    out[idx_have[inside]] = 0
    rest = idx_have[~inside]
    if len(rest):
        tree = shapely.STRtree(units.geometry.values)
        q = tree.query(shapely.points(lon[rest], lat[rest]), predicate="covered_by")
        hit = np.zeros(len(rest), dtype=bool)
        hit[np.unique(q[0])] = True
        out[rest] = np.where(hit, 2, 1)
    return out


# --------------------------------------------------------------------------
# Huvudflöde
# --------------------------------------------------------------------------

def run_one(units: gpd.GeoDataFrame, level: int, year: int, args, out_csv: Path) -> None:
    key = (args.source, year, str(args.raw_dir), str(args.calibration),
           str(args.flares) if args.flares else None)
    # Stora enheter först, så att inga långa jobb blir kvar i slutet.
    b = units.geometry.bounds
    order = np.argsort(-((b.maxx - b.minx) * (b.maxy - b.miny)).fillna(0).to_numpy())
    tasks = [(int(i), shapely.to_wkb(units.geometry.iloc[i]) if units.geometry.iloc[i] is not None else None, key)
             for i in order]
    results: list[dict | None] = [None] * len(units)
    t0 = time.time()
    with Progress(f"adm{level} {year}", total=len(tasks)) as p:
        if args.workers == 1:
            for t in tasks:
                i, r = _work(t)
                results[i] = r
                p.advance()
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                for i, r in pool.map(_work, tasks, chunksize=8):
                    results[i] = r
                    p.advance()

    df = pd.DataFrame(results)
    meta_cols = [c for c in units.columns if c != "geometry"]
    df = pd.concat([units[meta_cols].reset_index(drop=True), df], axis=1)
    df.insert(0, "year", year)
    df.insert(0, "level", level)
    df.insert(0, "source", args.source)
    df["aa_flag"] = flags(units, df["aa_lon"].to_numpy(), df["aa_lat"].to_numpy())
    df["bm_flag"] = flags(units, df["bm_lon"].to_numpy(), df["bm_lat"].to_numpy())

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_csv.with_suffix(".part")
    df.to_csv(tmp, index=False, float_format="%.8f")
    tmp.replace(out_csv)

    src = _get_source(key)
    meta = {"wpc_version": __version__, "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "seconds": round(time.time() - t0, 1), "level": level, "year": year,
            "boundaries": {"set": args.boundary_set, "path": str(args.boundaries)},
            "light_source": src.describe(), "units": int(len(df)),
            "with_centroid": int(df["aa_lon"].notna().sum()),
            "aa_flag_counts": {str(k): int(v) for k, v in df["aa_flag"].value_counts().sort_index().items()},
            "python": platform.python_version(), "platform": platform.platform(),
            "workers": args.workers}
    out_csv.with_suffix(".json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    log(f"adm{level} {year}: {meta['with_centroid']} av {meta['units']} enheter har centroid, "
        f"flaggor {meta['aa_flag_counts']} -> {out_csv}")


def parse_years(vals: list[str], available: list[int]) -> list[int]:
    if vals == ["all"]:
        return available
    years = []
    for v in vals:
        if "-" in v:
            a, b = v.split("-")
            years += list(range(int(a), int(b) + 1))
        else:
            years.append(int(v))
    return years


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Befolkningscentroider från nattljus (v2).")
    p.add_argument("--source", choices=SOURCES, default="dmsp-irqr",
                   help="Ljuskälla (standard dmsp-irqr).")
    p.add_argument("--years", nargs="+", default=["all"],
                   help="t.ex. all (alla år som finns för källan), 2010, 1992-1995")
    p.add_argument("--levels", nargs="+", type=int, default=[0, 1, 2])
    p.add_argument("--boundary-set", choices=list(BOUNDARIES), default="gadm41")
    p.add_argument("--boundaries", type=Path, default=Path("data/raw/gadm41/gadm_410-levels.gpkg"))
    p.add_argument("--raw-dir", type=Path, default=None,
                   help="Rådata (standard data/raw/dmsp, data/raw/ccnl eller data/raw/viirs).")
    p.add_argument("--calibration", type=Path, default=Path("config/calibration_v1.csv"))
    p.add_argument("--flares", type=Path, default=None,
                   help="Fackelmask. Standard: data/processed/flares_dmsp.gpkg (Elvidge m.fl. "
                        "2009) för dmsp-irqr och ccnl, data/processed/flares_viirs_r2km.gpkg "
                        "för viirs (se tools/prepare_flares_viirs.py). 'none' = ingen maskning.")
    p.add_argument("--out", type=Path, default=Path("data/results"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"),
                   help="Okomprimerade kopior av gränslagren (byggs automatiskt).")
    p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--limit", type=int, default=None, help="Bara de N första enheterna (för test).")
    args = p.parse_args(argv)

    if args.raw_dir is None:
        args.raw_dir = Path("data/raw") / ("dmsp" if args.source == "dmsp-irqr" else args.source)
    if args.flares is None:
        args.flares = Path("data/processed") / ("flares_viirs_r2km.gpkg" if args.source == "viirs"
                                                else "flares_dmsp.gpkg")
    if str(args.flares).lower() == "none":
        args.flares = None
    elif not args.flares.exists():
        sys.exit(f"Fackelmasken {args.flares} finns inte.")
    available = available_years(args.source, args.raw_dir, args.calibration)
    if not available:
        sys.exit(f"Inga filer för {args.source} i {args.raw_dir}.")
    years = parse_years(args.years, available)
    missing = sorted(set(years) - set(available))
    if missing:
        sys.exit(f"{args.source}: inga data för år {missing} i {args.raw_dir}.")
    log(f"wpc {__version__}: källa {args.source}, gränser {args.boundary_set}, nivåer {args.levels}, "
        f"år {years[0]}–{years[-1]} ({len(years)} st), {args.workers} parallella processer")
    t_all = time.time()
    for level in args.levels:
        units = None
        for year in years:
            out_csv = args.out / args.boundary_set / args.source / f"adm{level}" / f"{year}.csv"
            if args.limit is not None:
                out_csv = out_csv.with_name(f"{year}_test{args.limit}.csv")
            if out_csv.exists() and out_csv.with_suffix(".json").exists():
                log(f"adm{level} {year}: finns redan – hoppar över")
                continue
            if units is None:
                first = not cache_path(args.cache_dir, args.boundary_set, level).exists()
                text = (f"Förbereder gränser adm{level} (bara första gången, några minuter)"
                        if first else f"Läser gränser adm{level}")
                with Progress(text):
                    units = load_units(args.boundaries, args.boundary_set, level, args.cache_dir)
                    if args.limit is not None:
                        units = units.iloc[:args.limit].copy()
            run_one(units, level, year, args, out_csv)
    log(f"Klart på {fmt_seconds(time.time() - t_all)}.")


if __name__ == "__main__":
    main()
