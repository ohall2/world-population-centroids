#!/usr/bin/env python3
"""
World Population Centroids v2 – befolkningscentroider från folkräkningar
(facit för valideringen).

Data (tools/download_validation.py):
  census   U.S. Census Bureau, befolkning och mittpunkt per kvartersgrupp,
           2000, 2010, 2020 (USA och Puerto Rico).
  geostat  Eurostat GEOSTAT, befolkning per ruta om 1 km² (ETRS89-LAEA),
           2011 och 2021. Rutans mittpunkt används. För 2011 används bara
           rutor som bygger på aggregerade folkräkningsdata (METHD_CL = "A");
           modellerade rutor ("D", "M") och JRC-filen tas inte med.
           Befolkning som inte har kunnat placeras i en ruta (2021,
           "CC_unallocated") kan inte användas; andelen redovisas.

För varje GADM-enhet (nivå 0–2) beräknas befolkningsviktade centroider med
samma metoder som för nattljuset (wpc.centroids): Aboufadel–Austin
(3D-medelvärde) och Barmore (iterativ, med start i Aboufadel–Austin). Vikten
är befolkningen i punkten; ingen ytviktning behövs, eftersom vikten redan är
ett antal personer.

En enhet räknas som fullständigt täckt om ingen av dess adm2-enheter med
nattljus (DMSP 2010, sum_light > 0) saknar befolkning i facit. Det fångar
delar av länder som rutnätet inte täcker (t.ex. Kanarieöarna) och länder
där bara modellerade rutor finns. Kolumnen `complete` anger detta.

Utdata: data/validation/<data>/<år>/adm<nivå>.csv med kolumnerna
unit_id, level, year, dataset, pop, n_points, aa_lon, aa_lat, bm_lon,
bm_lat, complete, samt summary.json per år.

Exempel
  python tools/validation_centroids.py --dataset census --years 2000 2010 2020
  python tools/validation_centroids.py --dataset geostat --years 2011 2021
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import time
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.centroids import aboufadel_austin, barmore  # noqa: E402
from wpc.progress import Progress, log  # noqa: E402
from wpc.run import load_units  # noqa: E402

CENSUS_COLS = ["STATEFP", "COUNTYFP", "TRACTCE", "BLKGRPCE", "POPULATION", "LATITUDE", "LONGITUDE"]


def read_census(root: Path, year: int) -> pd.DataFrame:
    files = sorted((root / str(year)).glob("*.txt"))
    if not files:
        sys.exit(f"Inga Census-filer för {year} i {root}")
    parts = []
    for f in files:
        raw = f.read_bytes().decode("utf-8-sig", "replace")
        header = 0 if raw.startswith("STATEFP") else None
        d = pd.read_csv(io.StringIO(raw), header=header, names=None if header == 0 else CENSUS_COLS,
                        dtype={"STATEFP": str, "COUNTYFP": str, "TRACTCE": str})
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    return pd.DataFrame({"lon": d.LONGITUDE.astype(float), "lat": d.LATITUDE.astype(float),
                         "pop": d.POPULATION.astype(float)}), {"files": len(files)}


def _laea_centres(ids: pd.Series, pattern: str, scale: float) -> tuple[np.ndarray, np.ndarray]:
    m = ids.str.extract(pattern)
    n = m[0].astype(float) * scale + 500.0
    e = m[1].astype(float) * scale + 500.0
    tr = Transformer.from_crs(3035, 4326, always_xy=True)
    lon, lat = tr.transform(e.values, n.values)
    return np.asarray(lon), np.asarray(lat)


def read_geostat(root: Path, year: int):
    zips = sorted((root / str(year)).glob("*.zip"))
    if not zips:
        sys.exit(f"Ingen GEOSTAT-fil för {year} i {root}")
    zf = zipfile.ZipFile(zips[0])
    info = {"file": zips[0].name}
    if year == 2011:
        name = next(n for n in zf.namelist() if re.search(r"GEOSTAT_grid_POP_1K_2011_V2_0_1\.csv$", n))
        d = pd.read_csv(zf.open(name), dtype={"GRD_ID": str, "CNTR_CODE": str, "METHD_CL": str})
        info["method_by_country"] = {c: g.METHD_CL.value_counts().to_dict()
                                     for c, g in d.groupby("CNTR_CODE")}
        info["pop_total"] = float(d.TOT_P.sum())
        d = d[d.METHD_CL == "A"]
        info["pop_used"] = float(d.TOT_P.sum())
        lon, lat = _laea_centres(d.GRD_ID, r"N(\d+)E(\d+)", 1000.0)
        pop = d.TOT_P.astype(float).values
    else:
        name = next(n for n in zf.namelist() if re.search(r"ESTAT_Census_2021_V\d+\.csv$", n))
        d = pd.read_csv(zf.open(name), usecols=["GRD_ID", "T", "CNTR_ID"],
                        dtype={"GRD_ID": str, "CNTR_ID": str})
        d = d[d["T"] > 0]
        un = d.GRD_ID.str.endswith("_unallocated")
        info["pop_total"] = float(d["T"].sum())
        info["unallocated_by_country"] = d[un].set_index("GRD_ID")["T"].astype(float).to_dict()
        d = d[~un]
        info["pop_used"] = float(d["T"].sum())
        lon, lat = _laea_centres(d.GRD_ID, r"N(\d+)E(\d+)$", 1.0)
        pop = d["T"].astype(float).values
    return pd.DataFrame({"lon": lon, "lat": lat, "pop": pop}), info


def assign(points: pd.DataFrame, units: gpd.GeoDataFrame) -> np.ndarray:
    """Index i `units` för varje punkt (−1 om punkten ligger utanför alla)."""
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(points.lon, points.lat), crs=4326)
    j = gpd.sjoin(pts, units[["geometry"]].reset_index(drop=True), how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    return j["index_right"].fillna(-1).astype(int).values


def centroids(points: pd.DataFrame, idx: np.ndarray, units: gpd.GeoDataFrame, label: str):
    rows = []
    order = np.argsort(idx, kind="stable")
    sidx = idx[order]
    bounds = np.searchsorted(sidx, np.arange(len(units) + 1))
    lon, lat, w = (points[c].values[order] for c in ("lon", "lat", "pop"))
    with Progress(label, total=len(units)) as prog:
        for u in range(len(units)):
            a, b = bounds[u], bounds[u + 1]
            if b > a and w[a:b].sum() > 0:
                aa = aboufadel_austin(lon[a:b], lat[a:b], w[a:b])
                bm = barmore(lon[a:b], lat[a:b], w[a:b], start=aa)
                rows.append((units.unit_id.iloc[u], float(w[a:b].sum()), int(b - a), *aa, bm[0], bm[1]))
            prog.advance()
    return pd.DataFrame(rows, columns=["unit_id", "pop", "n_points", "aa_lon", "aa_lat", "bm_lon", "bm_lat"])


def main(argv=None):
    p = argparse.ArgumentParser(description="Befolkningscentroider från folkräkningar.")
    p.add_argument("--dataset", choices=["census", "geostat"], required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--raw", type=Path, default=Path("data/raw/validation"))
    p.add_argument("--out", type=Path, default=Path("data/validation"))
    p.add_argument("--boundaries", type=Path, default=Path("data/raw/gadm41/gadm_410-levels.gpkg"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--light-ref", type=Path,
                   default=Path("data/results_mask2km/gadm41/dmsp-irqr-reg/adm2/2010.csv"),
                   help="Nattljus per adm2 för kontrollen av täckning.")
    args = p.parse_args(argv)

    lit = pd.read_csv(args.light_ref, keep_default_na=False, usecols=["unit_id", "sum_light"])
    lit = set(lit.unit_id[pd.to_numeric(lit.sum_light, errors="coerce") > 0])
    units = {lv: load_units(args.boundaries, "gadm41", lv, args.cache_dir) for lv in (0, 1, 2)}
    parent = pd.DataFrame({"unit_id": units[2].unit_id.values, "adm1": units[2].GID_1.values,
                           "adm0": units[2].GID_0.values})

    for year in args.years:
        t0 = time.time()
        reader = read_census if args.dataset == "census" else read_geostat
        with Progress(f"Läser {args.dataset} {year}"):
            pts, info = reader(args.raw / args.dataset, year)
        pts = pts[(pts["pop"] > 0) & np.isfinite(pts.lon) & np.isfinite(pts.lat)].reset_index(drop=True)
        log(f"{args.dataset} {year}: {len(pts)} punkter, {pts['pop'].sum():,.0f} personer")
        od = args.out / args.dataset / str(year)
        od.mkdir(parents=True, exist_ok=True)
        summary = {"dataset": args.dataset, "year": year, "points": int(len(pts)),
                   "population": float(pts["pop"].sum()), **info, "levels": {}}
        for lv in (2, 1, 0):
            with Progress(f"adm{lv}: punkter till enheter"):
                idx = assign(pts, units[lv])
            outside = pts["pop"].values[idx < 0].sum()
            res = centroids(pts, idx, units[lv], f"adm{lv} {year}")
            if lv == 2:
                have = set(res.unit_id)
                bad = parent[parent.unit_id.isin(lit) & ~parent.unit_id.isin(have)]
                bad1, bad0 = set(bad.adm1), set(bad.adm0)
                res["complete"] = True
                summary["incomplete_adm2_by_country"] = {
                    k: int(v) for k, v in bad[bad.adm0.isin(res.unit_id.str[:3])].adm0.value_counts().items()}
            elif lv == 1:
                res["complete"] = ~res.unit_id.isin(bad1)
            else:
                res["complete"] = ~res.unit_id.isin(bad0)
            res.insert(1, "level", lv)
            res.insert(2, "year", year)
            res.insert(3, "dataset", args.dataset)
            res.to_csv(od / f"adm{lv}.csv", index=False, float_format="%.8f")
            summary["levels"][f"adm{lv}"] = {"units": int(len(res)), "complete": int(res.complete.sum()),
                                            "pop_outside_units": float(outside)}
            log(f"  adm{lv}: {len(res)} enheter ({int(res.complete.sum())} fullständiga), "
                f"{outside:,.0f} personer utanför GADM-enheterna")
        (od / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str) + "\n")
        log(f"{args.dataset} {year} klart på {time.time() - t0:.0f} s -> {od}")


if __name__ == "__main__":
    main()
