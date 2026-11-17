#!/usr/bin/env python3
"""
World Population Centroids v2 – befolkningscentroider från folkräkningar
(facit för valideringen).

Data (tools/download_validation.py):
  inegi    INEGI ITER, befolkning och läge per ort i Mexikos folkräkningar
           2010 och 2020. Summaraderna (land, delstat, kommun, och 9998/9999
           för orter med en eller två bostäder) tas bort; de orterna finns
           också var för sig med läge, så ingen befolkning saknas. Stora
           städer är ofta en enda ort per kommun, så facit är grovt för
           kommuner med stora städer.
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

Fullständig täckning (kolumnen `complete`, och `missing_light_share`):
andelen av enhetens nattljus (DMSP 2010, samregistrerad) som ligger i
underenheter på den finaste GADM-nivå som finns för landet och som saknar
befolkning i facit. En adm2-enhet är täckt om den har befolkning; för en
adm1-enhet utan adm2-underenheter (t.ex. Moldavien, Andorra) används
enheten själv, och likadant för länder utan adm1. En enhet räknas som
fullständig om andelen är under 1 %. Det fångar delar av länder som
rutnätet inte täcker (t.ex. Kanarieöarna), länder där bara modellerade
rutor finns och grannländer som bara fått några randrutor, men utesluter
inte ett land för att några få småenheter saknar invånare.

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

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.progress import Progress, log  # noqa: E402

# geopandas, pyproj och wpc.centroids/wpc.run importeras i de funktioner som
# behöver dem, så att läsarna kan testas där bara pandas finns.

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
    # 2000 års filer anger "+.,-." som läge för kvartersgrupper utan invånare.
    lon = pd.to_numeric(d.LONGITUDE, errors="coerce")
    lat = pd.to_numeric(d.LATITUDE, errors="coerce")
    pop = pd.to_numeric(d.POPULATION, errors="coerce")
    nopos = lon.isna() | lat.isna()
    if (pop[nopos] > 0).any():
        raise ValueError(f"{int((pop[nopos] > 0).sum())} kvartersgrupper med invånare saknar läge")
    return (pd.DataFrame({"lon": lon, "lat": lat, "pop": pop}),
            {"files": len(files), "block_groups": int(len(d)), "without_position_pop0": int(nopos.sum())})


def _dms(v) -> float:
    """INEGI-läge -> decimalgrader.

    2010: heltal GGGMMSS / GGMMSS (longitud väst och latitud nord utan tecken).
    2020: text som 102°17'45.768" W.
    """
    if pd.isna(v) or str(v).strip() == "":
        return np.nan
    t = str(v).strip()
    m = re.match(r"(\d+)°(\d+)'([\d.]+)\"?\s*([NSEWO])", t)
    if m:
        d = int(m[1]) + int(m[2]) / 60 + float(m[3]) / 3600
        return -d if m[4] in "SWO" else d
    if t.isdigit():
        n = int(t)
        return n // 10000 + (n // 100 % 100) / 60 + (n % 100) / 3600
    raise ValueError(f"Okänt INEGI-läge: {t!r}")


def read_inegi(root: Path, year: int):
    zips = sorted((root / str(year)).glob("*.zip"))
    if not zips:
        sys.exit(f"Ingen INEGI-fil för {year} i {root}")
    zf = zipfile.ZipFile(zips[0])
    name = next(n for n in zf.namelist() if "conjunto_de_datos/" in n and n.endswith(".csv"))
    d = pd.read_csv(zf.open(name), encoding="utf-8-sig", dtype=str, low_memory=False)
    d.columns = [c.lower() for c in d.columns]
    loc = d["loc"]                      # obs: d.loc är pandas indexering, inte kolumnen
    total = pd.to_numeric(d.pobtot[(d.entidad == "00") & (loc == "0000")], errors="coerce").iloc[0]
    agg = loc.isin(["0000", "9998", "9999"])
    small = pd.to_numeric(d.pobtot[(d.entidad != "00") & (d.mun != "000") & loc.isin(["9998", "9999"])],
                          errors="coerce").sum()
    d = d[~agg]
    lon = -d.longitud.map(_dms) if year == 2010 else d.longitud.map(_dms)
    lat = d.latitud.map(_dms)
    pop = pd.to_numeric(d.pobtot, errors="coerce")
    info = {"file": zips[0].name, "localities": int(len(d)), "pop_national_total": float(total),
            "pop_in_localities": float(pop.sum()),
            "pop_small_localities_summary_rows": float(small),
            "localities_without_position": int((lon.isna() | lat.isna()).sum())}
    return pd.DataFrame({"lon": lon.values, "lat": lat.values, "pop": pop.values}), info


def _laea_centres(ids: pd.Series, pattern: str, scale: float) -> tuple[np.ndarray, np.ndarray]:
    from pyproj import Transformer
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


def assign(points: pd.DataFrame, units) -> np.ndarray:
    """Index i `units` för varje punkt (−1 om punkten ligger utanför alla)."""
    import geopandas as gpd
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(points.lon, points.lat), crs=4326)
    j = gpd.sjoin(pts, units[["geometry"]].reset_index(drop=True), how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    return j["index_right"].fillna(-1).astype(int).values


def centroids(points: pd.DataFrame, idx: np.ndarray, units, label: str):
    from wpc.centroids import aboufadel_austin, barmore
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


def missing_light_share(have: dict, light: dict, u1: pd.DataFrame, u2: pd.DataFrame) -> dict:
    """Andel av enhetens ljus i underenheter utan befolkning i facit (se ovan)."""
    l2 = u2.assign(light=u2.unit_id.map(light[2]).fillna(0.0).values,
                   miss=~u2.unit_id.isin(have[2]))
    m2 = pd.Series(np.where(l2.miss, 1.0, 0.0), index=l2.unit_id)
    g = l2.groupby("adm1")
    tot1 = g.light.sum()
    mis1 = l2[l2.miss].groupby("adm1").light.sum().reindex(tot1.index).fillna(0.0)
    share1 = (mis1 / tot1.where(tot1 > 0)).fillna(0.0)
    m1 = pd.Series(np.where(~u1.unit_id.isin(have[1]), 1.0, 0.0), index=u1.unit_id)
    m1.update(share1)                    # adm1 med adm2-underenheter: andel av ljuset
    l1 = u1.assign(light=u1.unit_id.map(light[1]).fillna(0.0).values, share=m1.values)
    tot0 = l1.groupby("adm0").light.sum()
    mis0 = (l1.light * l1.share).groupby(l1.adm0).sum()
    share0 = (mis0 / tot0.where(tot0 > 0)).fillna(0.0)
    m0 = pd.Series({u: (0.0 if u in have[0] else 1.0) for u in light[0].index})
    m0.update(share0)
    return {2: m2, 1: m1, 0: m0}


def main(argv=None):
    p = argparse.ArgumentParser(description="Befolkningscentroider från folkräkningar.")
    p.add_argument("--dataset", choices=["census", "geostat", "inegi"], required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--raw", type=Path, default=Path("data/raw/validation"))
    p.add_argument("--out", type=Path, default=Path("data/validation"))
    p.add_argument("--boundaries", type=Path, default=Path("data/raw/gadm41/gadm_410-levels.gpkg"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--light-ref", type=Path,
                   default=Path("data/results_mask2km/gadm41/dmsp-irqr-reg"),
                   help="Nattljus per enhet (adm0–2, år 2010) för kontrollen av täckning.")
    p.add_argument("--max-missing", type=float, default=0.01)
    args = p.parse_args(argv)

    from wpc.run import load_units
    units = {lv: load_units(args.boundaries, "gadm41", lv, args.cache_dir) for lv in (0, 1, 2)}
    light = {}
    for lv in (0, 1, 2):
        d = pd.read_csv(args.light_ref / f"adm{lv}" / "2010.csv", keep_default_na=False,
                        usecols=["unit_id", "sum_light"])
        light[lv] = pd.to_numeric(d.set_index("unit_id").sum_light, errors="coerce").fillna(0.0)
    u1 = pd.DataFrame({"unit_id": units[1].unit_id.values, "adm0": units[1].GID_0.values})
    u2 = pd.DataFrame({"unit_id": units[2].unit_id.values, "adm1": units[2].GID_1.values,
                       "adm0": units[2].GID_0.values})

    for year in args.years:
        t0 = time.time()
        reader = {"census": read_census, "geostat": read_geostat, "inegi": read_inegi}[args.dataset]
        with Progress(f"Läser {args.dataset} {year}"):
            pts, info = reader(args.raw / args.dataset, year)
        pts = pts[(pts["pop"] > 0) & np.isfinite(pts.lon) & np.isfinite(pts.lat)].reset_index(drop=True)
        log(f"{args.dataset} {year}: {len(pts)} punkter, {pts['pop'].sum():,.0f} personer")
        od = args.out / args.dataset / str(year)
        od.mkdir(parents=True, exist_ok=True)
        summary = {"dataset": args.dataset, "year": year, "points": int(len(pts)),
                   "population": float(pts["pop"].sum()), **info, "levels": {}}
        have, results = {}, {}
        for lv in (2, 1, 0):
            with Progress(f"adm{lv}: punkter till enheter"):
                idx = assign(pts, units[lv])
            outside = pts["pop"].values[idx < 0].sum()
            res = centroids(pts, idx, units[lv], f"adm{lv} {year}")
            have[lv] = set(res.unit_id)
            results[lv] = res
            summary["levels"][f"adm{lv}"] = {"units": int(len(res)), "pop_outside_units": float(outside)}
        miss = missing_light_share(have, light, u1, u2)
        for lv in (2, 1, 0):
            res = results[lv]
            res["missing_light_share"] = res.unit_id.map(miss[lv]).fillna(0.0).values
            res["complete"] = res.missing_light_share < args.max_missing
            res.insert(1, "level", lv)
            res.insert(2, "year", year)
            res.insert(3, "dataset", args.dataset)
            res.to_csv(od / f"adm{lv}.csv", index=False, float_format="%.8f")
            summary["levels"][f"adm{lv}"]["complete"] = int(res.complete.sum())
            log(f"  adm{lv}: {len(res)} enheter ({int(res.complete.sum())} fullständiga), "
                f"{summary['levels'][f'adm{lv}']['pop_outside_units']:,.0f} personer utanför GADM-enheterna")
        summary["incomplete_adm0"] = {u: round(float(v), 4) for u, v in miss[0].items()
                                      if u in have[0] and v >= args.max_missing}
        (od / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str) + "\n")
        log(f"{args.dataset} {year} klart på {time.time() - t0:.0f} s -> {od}")


if __name__ == "__main__":
    main()
