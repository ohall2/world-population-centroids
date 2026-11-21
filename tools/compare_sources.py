#!/usr/bin/env python3
"""Centroider ur en jämförelsekälla för facitenheterna (analys B, planen §6).

Kör samma centroidkod (wpc.centroids.compute_unit) som utgåvan, men bara för de
enheter som finns i facit (complete, minst tio punkter) på nivå 1 och 2, och
bara för facitländerna. Varje facitår knyts till källans närmaste epok
(GHS-POP: 2000/2010/2020). Utdata i utgåvans format, en fil per epok:
  <out>/<källa>/adm<N>/<epok>.csv  med unit_id, year, bm_lon, bm_lat, aa_lon, aa_lat

Exempel
  python tools/compare_sources.py --source ghspop
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from wpc.lights import available_years, make_source
from wpc.centroids import compute_unit
from tools.change_detection import truth

# Facitår per facit; knyts till källans närmaste epok.
FACIT_YEARS = {"census": [2000, 2010, 2020], "inegi": [2010, 2020],
               "ibge": [2022], "geostat": [2011, 2021]}

_SRC = None
_KEY = None


def _get(key):
    global _SRC, _KEY
    if _KEY != key:
        if _SRC is not None:
            _SRC.close()
        name, epoch, raw_dir = key
        _SRC = make_source(name, epoch, Path(raw_dir), None, None)
        _KEY = key
    return _SRC


def _work(task):
    idx, wkb, key = task
    r = compute_unit(shapely.from_wkb(wkb), _get(key))
    return idx, r.bm_lon, r.bm_lat, r.aa_lon, r.aa_lat, r.n_lit


def nearest(epochs: list[int], year: int) -> int:
    return min(epochs, key=lambda e: abs(e - year))


def needed_units(truth_root: Path, level: int, epochs: list[int]) -> dict[int, set]:
    """epok -> facitenheter (unit_id) som behöver den epoken på denna nivå."""
    out: dict[int, set] = {e: set() for e in epochs}
    for ds, years in FACIT_YEARS.items():
        for y in years:
            e = nearest(epochs, y)
            out[e] |= set(truth(truth_root, ds, y, level, 10).index)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Centroider ur en jämförelsekälla för facitenheterna.")
    p.add_argument("--source", default="ghspop")
    p.add_argument("--raw-dir", type=Path, default=None)
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--levels", type=int, nargs="+", default=[1, 2])
    p.add_argument("--out", type=Path, default=Path("data/results/gadm41"))
    p.add_argument("--workers", type=int, default=max(1, __import__("os").cpu_count() - 1))
    a = p.parse_args(argv)
    raw_dir = a.raw_dir or Path("data/raw") / a.source
    epochs = available_years(a.source, raw_dir, None)
    print(f"{a.source}: epoker {epochs}, {a.workers} processer")

    for lv in a.levels:
        gcol = f"GID_{lv}"
        g = gpd.read_file(a.cache_dir / f"{a.boundary_set}_adm{lv}.gpkg", engine="pyogrio")
        g = g.rename(columns={gcol: "unit_id"}).set_index("unit_id")
        need = needed_units(a.truth, lv, epochs)
        for epoch in epochs:
            ids = sorted(need[epoch] & set(g.index))
            if not ids:
                continue
            out_csv = a.out / a.source / f"adm{lv}" / f"{epoch}.csv"
            if out_csv.exists():
                print(f"adm{lv} {epoch}: finns redan ({len(ids)} enheter) – hoppar över")
                continue
            geoms = g.geometry.loc[ids]
            key = (a.source, epoch, str(raw_dir))
            # stora enheter först (omslutande rektangelns yta, undviker CRS-varning)
            b = geoms.bounds
            order = np.argsort(-((b.maxx - b.minx) * (b.maxy - b.miny)).to_numpy())
            tasks = [(ids[i], shapely.to_wkb(geoms.iloc[i]), key) for i in order
                     if geoms.iloc[i] is not None]
            res = {}
            with ProcessPoolExecutor(max_workers=a.workers) as pool:
                for uid, blon, blat, alon, alat, nlit in pool.map(_work, tasks, chunksize=16):
                    res[uid] = (blon, blat, alon, alat)
            df = pd.DataFrame([{"unit_id": u, "year": epoch, "bm_lon": v[0], "bm_lat": v[1],
                                "aa_lon": v[2], "aa_lat": v[3]} for u, v in res.items()])
            out_csv.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out_csv, index=False, float_format="%.8f")
            print(f"adm{lv} {epoch}: {df.bm_lon.notna().sum()} av {len(df)} enheter med centroid -> {out_csv}")


if __name__ == "__main__":
    main()
