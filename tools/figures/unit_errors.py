"""Lägesfel per enhet (ljus mot facit, Barmore) för utgåvan v2.0.0.
Samma urval som tools/change_detection.py: complete == True och n_points >= 10.
Utdata: docs/figures/data/unit_errors.csv.gz (dataset, year, level, unit_id, gid0, err_km).
Kör från repots rot: python tools/figures/unit_errors.py"""
import numpy as np, pandas as pd
from pathlib import Path
R = Path.cwd(); REL = R / "data/release/v2.0.0"; VAL = R / "data/validation"
OUT = R / "docs/figures/data"; OUT.mkdir(parents=True, exist_ok=True)
CASES = {"census": [2000, 2010, 2020], "geostat": [2011, 2021], "inegi": [2010, 2020], "ibge": [2022]}

def gc(lon0, lat0, lon1, lat1):
    l0, p0, l1, p1 = map(np.radians, (lon0, lat0, lon1, lat1))
    a = np.sin((p1-p0)/2)**2 + np.cos(p0)*np.cos(p1)*np.sin((l1-l0)/2)**2
    return 2*6371.0*np.arcsin(np.sqrt(a))

rows = []
for lv in (0, 1, 2):
    rel = pd.read_csv(REL / f"wpc_v2_gadm41_adm{lv}.csv.gz", keep_default_na=False,
                      usecols=["unit_id", "year", "bm_lon", "bm_lat"])
    for c in ("bm_lon", "bm_lat"): rel[c] = pd.to_numeric(rel[c], errors="coerce")
    for ds, years in CASES.items():
        for y in years:
            t = pd.read_csv(VAL / ds / str(y) / f"adm{lv}.csv", keep_default_na=False)
            t = t[(t.complete.astype(str).str.lower() == "true") & (t.n_points >= 10)]
            m = t.merge(rel[rel.year == y], on="unit_id", suffixes=("_t", "_l")).dropna(subset=["bm_lon_l"])
            e = gc(m.bm_lon_t.astype(float), m.bm_lat_t.astype(float), m.bm_lon_l, m.bm_lat_l)
            rows.append(pd.DataFrame({"dataset": ds, "year": y, "level": lv, "unit_id": m.unit_id,
                                      "gid0": m.unit_id.str.split(".").str[0], "err_km": e.round(4)}))
d = pd.concat(rows); d.to_csv(OUT / "unit_errors.csv.gz", index=False)
print(d.groupby(["dataset", "level", "year"]).err_km.agg(["size", "median"]).round(3).to_string())
