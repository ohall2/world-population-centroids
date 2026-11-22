"""Förflyttning per enhet (nivå 2), sann (facit) och mätt (ljus, utgåvan v2.0.0), i km öst/nord.
Samma urval som tools/change_detection.py (complete, n_points >= 10 båda åren).
Kanada (statcan) och Japan (estat_jp) läses ur data/validation_canada och
data/validation_japan och begränsas till landet, som i tools/canada_confirm.py.
Utdata: docs/figures/data/unit_shifts.csv.gz. Kör från repots rot."""
import numpy as np, pandas as pd
from pathlib import Path
R = Path.cwd(); REL = R / "data/release/v2.0.0"; VAL = R / "data/validation"
OUT = R / "docs/figures/data"; OUT.mkdir(parents=True, exist_ok=True)
PAIRS = [("census", 2000, 2020), ("census", 2000, 2010), ("census", 2010, 2020),
         ("inegi", 2010, 2020), ("geostat", 2011, 2021), ("statcan", 2016, 2021), ("estat_jp", 2015, 2020)]
ROOT = {"statcan": (R / "data/validation_canada", "CAN"), "estat_jp": (R / "data/validation_japan", "JPN")}
def en(lon0, lat0, lon1, lat1):
    dlon = (np.asarray(lon1) - np.asarray(lon0) + 180.0) % 360.0 - 180.0
    return 6371.0*np.radians(dlon)*np.cos(np.radians(lat0)), 6371.0*np.radians(np.asarray(lat1) - np.asarray(lat0))
def truth(ds, y):
    root, iso = ROOT.get(ds, (VAL, None))
    t = pd.read_csv(root / ds / str(y) / "adm2.csv", keep_default_na=False)
    t = t[(t.complete.astype(str).str.lower() == "true") & (t.n_points >= 10)]
    if iso is not None: t = t[t.unit_id.str[:3] == iso]
    return t.set_index("unit_id")[["bm_lon", "bm_lat"]].astype(float)
rel = pd.read_csv(REL / "wpc_v2_gadm41_adm2.csv.gz", keep_default_na=False, usecols=["unit_id", "year", "bm_lon", "bm_lat"])
for c in ("bm_lon", "bm_lat"): rel[c] = pd.to_numeric(rel[c], errors="coerce")
rows = []
for ds, a, b in PAIRS:
    t0, t1 = truth(ds, a), truth(ds, b)
    l0 = rel[rel.year == a].set_index("unit_id")[["bm_lon", "bm_lat"]]
    l1 = rel[rel.year == b].set_index("unit_id")[["bm_lon", "bm_lat"]]
    ids = t0.index.intersection(t1.index).intersection(l0.dropna().index).intersection(l1.dropna().index)
    te, tn = en(t0.loc[ids, "bm_lon"], t0.loc[ids, "bm_lat"], t1.loc[ids, "bm_lon"], t1.loc[ids, "bm_lat"])
    me, mn = en(l0.loc[ids, "bm_lon"], l0.loc[ids, "bm_lat"], l1.loc[ids, "bm_lon"], l1.loc[ids, "bm_lat"])
    rows.append(pd.DataFrame({"dataset": ds, "from": a, "to": b, "unit_id": ids, "true_e": te, "true_n": tn,
                              "meas_e": me, "meas_n": mn}).round(4))
d = pd.concat(rows); d.to_csv(OUT / "unit_shifts.csv.gz", index=False)
tm = np.hypot(d.true_e, d.true_n); err = np.hypot(d.meas_e - d.true_e, d.meas_n - d.true_n)
cosang = (d.true_e*d.meas_e + d.true_n*d.meas_n) > 0
s = d.assign(true=tm, err=err, right=cosang).groupby(["dataset", "from", "to"]).agg(
    n=("true", "size"), true_med=("true", "median"), err_med=("err", "median"), dir_right=("right", "mean"))
print(s.round(4).to_string())
