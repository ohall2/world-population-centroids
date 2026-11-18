#!/usr/bin/env python3
"""
World Population Centroids v2 – kan nattljusets centroider mäta hur
befolkningen förflyttar sig? Statistisk prövning mot folkräkningar.

För varje enhet i, och period y1 → y2:
  t_i  verklig förflyttning: facitcentroid y2 − facitcentroid y1 (km öster, norr)
  e_i  uppmätt förflyttning: ljuscentroid y2 − ljuscentroid y1
Mått (alla Barmore):
  n                   antal enheter (fullständiga i facit båda åren, minst
                      --min-points facitpunkter, centroid i serien båda åren)
  true_med, err_med   median |t| och median |e − t| (km)
  ratio               err_med / true_med (samma mått som validate_series.py)
  r                   korrelation mellan e och t, komponenterna öst och nord
                      sammanslagna (2n värden)
  rho                 Spearmans rangkorrelation, samma värden (robust mot
                      enstaka stora förflyttningar)
  slope               OLS-lutning för e på t (med intercept), sammanslagna
                      komponenter; 1 = rätt storlek, 0 = inget samband
  dir_right           andel enheter där e pekar åt rätt håll (cos(e, t) > 0);
                      0,5 = slump
  snr_gt1             andel enheter där |t| > |e − t|
För alla mått ges 95 %-intervall med blockbootstrap (1 000 dragningar).
Block: för nivå 2 överordnad nivå 1-enhet, för nivå 1 landet om facit har
minst fem länder (Europa), annars enheten själv.

Därtill ett mått på hur stor en förflyttning måste vara för att synas:
enheterna delas i klasser efter |t|, och för varje klass redovisas andelen
med rätt riktning och median |e − t|.

Utdata: --out/<dataset>_change_detection.csv och _change_detection_bins.csv.

Exempel
  python tools/change_detection.py --dataset census --pairs 2000-2010 2010-2020 2000-2020 \\
      --release data/release/v2.0.0 --out docs/change_detection
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R_KM = 6371.0
RNG = np.random.default_rng(20261118)


def truth(root: Path, ds: str, year: int, level: int, min_points: int) -> pd.DataFrame:
    d = pd.read_csv(root / ds / str(year) / f"adm{level}.csv", keep_default_na=False)
    d = d[(d.complete.astype(str).str.lower() == "true") & (d.n_points >= min_points)]
    return d.set_index("unit_id")[["bm_lon", "bm_lat", "pop"]].astype(float)


def release(path: Path, level: int, years: set[int]) -> pd.DataFrame:
    d = pd.read_csv(path / f"wpc_v2_gadm41_adm{level}.csv.gz", keep_default_na=False,
                    usecols=["unit_id", "year", "bm_lon", "bm_lat"])
    d = d[d.year.isin(years)]
    for c in ("bm_lon", "bm_lat"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def results(path: Path, level: int, years: set[int]) -> pd.DataFrame:
    parts = []
    for y in years:
        d = pd.read_csv(path / f"adm{level}" / f"{y}.csv", keep_default_na=False,
                        usecols=["unit_id", "year", "bm_lon", "bm_lat"])
        for c in ("bm_lon", "bm_lat"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
        parts.append(d)
    return pd.concat(parts)


def en(lon0, lat0, lon1, lat1):
    dlon = (np.asarray(lon1) - np.asarray(lon0) + 180.0) % 360.0 - 180.0
    return (R_KM * np.radians(dlon) * np.cos(np.radians(lat0)),
            R_KM * np.radians(np.asarray(lat1) - np.asarray(lat0)))


def metrics(t: np.ndarray, e: np.ndarray) -> dict:
    """t, e: (n, 2) i km."""
    tn = np.hypot(*t.T)
    err = np.hypot(*(e - t).T)
    tt, ee = t.ravel(), e.ravel()
    r = np.corrcoef(tt, ee)[0, 1] if tt.std() > 0 and ee.std() > 0 else np.nan
    slope = np.polyfit(tt, ee, 1)[0] if tt.std() > 0 else np.nan
    rho = pd.Series(tt).corr(pd.Series(ee), method="spearman")
    cos = (t * e).sum(1) / np.maximum(tn * np.hypot(*e.T), 1e-12)
    ok = tn > 0
    return {"true_med": np.median(tn), "err_med": np.median(err), "ratio": np.median(err) / np.median(tn),
            "r": r, "rho": rho, "slope": slope, "dir_right": np.mean(cos[ok] > 0), "snr_gt1": np.mean(tn > err)}


def block_ids(units: pd.Index, level: int, multi_country: bool) -> np.ndarray:
    u = units.to_series()
    if level == 2:
        # USA.5.13_1 -> USA.5
        return u.str.extract(r"^([A-Z]{3}\.\d+)")[0].fillna(u).values
    if level == 1 and multi_country:
        return u.str[:3].values
    return u.values


def bootstrap(t, e, blocks, n_boot=1000) -> dict:
    ub, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(len(ub))]
    out = []
    for _ in range(n_boot):
        pick = RNG.integers(0, len(ub), len(ub))
        idx = np.concatenate([members[k] for k in pick])
        out.append(metrics(t[idx], e[idx]))
    df = pd.DataFrame(out)
    return {f"{k}_lo": df[k].quantile(0.025) for k in df} | {f"{k}_hi": df[k].quantile(0.975) for k in df}


def main(argv=None):
    p = argparse.ArgumentParser(description="Prövning av om ljusets centroider mäter förflyttning.")
    p.add_argument("--dataset", choices=["census", "geostat", "inegi", "ibge"], required=True)
    p.add_argument("--pairs", nargs="+", required=True, help="Perioder som 2000-2010.")
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--series", nargs="*", default=[],
                   help="Fler serier som namn=mapp med adm<N>/<år>.csv, t.ex. "
                        "chiov_bltc=data/results_mask2km/gadm41/chiovelli-bltc")
    p.add_argument("--levels", type=int, nargs="+", default=[1, 2])
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--min-points", type=int, default=10)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/change_detection"))
    a = p.parse_args(argv)

    pairs = [tuple(int(x) for x in s.split("-")) for s in a.pairs]
    years = {y for pr in pairs for y in pr}
    series = [("release", None)] + [tuple(s.split("=", 1)) for s in a.series]
    rows, bins = [], []
    for lv in a.levels:
        tru = {y: truth(a.truth, a.dataset, y, lv, a.min_points) for y in years}
        for name, path in series:
            s = release(a.release, lv, years) if path is None else results(Path(path), lv, years)
            s = s.set_index(["unit_id", "year"])
            for y1, y2 in pairs:
                ids = tru[y1].index.intersection(tru[y2].index)
                k1 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y1)])
                k2 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y2)])
                s1, s2 = s.reindex(k1), s.reindex(k2)
                ok = (s1.bm_lon.notna().values & s2.bm_lon.notna().values)
                ids = ids[ok]
                T1, T2 = tru[y1].loc[ids], tru[y2].loc[ids]
                S1, S2 = s1[ok], s2[ok]
                t = np.column_stack(en(T1.bm_lon, T1.bm_lat, T2.bm_lon, T2.bm_lat))
                e = np.column_stack(en(S1.bm_lon.values, S1.bm_lat.values, S2.bm_lon.values, S2.bm_lat.values))
                multi = ids.str[:3].nunique() >= 5   # Europa; USA + Puerto Rico räknas som ett
                m = metrics(t, e)
                ci = bootstrap(t, e, block_ids(ids, lv, multi), a.n_boot)
                nb = len(np.unique(block_ids(ids, lv, multi)))
                rows.append({"dataset": a.dataset, "series": name, "level": lv, "from": y1, "to": y2,
                             "n": len(ids), "n_blocks": nb, **m, **ci})
                tn = np.hypot(*t.T)
                err = np.hypot(*(e - t).T)
                cos = (t * e).sum(1) / np.maximum(tn * np.hypot(*e.T), 1e-12)
                edges = [0, 0.1, 0.25, 0.5, 1, 2, 5, 10, np.inf]
                cl = pd.cut(tn, edges, right=False)
                for c, g in pd.DataFrame({"cl": cl, "ok": cos > 0, "err": err, "tn": tn}).groupby("cl", observed=True):
                    bins.append({"dataset": a.dataset, "series": name, "level": lv, "from": y1, "to": y2,
                                 "true_move_km": str(c), "n": len(g), "dir_right": g.ok.mean(),
                                 "true_med": g.tn.median(), "err_med": g.err.median()})
    a.out.mkdir(parents=True, exist_ok=True)
    res = pd.DataFrame(rows)
    res.to_csv(a.out / f"{a.dataset}_change_detection.csv", index=False, float_format="%.4f")
    pd.DataFrame(bins).to_csv(a.out / f"{a.dataset}_change_detection_bins.csv", index=False, float_format="%.4f")
    show = res[["series", "level", "from", "to", "n", "true_med", "err_med", "ratio", "rho", "rho_lo", "rho_hi",
                "slope", "slope_lo", "slope_hi", "dir_right", "dir_right_lo", "dir_right_hi", "snr_gt1"]]
    print(show.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
