#!/usr/bin/env python3
"""A3 – Detektionsmodell (docs/ANALYSIS_PLAN.md §5).

Fråga: går det att ange per enhet hur stor en förflyttning måste vara för att
synas?

Brusmått σ per enhet ur utgåvan: noise_dmsp_km för perioder inom DMSP,
√(noise_dmsp_km² + noise_viirs_km²) över skarven.

Steg (nivå 2):
  1. Spearmans ρ mellan σ och |e − t| per facit och period.
  2. Logistisk regression: rätt riktning ~ log(|t| / σ). Jämförelse: ~ log|t|.
  3. Anpassning på USA counties 2000–2020. Prövning på Mexiko 2010–2020 och USA
     2010–2020: förutsagd mot observerad andel i tio lika stora klasser,
     Brier-värde mot konstant andel, AUC för båda modellerna.
  4. Minsta upptäckbara förflyttning per enhet: det |t| där förutsagd
     sannolikhet är 0,8.

e är ändpunktsskattningen (utgåvans ljuscentroider y1→y2). Måtten (en, metrics)
och facitläsningen återanvänds ur tools/change_detection.py. Logistisk
regression och AUC görs i numpy (inget sklearn/statsmodels i miljön).

Hypoteser: H3a ρ > 0 i alla facit; H3b σ-modellen har högre AUC än |t|-modellen
i prövningen på Mexiko; H3c förutsagd och observerad andel skiljer sig med högst
0,05 i minst åtta av tio klasser i Mexiko.

Perioder (nivå 2): USA 2000–2010, 2010–2020, 2000–2020; Mexiko 2010–2020; Europa
2011–2021. Kanadas facit rörs inte.

Utdata i --out: detection_rho.csv, detection_models.csv, detection_test.csv,
detection_calibration.csv, detection_min_move.csv.

Exempel
  python tools/detection_model.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from tools.change_detection import en, truth, block_ids

SEAM = 2012
R_KM = 6371.0
RNG = np.random.default_rng(20261121)

PAIRS = {"census": [(2000, 2010), (2010, 2020), (2000, 2020)],
         "inegi": [(2010, 2020)], "geostat": [(2011, 2021)]}
FIT = ("census", 2000, 2020)                     # USA counties, anpassning
TESTS = [("inegi", 2010, 2020), ("census", 2010, 2020)]   # Mexiko och USA, prövning


# --- rena funktioner (prövas i tests/test_detection_model.py) --------------

def across_seam(y1: int, y2: int) -> bool:
    return y1 <= SEAM < y2


def unit_sigma(noise_dmsp, noise_viirs, across_seam: bool):
    """σ per enhet: DMSP-brus inom perioden, kvadratsumman över skarven."""
    nd = np.asarray(noise_dmsp, float)
    nv = np.asarray(noise_viirs, float)
    return np.sqrt(nd ** 2 + nv ** 2) if across_seam else nd


def logit(p):
    return np.log(p / (1 - p))


def logistic_fit(X, y, max_iter: int = 100, tol: float = 1e-10) -> np.ndarray:
    """Logistisk regression med IRLS/Newton. X utan interceptkolumn; den läggs
    till. Returnerar koefficienterna (intercept först)."""
    X = np.column_stack([np.ones(len(y)), np.asarray(X, float).reshape(len(y), -1)])
    y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    for _ in range(max_iter):
        p = 1.0 / (1.0 + np.exp(-(X @ beta)))
        w = np.clip(p * (1 - p), 1e-12, None)
        grad = X.T @ (y - p)
        hess = (X * w[:, None]).T @ X
        step = np.linalg.solve(hess, grad)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta


def predict_prob(beta: np.ndarray, X) -> np.ndarray:
    X = np.column_stack([np.ones(np.asarray(X).reshape(-1, beta.size - 1).shape[0]),
                         np.asarray(X, float).reshape(-1, beta.size - 1)])
    return 1.0 / (1.0 + np.exp(-(X @ beta)))


def auc(scores, y) -> float:
    """AUC via Mann–Whitney-U med medelrang för lika värden."""
    y = np.asarray(y).astype(int)
    s = np.asarray(scores, float)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def brier(p, y) -> float:
    return float(np.mean((np.asarray(p, float) - np.asarray(y, float)) ** 2))


def calibration_bins(p, y, nbins: int = 10) -> pd.DataFrame:
    """Tio (nbins) lika stora klasser efter förutsagd sannolikhet; förutsagt
    medel mot observerad andel per klass."""
    p = np.asarray(p, float)
    y = np.asarray(y, float)
    order = np.argsort(p, kind="mergesort")
    rows = []
    for b, idx in enumerate(np.array_split(order, nbins)):
        rows.append({"bin": b + 1, "n": len(idx),
                     "pred_mean": float(p[idx].mean()), "obs_prop": float(y[idx].mean())})
    return pd.DataFrame(rows)


def min_detectable_move(sigma, beta0: float, beta1: float, target: float = 0.8):
    """|t| där modellen rätt riktning ~ β0 + β1·log(|t|/σ) ger sannolikheten
    target. β1 > 0 krävs för att ett sådant |t| ska finnas."""
    return np.asarray(sigma, float) * np.exp((logit(target) - beta0) / beta1)


# --- per facit och period: enhetstabell på nivå 2 --------------------------

def load_sigma(release: Path) -> pd.DataFrame:
    d = pd.read_csv(release / "wpc_v2_gadm41_adm2.csv.gz", keep_default_na=False,
                    low_memory=False, usecols=["unit_id", "noise_dmsp_km", "noise_viirs_km"])
    for c in ("noise_dmsp_km", "noise_viirs_km"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna().drop_duplicates("unit_id").set_index("unit_id")


def load_light(release: Path) -> pd.DataFrame:
    d = pd.read_csv(release / "wpc_v2_gadm41_adm2.csv.gz", keep_default_na=False,
                    low_memory=False, usecols=["unit_id", "year", "bm_lon", "bm_lat"])
    for c in ("bm_lon", "bm_lat"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna().set_index(["unit_id", "year"])


def unit_table(truth_root, ds, y1, y2, sigma, light, min_points) -> pd.DataFrame:
    """Per nivå 2-enhet: |t|, |e−t|, rätt riktning och σ."""
    tru = {y: truth(truth_root, ds, y, 2, min_points) for y in (y1, y2)}
    ids = tru[y1].index.intersection(tru[y2].index)
    k1 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y1)])
    k2 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y2)])
    L1, L2 = light.reindex(k1), light.reindex(k2)
    sig = sigma.reindex(ids)
    ok = L1.bm_lon.notna().values & L2.bm_lon.notna().values & sig.noise_dmsp_km.notna().values
    ids = ids[ok]
    T1, T2, L1, L2, sig = tru[y1].loc[ids], tru[y2].loc[ids], L1[ok], L2[ok], sig[ok]
    t = np.column_stack(en(T1.bm_lon, T1.bm_lat, T2.bm_lon, T2.bm_lat))
    e = np.column_stack(en(L1.bm_lon.values, L1.bm_lat.values, L2.bm_lon.values, L2.bm_lat.values))
    tn = np.hypot(*t.T)
    err = np.hypot(*(e - t).T)
    cos = (t * e).sum(1) / np.maximum(tn * np.hypot(*e.T), 1e-12)
    s = unit_sigma(sig.noise_dmsp_km.values, sig.noise_viirs_km.values, across_seam(y1, y2))
    return pd.DataFrame({"unit_id": ids, "t_abs": tn, "err": err,
                         "dir_right": (cos > 0).astype(int), "sigma": s}).set_index("unit_id")


def boot_rho_ci(sigma_v, err_v, blocks, n_boot) -> tuple[float, float]:
    ub, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(len(ub))]
    vals = []
    for _ in range(n_boot):
        idx = np.concatenate([members[k] for k in RNG.integers(0, len(ub), len(ub))])
        vals.append(pd.Series(sigma_v[idx]).corr(pd.Series(err_v[idx]), method="spearman"))
    return float(np.nanquantile(vals, 0.025)), float(np.nanquantile(vals, 0.975))


def model_frame(tab: pd.DataFrame) -> pd.DataFrame:
    """Enheter med |t| > 0 och σ > 0, med prediktorerna log|t| och log(|t|/σ)."""
    m = tab[(tab.t_abs > 0) & (tab.sigma > 0)].copy()
    m["log_t"] = np.log(m.t_abs)
    m["log_snr"] = np.log(m.t_abs / m.sigma)
    return m


def main(argv=None):
    p = argparse.ArgumentParser(description="A3 – detektionsmodell.")
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--min-points", type=int, default=10)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    sigma, light = load_sigma(a.release), load_light(a.release)
    tabs = {(ds, y1, y2): unit_table(a.truth, ds, y1, y2, sigma, light, a.min_points)
            for ds, prs in PAIRS.items() for (y1, y2) in prs}

    # Steg 1: Spearmans ρ(σ, |e−t|)
    rho_rows = []
    for (ds, y1, y2), tab in tabs.items():
        rho = pd.Series(tab.sigma.values).corr(pd.Series(tab.err.values), method="spearman")
        multi = tab.index.str[:3].nunique() >= 5
        lo, hi = boot_rho_ci(tab.sigma.values, tab.err.values, block_ids(tab.index, 2, multi), a.n_boot)
        rho_rows.append({"dataset": ds, "from": y1, "to": y2, "level": 2, "n": len(tab),
                         "rho": rho, "rho_lo": lo, "rho_hi": hi})
    pd.DataFrame(rho_rows).to_csv(a.out / "detection_rho.csv", index=False, float_format="%.4f")

    # Steg 2–3: anpassning på USA counties 2000–2020, båda modellerna
    fit = model_frame(tabs[FIT])
    b_snr = logistic_fit(fit[["log_snr"]].values, fit.dir_right.values)
    b_t = logistic_fit(fit[["log_t"]].values, fit.dir_right.values)
    pd.DataFrame([
        {"model": "snr", "predictor": "log(|t|/sigma)", "fit_dataset": FIT[0], "from": FIT[1], "to": FIT[2],
         "n": len(fit), "beta0": b_snr[0], "beta1": b_snr[1]},
        {"model": "abs_t", "predictor": "log|t|", "fit_dataset": FIT[0], "from": FIT[1], "to": FIT[2],
         "n": len(fit), "beta0": b_t[0], "beta1": b_t[1]},
    ]).to_csv(a.out / "detection_models.csv", index=False, float_format="%.6f")

    # Prövning: AUC, Brier mot konstant, tio-klasskalibrering (σ-modellen)
    test_rows, cal_rows = [], []
    for ds, y1, y2 in TESTS:
        te = model_frame(tabs[(ds, y1, y2)])
        y = te.dir_right.values
        p_snr = predict_prob(b_snr, te[["log_snr"]].values)
        p_t = predict_prob(b_t, te[["log_t"]].values)
        base = y.mean()
        for name, pr in (("snr", p_snr), ("abs_t", p_t)):
            test_rows.append({"test_dataset": ds, "from": y1, "to": y2, "model": name, "n": len(te),
                              "base_rate": base, "auc": auc(pr, y), "brier": brier(pr, y),
                              "brier_const": brier(np.full(len(y), base), y)})
        cal = calibration_bins(p_snr, y, nbins=10)
        cal.insert(0, "to", y2); cal.insert(0, "from", y1); cal.insert(0, "test_dataset", ds)
        cal["abs_gap"] = (cal.pred_mean - cal.obs_prop).abs()
        cal_rows.append(cal)
    pd.DataFrame(test_rows).to_csv(a.out / "detection_test.csv", index=False, float_format="%.4f")
    pd.concat(cal_rows, ignore_index=True).to_csv(a.out / "detection_calibration.csv",
                                                  index=False, float_format="%.4f")

    # Steg 4: minsta upptäckbara förflyttning per enhet (σ-modellen anpassad på USA 2000–2020)
    mm_rows = []
    for (ds, y1, y2), tab in tabs.items():
        m = model_frame(tab)
        mdm = min_detectable_move(m.sigma.values, b_snr[0], b_snr[1], target=0.8)
        mm_rows.append({"dataset": ds, "from": y1, "to": y2, "n": len(m),
                        "sigma_med": float(np.median(m.sigma.values)),
                        "min_move_q25": float(np.quantile(mdm, 0.25)),
                        "min_move_med": float(np.median(mdm)),
                        "min_move_q75": float(np.quantile(mdm, 0.75))})
    pd.DataFrame(mm_rows).to_csv(a.out / "detection_min_move.csv", index=False, float_format="%.4f")

    print("ρ(σ, |e−t|) nivå 2:")
    print(pd.DataFrame(rho_rows)[["dataset", "from", "to", "n", "rho", "rho_lo", "rho_hi"]].round(3).to_string(index=False))
    print("\nmodeller (USA counties 2000–2020):")
    print(f"  snr:   beta0={b_snr[0]:.3f} beta1={b_snr[1]:.3f}")
    print(f"  abs_t: beta0={b_t[0]:.3f} beta1={b_t[1]:.3f}")
    print("\nprövning:")
    print(pd.DataFrame(test_rows)[["test_dataset", "from", "to", "model", "n", "base_rate", "auc", "brier", "brier_const"]].round(3).to_string(index=False))
    print("\nminsta upptäckbara förflyttning (km):")
    print(pd.DataFrame(mm_rows)[["dataset", "from", "to", "n", "sigma_med", "min_move_med", "min_move_q25", "min_move_q75"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
