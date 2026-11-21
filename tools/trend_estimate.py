#!/usr/bin/env python3
"""A2 – Trendskattning (docs/ANALYSIS_PLAN.md §4).

Fråga: minskar felet om hela tidsserien används i stället för två ändpunkter?

Per enhet och koordinat (km öster och norr i ett lokalt plan kring enhetens
ljuscentroid 2013) skattas förflyttningen ur en Theil–Sen-lutning inom varje
sensorperiod (DMSP t.o.m. 2012, VIIRS fr.o.m. 2013):
  - Period helt inom DMSP (y1, y2 ≤ 2012): lutning över åren y1−2 … y2+2, högst
    till 2012; flytt = lutning × (y2 − y1).
  - Period över skarven: s_D över y1−2 … 2012 och s_V över 2013 … y2+2, högst
    till 2025; flytt = s_D × (2013 − y1) + s_V × (y2 − 2013). Skarvens
    förskjutning ingår inte.
  - Minst fyra år krävs i varje del, annars utgår enheten.

Jämförs med ändpunkter (nuvarande, ettårsfönster) och treårsmedel (3D-medel av
tre år), på samma enheter. Ändpunkts- och treårsskattningarna och måtten
återanvänds oförändrade ur tools/change_detection.py och tools/robustness_checks.py.

Mått: median |e − t| och rätt riktning, nivå 1 och 2, med blockbootstrap.

Hypotes H2: trendskattningen minskar median |e − t| på nivå 2 med minst 20 % mot
ändpunkter i minst fyra av de fem perioderna (treårsmedel gav 6–15 %).

Perioder (planen §2): USA 2000–2010, 2010–2020, 2000–2020; Mexiko 2010–2020;
Europa 2011–2021. Kanadas facit rörs inte.

Utdata i --out: <ds>_trend_change.csv per facit.

Exempel
  python tools/trend_estimate.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tools.change_detection import en, metrics, block_ids, bootstrap, truth
from tools.robustness_checks import load_release, light_pos

SEAM = 2012          # sista DMSP-året; VIIRS fr.o.m. 2013
DMSP_MIN = 1992      # första året i serien
VIIRS_MAX = 2025     # sista året i serien
REF_YEAR = 2013      # lokalt plan kring enhetens centroid detta år
R_KM = 6371.0

PAIRS = {"census": [(2000, 2010), (2010, 2020), (2000, 2020)],
         "inegi": [(2010, 2020)], "geostat": [(2011, 2021)]}
SERIES = ["endpoint", "three_year", "trend"]


# --- rena funktioner (prövas i tests/test_trend_estimate.py) ---------------

def theil_sen_slope(years, values) -> float:
    """Theil–Sen-lutning: medianen av alla parvisa lutningar."""
    years = np.asarray(years, float)
    values = np.asarray(values, float)
    i, j = np.triu_indices(len(years), k=1)
    dt = years[j] - years[i]
    ok = dt != 0
    return float(np.median((values[j] - values[i])[ok] / dt[ok]))


def trend_spec(y1: int, y2: int) -> list[tuple[int, int, int]]:
    """Delar att skatta lutning i, som (år_lo, år_hi, faktor).
    Flytt = Σ lutning_del × faktor."""
    if y2 <= SEAM:
        return [(max(y1 - 2, DMSP_MIN), min(y2 + 2, SEAM), y2 - y1)]
    return [(max(y1 - 2, DMSP_MIN), SEAM, REF_YEAR - y1),
            (REF_YEAR, min(y2 + 2, VIIRS_MAX), y2 - REF_YEAR)]


def trend_movement_series(years, east, north, y1: int, y2: int) -> tuple[float, float]:
    """Flytt (km öster, norr) ur Theil–Sen-lutningar per sensordel. NaN om någon
    del har färre än fyra år med centroid."""
    years = np.asarray(years)
    east = np.asarray(east, float)
    north = np.asarray(north, float)
    me = mn = 0.0
    for lo, hi, factor in trend_spec(y1, y2):
        m = (years >= lo) & (years <= hi) & np.isfinite(east) & np.isfinite(north)
        if m.sum() < 4:
            return np.nan, np.nan
        me += theil_sen_slope(years[m], east[m]) * factor
        mn += theil_sen_slope(years[m], north[m]) * factor
    return me, mn


# --- trendskattning över ett enhetsurval -----------------------------------

def trend_e(s: pd.DataFrame, ids: pd.Index, y1: int, y2: int) -> np.ndarray:
    """(n, 2)-matris med trendflytt per enhet; NaN-rad om enheten utgår.

    Alla år i periodens delar läggs i ett lokalt plan kring enhetens
    ljuscentroid REF_YEAR (en() med referensen som origo)."""
    spec = trend_spec(y1, y2)
    yrs = np.array(sorted({y for lo, hi, _ in spec for y in range(lo, hi + 1)}))
    mi = pd.MultiIndex.from_product([ids, yrs])
    r = s.reindex(mi)
    lon = r.bm_lon.to_numpy().reshape(len(ids), len(yrs))
    lat = r.bm_lat.to_numpy().reshape(len(ids), len(yrs))
    ref = s.reindex(pd.MultiIndex.from_arrays([ids, np.full(len(ids), REF_YEAR)]))
    rlon = ref.bm_lon.to_numpy()[:, None]
    rlat = ref.bm_lat.to_numpy()[:, None]
    east, north = en(rlon, rlat, lon, lat)    # lokalt plan kring 2013-centroiden
    out = np.full((len(ids), 2), np.nan)
    for k in range(len(ids)):
        out[k] = trend_movement_series(yrs, east[k], north[k], y1, y2)
    return out


def estimators(s, tru, y1, y2, keep):
    """ids och t (facit) samt e för varje skattning, på de enheter där alla tre
    skattningarna och facit finns (planen: 'på samma enheter')."""
    ids = tru[y1].index.intersection(tru[y2].index)
    if keep is not None:
        ids = ids[keep(ids)]
    ep1, ep2 = light_pos(s, ids, [y1]), light_pos(s, ids, [y2])
    t31 = light_pos(s, ids, [y1 - 1, y1, y1 + 1])
    t32 = light_pos(s, ids, [y2 - 1, y2, y2 + 1])
    e = {"endpoint": np.column_stack(en(ep1.lon.values, ep1.lat.values, ep2.lon.values, ep2.lat.values)),
         "three_year": np.column_stack(en(t31.lon.values, t31.lat.values, t32.lon.values, t32.lat.values)),
         "trend": trend_e(s, ids, y1, y2)}
    T1, T2 = tru[y1].loc[ids], tru[y2].loc[ids]
    t = np.column_stack(en(T1.bm_lon, T1.bm_lat, T2.bm_lon, T2.bm_lat))
    ok = np.isfinite(t).all(1)
    for v in e.values():
        ok &= np.isfinite(v).all(1)
    ids = ids[ok]
    return ids, t[ok], {k: v[ok] for k, v in e.items()}


def run_dataset(ds, truth_root, rel, min_points, n_boot):
    rows = []
    for lv in (1, 2):
        keep = (lambda i: i.str.startswith("USA.")) if ds == "census" and lv == 1 else None
        ys = {y for pr in PAIRS[ds] for y in pr}
        tru = {y: truth(truth_root, ds, y, lv, min_points) for y in ys}
        for y1, y2 in PAIRS[ds]:
            ids, t, e = estimators(rel[lv], tru, y1, y2, keep)
            if len(ids) == 0:
                continue
            multi = ids.str[:3].nunique() >= 5
            blk = block_ids(ids, lv, multi)
            nb = len(np.unique(blk))
            for name in SERIES:
                m = metrics(t, e[name])
                ci = bootstrap(t, e[name], blk, n_boot)
                rows.append({"dataset": ds, "series": name, "level": lv, "from": y1, "to": y2,
                             "n": len(ids), "n_blocks": nb, **m, **ci})
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description="A2 – trendskattning av förflyttning.")
    p.add_argument("--datasets", nargs="+", default=list(PAIRS), choices=list(PAIRS))
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--min-points", type=int, default=10)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)

    a.out.mkdir(parents=True, exist_ok=True)
    rel = {lv: load_release(a.release, lv) for lv in (1, 2)}
    for ds in a.datasets:
        print(f"== {ds} ==")
        res = run_dataset(ds, a.truth, rel, a.min_points, a.n_boot)
        res.to_csv(a.out / f"{ds}_trend_change.csv", index=False, float_format="%.4f")
        cols = ["series", "level", "from", "to", "n", "err_med", "err_med_lo", "err_med_hi",
                "dir_right", "dir_right_lo", "dir_right_hi"]
        print(res[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
