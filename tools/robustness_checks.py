#!/usr/bin/env python3
"""
Fyra kontroller av testet av förändring och av landsfelet, begärda i
granskningen av manuset (docs/review/review-manuscript-2026-11-18.md).
Bara befintliga data: data/validation (facit) och data/release/v2.0.0.

1. USA:s delstater (nivå 1, utan Puerto Rico): förändring 2000–2010,
   2010–2020, 2000–2020, med samma mått som tools/change_detection.py.
2. Treårsmedel: ljusets centroid som medel av tre år kring varje
   folkräkningsår (3D-medel av Barmore-centroiderna), jämfört med ett år.
3. Permutationstest: den uppmätta förflyttningen blandas mellan enheter
   inom samma block (nivå 2: överordnad region; nivå 1: land eller inget),
   1 000 gånger. Ger nollfördelningen för rätt riktning och rangkorrelation
   i stället för antagandet att slumpen ger 0,5 och 0.
4. Korrektion av landsfelet: landets centroid som befolkningsviktat medel
   av regionernas ljuscentroider (vikt = regionens befolkning i facit), mot
   den ljusviktade nationella centroiden i utgåvan. Följer av identiteten:
   befolkningsvikter mellan regioner tar bort mellanregionsdelen Bb.

Utdata i --out: robust_states.csv, robust_3yr.csv, robust_permutation.csv,
robust_national_correction.csv.

Exempel (från repots rot)
  python tools/robustness_checks.py --out docs/robustness
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from change_detection import block_ids, bootstrap, en, metrics, truth  # noqa: E402

RNG = np.random.default_rng(20261118)
R_KM = 6371.0
PAIRS = {"census": [(2000, 2010), (2010, 2020), (2000, 2020)], "inegi": [(2010, 2020)], "geostat": [(2011, 2021)]}


def gc(lon0, lat0, lon1, lat1):
    l0, p0, l1, p1 = map(np.radians, (lon0, lat0, lon1, lat1))
    a = np.sin((p1 - p0) / 2) ** 2 + np.cos(p0) * np.cos(p1) * np.sin((l1 - l0) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(a))


def vec(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def unvec(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arctan2(v[..., 1], v[..., 0])), np.degrees(np.arcsin(v[..., 2]))


def load_release(rel: Path, level: int) -> pd.DataFrame:
    d = pd.read_csv(rel / f"wpc_v2_gadm41_adm{level}.csv.gz", keep_default_na=False,
                    usecols=["unit_id", "year", "bm_lon", "bm_lat"])
    for c in ("bm_lon", "bm_lat"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna().set_index(["unit_id", "year"])


def light_pos(s: pd.DataFrame, ids, years) -> pd.DataFrame:
    """Medel (3D) av ljuscentroiden över years; NaN om något år saknas."""
    vs = []
    for y in years:
        k = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y)])
        r = s.reindex(k)
        vs.append(vec(r.bm_lon.values, r.bm_lat.values))
    v = np.nanmean(np.stack(vs), axis=0) if len(vs) > 1 else vs[0]
    bad = np.any([np.isnan(x).any(1) for x in vs], axis=0)
    lon, lat = unvec(v)
    lon[bad] = np.nan
    return pd.DataFrame({"lon": lon, "lat": lat}, index=ids)


def change(tru, s, ds, lv, y1, y2, window, keep=None):
    ids = tru[y1].index.intersection(tru[y2].index)
    if keep is not None:
        ids = ids[keep(ids)]
    w = window // 2
    L1 = light_pos(s, ids, range(y1 - w, y1 + w + 1))
    L2 = light_pos(s, ids, range(y2 - w, y2 + w + 1))
    ok = L1.lon.notna().values & L2.lon.notna().values
    ids = ids[ok]
    T1, T2, L1, L2 = tru[y1].loc[ids], tru[y2].loc[ids], L1[ok], L2[ok]
    t = np.column_stack(en(T1.bm_lon, T1.bm_lat, T2.bm_lon, T2.bm_lat))
    e = np.column_stack(en(L1.lon.values, L1.lat.values, L2.lon.values, L2.lat.values))
    multi = ids.str[:3].nunique() >= 5
    return ids, t, e, block_ids(ids, lv, multi)


def row(ds, lv, y1, y2, ids, t, e, blocks, **extra):
    return {"dataset": ds, "level": lv, "from": y1, "to": y2, "n": len(ids), **extra,
            **metrics(t, e), **bootstrap(t, e, blocks)}


def permutation(t, e, blocks, n=1000):
    obs = metrics(t, e)
    ub, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(len(ub))]
    dr, rho = [], []
    for _ in range(n):
        idx = np.empty(len(e), int)
        for m in members:
            idx[m] = RNG.permutation(m)
        mm = metrics(t, e[idx])
        dr.append(mm["dir_right"])
        rho.append(mm["rho"])
    dr, rho = np.array(dr), np.array(rho)
    return {"dir_right": obs["dir_right"], "dir_right_null_mean": dr.mean(),
            "dir_right_null_lo": np.quantile(dr, 0.025), "dir_right_null_hi": np.quantile(dr, 0.975),
            "p_dir_right_upper": (np.sum(dr >= obs["dir_right"]) + 1) / (n + 1),
            "p_dir_right_lower": (np.sum(dr <= obs["dir_right"]) + 1) / (n + 1),
            "rho": obs["rho"], "rho_null_lo": np.quantile(rho, 0.025), "rho_null_hi": np.quantile(rho, 0.975),
            "p_rho_upper": (np.sum(rho >= obs["rho"]) + 1) / (n + 1),
            "p_rho_lower": (np.sum(rho <= obs["rho"]) + 1) / (n + 1)}


def national_correction(val: Path, rel: Path, ds: str, years) -> list[dict]:
    """Befolkningsviktat medel av regionernas ljuscentroider mot facit och utgåvans landscentroid."""
    s1, s0 = load_release(rel, 1), load_release(rel, 0)
    out = []
    for y in years:
        t0 = pd.read_csv(val / ds / str(y) / "adm0.csv", keep_default_na=False)
        t0 = t0[t0.complete.astype(str).str.lower() == "true"].set_index("unit_id")
        t1 = pd.read_csv(val / ds / str(y) / "adm1.csv", keep_default_na=False)
        t1["pop"] = pd.to_numeric(t1["pop"], errors="coerce")
        t1 = t1[t1["pop"] > 0].set_index("unit_id")
        for c in t0.index:
            reg = t1[t1.index.str.split(".").str[0] == c]
            k = pd.MultiIndex.from_arrays([reg.index, np.full(len(reg), y)])
            L = s1.reindex(k)
            ok = L.bm_lon.notna().values
            if ok.sum() < 2:
                continue
            v = (vec(L.bm_lon.values[ok], L.bm_lat.values[ok]) * reg["pop"].values[ok, None]).sum(0)
            plon, plat = unvec(v)
            share = reg["pop"].values[ok].sum() / reg["pop"].sum()
            tl, tla = float(t0.loc[c, "bm_lon"]), float(t0.loc[c, "bm_lat"])
            r0 = s0.reindex(pd.MultiIndex.from_tuples([(c, y)]))
            out.append({"dataset": ds, "year": y, "country": c, "regions": int(ok.sum()),
                        "pop_share_regions_with_light": round(share, 4),
                        "err_light_weighted_km": float(gc(r0.bm_lon.iloc[0], r0.bm_lat.iloc[0], tl, tla)),
                        "err_pop_weighted_regions_km": float(gc(plon, plat, tl, tla))})
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Fyra kontroller efter granskningen av manuset.")
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--out", type=Path, default=Path("docs/robustness"))
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    rel = {lv: load_release(a.release, lv) for lv in (1, 2)}

    # 1. USA:s delstater
    years = {2000, 2010, 2020}
    tru = {y: truth(a.truth, "census", y, 1, 10) for y in years}
    states = []
    for y1, y2 in PAIRS["census"]:
        ids, t, e, b = change(tru, rel[1], "census", 1, y1, y2, 1, keep=lambda i: i.str.startswith("USA."))
        states.append(row("census", 1, y1, y2, ids, t, e, b, units="US states"))
    pd.DataFrame(states).to_csv(a.out / "robust_states.csv", index=False, float_format="%.4f")
    print(pd.DataFrame(states)[["from", "to", "n", "true_med", "dir_right", "dir_right_lo", "dir_right_hi",
                                "rho", "rho_lo", "rho_hi", "snr_gt1"]].round(3).to_string(index=False))

    # 2. Treårsmedel, 3. permutation (nivå 2; nivå 1 för Mexiko, Europa och USA:s delstater)
    r3, rp = [], []
    for ds, pairs in PAIRS.items():
        ys = {y for pr in pairs for y in pr}
        for lv in (1, 2):
            if ds == "census" and lv == 1:
                keep = lambda i: i.str.startswith("USA.")
            else:
                keep = None
            tru = {y: truth(a.truth, ds, y, lv, 10) for y in ys}
            for y1, y2 in pairs:
                for win in (1, 3):
                    ids, t, e, b = change(tru, rel[lv], ds, lv, y1, y2, win, keep)
                    r3.append(row(ds, lv, y1, y2, ids, t, e, b, window_years=win))
                    if win == 1:
                        rp.append({"dataset": ds, "level": lv, "from": y1, "to": y2, "n": len(ids),
                                   **permutation(t, e, b)})
    pd.DataFrame(r3).to_csv(a.out / "robust_3yr.csv", index=False, float_format="%.4f")
    pd.DataFrame(rp).to_csv(a.out / "robust_permutation.csv", index=False, float_format="%.4f")
    print(pd.DataFrame(r3)[["dataset", "level", "from", "to", "window_years", "n", "err_med", "dir_right",
                            "rho", "snr_gt1"]].round(3).to_string(index=False))
    print(pd.DataFrame(rp).round(3).to_string(index=False))

    # 4. Korrektion av landsfelet
    nc = []
    for ds, ys in {"census": [2000, 2010, 2020], "inegi": [2010, 2020], "ibge": [2022],
                   "geostat": [2011, 2021]}.items():
        nc += national_correction(a.truth, a.release, ds, ys)
    nc = pd.DataFrame(nc)
    nc.to_csv(a.out / "robust_national_correction.csv", index=False, float_format="%.3f")
    print(nc.round(2).to_string(index=False))
    print(nc.groupby("dataset")[["err_light_weighted_km", "err_pop_weighted_regions_km"]].median().round(2))


if __name__ == "__main__":
    main()
