#!/usr/bin/env python3
"""C – Kanada, bekräftande (docs/ANALYSIS_PLAN.md §7).

Separat drivrutin som IMPORTERAR de frysta funktionerna och tillämpar dem
oförändrade på Kanada. Inga frysta skript ändras. Att drivrutinen räknar som de
frysta skripten prövas i tests/test_canada_confirm.py (census nivå 2 2010–2020
återges exakt: n=2383, rätt riktning 0,582, median |e−t| 0,7947).

Huvudmått (§7): andel rätt riktning på nivå 2, 2016–2021, med blockbootstrap
(block = provins, automatiskt ur block_ids på nivå 2) och permutationstest.

Bimått (§7): lägesfel 2021 (median, nivå 1–2); median |e−t| mot median |t|
(2016–2021); hybridskattningen på provinsnivå (A1, beskrivande); jämförelsemåtten
(B). A3 hanteras särskilt (se dokumentet och drivrutinens utskrift): den frysta
A3-modellens brusmått σ är inte definierat för en period helt inom VIIRS
(2016–2021), och A3-modellen klarade inte sin egen prövning (H3b, H3c), så A3
redovisas inte här utan en uttrycklig σ-definition.

Facit i data/validation_canada/ (inte den frysta data/validation/). Urvalet är
som för övriga facit: complete och minst tio punkter (DA) båda åren.

Utdata i docs/extended/: C_canada_{main,position,change,hybrid,benchmarks,
benchmark_pairs}.csv.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import tools.change_detection as cd
import tools.hybrid_estimate as he
from tools.change_detection import en, metrics, block_ids, bootstrap, truth, release
from tools.robustness_checks import permutation
from tools.hybrid_estimate import boot_ci, great_circle_km, geometry_table, per_year_table, geometric_centroids
from tools.compare_benchmarks import pairwise

BOOT_SEED = 20261121
ISO = "CAN"
Y1, Y2 = 2016, 2021
# Jämförelsekällor (bimått). Facitår knyts till närmaste epok; se Avvikelser §6/§7.
BENCH = ["ghspop", "ghsbuilts", "worldpop"]


def _reseed():
    cd.RNG = np.random.default_rng(BOOT_SEED)
    he.RNG = np.random.default_rng(BOOT_SEED)


# --- förändring (som change_detection) -------------------------------------

def change_inputs(truth_root, release_path, ds, y1, y2, level, min_points=10, iso=None):
    """ids, t (facit-flytt), e (ljusflytt ur utgåvan) och block, exakt som
    tools/change_detection.py: complete & minst min_points punkter båda åren,
    enheter som finns i utgåvan båda åren, t/e via en()."""
    tru = {y: truth(truth_root, ds, y, level, min_points) for y in (y1, y2)}
    if iso is not None:
        tru = {y: t[t.index.str[:3] == iso] for y, t in tru.items()}
    ids = tru[y1].index.intersection(tru[y2].index)
    s = release(release_path, level, {y1, y2}).set_index(["unit_id", "year"])
    k1 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y1)])
    k2 = pd.MultiIndex.from_arrays([ids, np.full(len(ids), y2)])
    s1, s2 = s.reindex(k1), s.reindex(k2)
    ok = s1.bm_lon.notna().values & s2.bm_lon.notna().values
    ids = ids[ok]
    T1, T2, S1, S2 = tru[y1].loc[ids], tru[y2].loc[ids], s1[ok], s2[ok]
    t = np.column_stack(en(T1.bm_lon, T1.bm_lat, T2.bm_lon, T2.bm_lat))
    e = np.column_stack(en(S1.bm_lon.values, S1.bm_lat.values, S2.bm_lon.values, S2.bm_lat.values))
    multi = ids.str[:3].nunique() >= 5
    return ids, t, e, block_ids(ids, level, multi)


def change_measure(truth_root, release_path, ds, y1, y2, level, min_points=10,
                   iso=None, n_boot=1000):
    ids, t, e, blocks = change_inputs(truth_root, release_path, ds, y1, y2, level, min_points, iso)
    m = metrics(t, e)
    _reseed()
    ci = bootstrap(t, e, blocks, n_boot)
    return {"ids": ids, "t": t, "e": e, "blocks": blocks, "n": len(ids),
            "n_blocks": int(len(np.unique(blocks))), "metrics": m, "ci": ci}


# --- läge (median great circle-fel mot facit) ------------------------------

def position_measure(truth_root, release_path, source_year_pts, year, level, n_boot, source="light"):
    """Median lägesfel mot facit för en källa. source_year_pts: None för ljuset
    (läses ur utgåvan) eller {unit_id:(lon,lat)} för en jämförelsekälla."""
    t = truth(truth_root, "statcan", year, level, 10)
    t = t[t.index.str[:3] == ISO]
    if source == "light":
        rel = release(release_path, level, {year})
        rel = rel[rel.year == year].set_index("unit_id")
        est = rel.reindex(t.index)[["bm_lon", "bm_lat"]]
    else:
        est = pd.DataFrame({"bm_lon": pd.Series({u: p[0] for u, p in source_year_pts.items()}),
                            "bm_lat": pd.Series({u: p[1] for u, p in source_year_pts.items()})}).reindex(t.index)
    err = great_circle_km(est.bm_lon.values, est.bm_lat.values, t.bm_lon.values, t.bm_lat.values)
    ok = np.isfinite(err)
    ids, err = t.index[ok], err[ok]
    multi = ids.str[:3].nunique() >= 5
    blk = block_ids(ids, level, multi)
    _reseed()
    lo, hi = boot_ci(lambda idx, e=err: np.median(e[idx]), blk, n_boot)
    return {"source": source, "level": level, "year": year, "n": int(len(ids)),
            "err_med": float(np.median(err)), "err_med_lo": lo, "err_med_hi": hi,
            "err": err, "ids": ids, "blk": blk}


# --- jämförelsekällornas centroider för Kanada (B) -------------------------

def canada_source_points(source, level, years, cache_dir, raw_dir):
    """{facitår: {unit_id: (lon, lat)}} ur en jämförelsekälla, för Kanadas
    GADM-enheter, via samma centroidkod (compute_unit). Facitår -> närmaste epok."""
    import geopandas as gpd
    from wpc.lights import make_source, available_years
    from wpc.centroids import compute_unit
    from tools.compare_sources import nearest
    g = gpd.read_file(cache_dir / f"gadm41_adm{level}.gpkg", engine="pyogrio")
    g = g.rename(columns={f"GID_{level}": "unit_id"})
    g = g[g.unit_id.str[:3] == ISO].set_index("unit_id")
    epochs = available_years(source, raw_dir, None)
    year_epoch = {y: nearest(epochs, y) for y in years}
    by_epoch = {}
    for e in sorted(set(year_epoch.values())):
        src = make_source(source, e, raw_dir, None, None)
        pts = {}
        for uid, geom in g.geometry.items():
            r = compute_unit(geom, src)
            pts[uid] = (r.bm_lon, r.bm_lat)
        src.close()
        by_epoch[e] = pts
    return {y: by_epoch[year_epoch[y]] for y in years}, year_epoch


def main(argv=None):
    p = argparse.ArgumentParser(description="C – Kanada, bekräftande (drivrutin).")
    p.add_argument("--truth", type=Path, default=Path("data/validation_canada"))
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--gadm2", type=Path, default=Path("data/processed/boundaries/gadm41_adm2.gpkg"))
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    # ---- Huvudmått: rätt riktning nivå 2, 2016–2021, block=provins, + permutation
    r2 = change_measure(a.truth, a.release, "statcan", Y1, Y2, 2, iso=ISO, n_boot=a.n_boot)
    _reseed()
    perm = permutation(r2["t"], r2["e"], r2["blocks"])
    ci = r2["ci"]
    width = ci["dir_right_hi"] - ci["dir_right_lo"]
    main_row = {"measure": "dir_right_adm2_2016_2021", "n": r2["n"], "n_blocks": r2["n_blocks"],
                "dir_right": r2["metrics"]["dir_right"], "dir_right_lo": ci["dir_right_lo"],
                "dir_right_hi": ci["dir_right_hi"], "ci_width": width,
                "perm_null_mean": perm["dir_right_null_mean"], "perm_null_lo": perm["dir_right_null_lo"],
                "perm_null_hi": perm["dir_right_null_hi"], "perm_p_upper": perm["p_dir_right_upper"]}
    pd.DataFrame([main_row]).to_csv(a.out / "C_canada_main.csv", index=False, float_format="%.4f")
    print(f"HUVUDMÅTT nivå 2 2016–2021: n={r2['n']} block(provins)={r2['n_blocks']}  "
          f"rätt riktning {r2['metrics']['dir_right']:.4f} "
          f"[{ci['dir_right_lo']:.4f}, {ci['dir_right_hi']:.4f}] bredd {width:.4f}  "
          f"perm null {perm['dir_right_null_mean']:.4f} p(övre)={perm['p_dir_right_upper']:.4f}")

    # ---- Förändring nivå 1–2: median |e−t| mot median |t|, rätt riktning
    chg_rows = []
    for lv in (1, 2):
        r = r2 if lv == 2 else change_measure(a.truth, a.release, "statcan", Y1, Y2, 1, iso=ISO, n_boot=a.n_boot)
        m, c = r["metrics"], r["ci"]
        chg_rows.append({"level": lv, "from": Y1, "to": Y2, "n": r["n"], "true_med": m["true_med"],
                         "err_med": m["err_med"], "ratio": m["ratio"], "dir_right": m["dir_right"],
                         "dir_right_lo": c["dir_right_lo"], "dir_right_hi": c["dir_right_hi"]})
    pd.DataFrame(chg_rows).to_csv(a.out / "C_canada_change.csv", index=False, float_format="%.4f")

    # ---- A1 hybrid på provinsnivå (beskrivande), 2016 och 2021
    geo = geometry_table(a.gadm2, [ISO])
    cfg = {"iso": [ISO], "years": [Y1, Y2], "pairs": [(Y1, Y2)]}
    hyb_rows = []
    for year in (Y1, Y2):
        E, _cov = per_year_table(cfg, "statcan", year, a.truth, a.release, geo, 10)
        for est, (lo, la) in {"light": ("L_lon", "L_lat"), "hybrid_light": ("hyl_lon", "hyl_lat"),
                              "hybrid_geom": ("hyg_lon", "hyg_lat")}.items():
            err = great_circle_km(E[lo].to_numpy(), E[la].to_numpy(), E.tru_lon.to_numpy(), E.tru_lat.to_numpy())
            err = err[np.isfinite(err)]
            hyb_rows.append({"year": year, "level": 1, "estimate": est, "n": int(len(err)),
                             "err_med": float(np.median(err))})
    pd.DataFrame(hyb_rows).to_csv(a.out / "C_canada_hybrid.csv", index=False, float_format="%.4f")

    # ---- Lägesfel 2021 (bimått): ljuset + jämförelsekällorna (B), nivå 1–2,
    #      och parvis jämförelse per enhet mot ljuset.
    import geopandas as gpd
    n_boot = a.n_boot
    pos_rows, pair_rows, epochs_note = [], [], {}
    for lv in (1, 2):
        light = position_measure(a.truth, a.release, None, Y2, lv, n_boot, "light")
        pos_rows.append({k: light[k] for k in ("source", "level", "year", "n", "err_med", "err_med_lo", "err_med_hi")})
        el, lids = light["err"], light["ids"]
        # källpunkter för år 2021: jämförelsekällorna + geometrisk mittpunkt (B:s fjärde)
        src_pts = {}
        for source in BENCH:
            raw_dir = Path("data/raw") / source
            if not raw_dir.exists():
                continue
            pts, year_epoch = canada_source_points(source, lv, [Y1, Y2], a.cache_dir, raw_dir)
            epochs_note[source] = year_epoch
            src_pts[source] = pts[Y2]
        g = gpd.read_file(a.cache_dir / f"gadm41_adm{lv}.gpkg", engine="pyogrio").rename(columns={f"GID_{lv}": "unit_id"})
        g = g[g.unit_id.str[:3] == ISO]
        glon, glat = geometric_centroids(g, ISO)
        src_pts["geometric"] = {u: (lo, la) for u, lo, la in zip(g.unit_id.to_numpy(), glon, glat)}

        for source, pts2 in src_pts.items():
            pm = position_measure(a.truth, a.release, pts2, Y2, lv, n_boot, source)
            pos_rows.append({k: pm[k] for k in ("source", "level", "year", "n", "err_med", "err_med_lo", "err_med_hi")})
            common = lids.intersection(pm["ids"])
            if len(common):
                ei = pd.Series(el, index=lids).loc[common].to_numpy()
                es = pd.Series(pm["err"], index=pm["ids"]).loc[common].to_numpy()
                closer, diff = pairwise(ei, es)
                multi = common.str[:3].nunique() >= 5
                blk = block_ids(common, lv, multi)
                _reseed()
                cs_lo, cs_hi = boot_ci(lambda idx, s=es, l=ei: float(np.mean(s[idx] < l[idx])), blk, n_boot)
                _reseed()
                dm_lo, dm_hi = boot_ci(lambda idx, s=es, l=ei: float(np.median(l[idx] - s[idx])), blk, n_boot)
                pair_rows.append({"source": source, "level": lv, "year": Y2, "n": int(len(common)),
                                  "closer_share": closer, "closer_lo": cs_lo, "closer_hi": cs_hi,
                                  "diff_med": diff, "diff_lo": dm_lo, "diff_hi": dm_hi})
    pd.DataFrame(pos_rows).to_csv(a.out / "C_canada_position.csv", index=False, float_format="%.4f")
    pd.DataFrame(pair_rows).to_csv(a.out / "C_canada_benchmark_pairs.csv", index=False, float_format="%.4f")
    print("epokval (facitår -> epok):", epochs_note)
    print("klart ->", a.out)


if __name__ == "__main__":
    main()
