#!/usr/bin/env python3
"""B – Jämförelsemått (docs/ANALYSIS_PLAN.md §6).

Hur står sig utgåvans ljuscentroid mot det en användare kan få i dag? Varje
jämförelsekälla körs genom samma centroidkod (tools/compare_sources.py) och
jämförs mot facit med samma läges- och förändringsmått som för ljuset, på nivå
1 och 2, bara för facitländerna. Geometrisk mittpunkt är med som jämförelse.

Källor i denna körning: ljus (utgåvan), GHS-POP (JRC GHSL, R2023A, 30"),
geometrisk mittpunkt. WorldPop och GHS-BUILT-S läggs till senare.

Förbehåll (planen §6, utökat under Avvikelser): GHS-POP bygger på folkräkning
och bebyggelse, inte på nattljus; den är inte oberoende av facit och förändring
mellan epoker är delvis interpolerad. WorldPop (senare) använder nattljus som
förklaringsvariabel och redovisas skilt. Jämförelsen visar vad en användare kan
få, inte ett oberoende test.

Mått: läge = median great circle-fel mot facit (km) med blockbootstrap;
förändring = de frysta måtten i tools/change_detection.py. Facitår knyts till
GHS-POP:s närmaste epok (2011→2010, 2021→2020, 2022→2020).

Utdata i --out: <ds>_benchmark_position.csv och <ds>_benchmark_change.csv.

Exempel
  python tools/compare_benchmarks.py --sources ghspop
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

import tools.change_detection as cd
import tools.hybrid_estimate as he
from tools.change_detection import en, metrics, block_ids, bootstrap, truth, release
from tools.hybrid_estimate import geometric_centroids, great_circle_km, boot_ci
from tools.compare_sources import nearest

# Blockbootstrapen återställs till samma frö före varje mått, så intervallen blir
# reproducerbara och oberoende av vilka och hur många källor som ingår (annars
# skulle en tillagd källa flytta de andras intervall via RNG-sekvensen).
BOOT_SEED = 20261121


def _reseed():
    import numpy as _np
    cd.RNG = _np.random.default_rng(BOOT_SEED)
    he.RNG = _np.random.default_rng(BOOT_SEED)

DS = {"census":  {"iso": ["USA"], "years": [2000, 2010, 2020],
                  "pairs": [(2000, 2010), (2010, 2020), (2000, 2020)]},
      "inegi":   {"iso": ["MEX"], "years": [2010, 2020], "pairs": [(2010, 2020)]},
      "ibge":    {"iso": ["BRA"], "years": [2022], "pairs": []},
      "geostat": {"iso": None,    "years": [2011, 2021], "pairs": [(2011, 2021)]}}


def facit_iso(cfg, truth_root, ds, lv):
    if cfg["iso"] is not None:
        return cfg["iso"]
    return sorted({u[:3] for u in truth(truth_root, ds, cfg["years"][0], lv, 10).index})


def light_series(release_path, lv, years):
    d = release(release_path, lv, set(years)).set_index(["unit_id", "year"])
    return {y: d.xs(y, level="year")[["bm_lon", "bm_lat"]] for y in years}


def source_series(results_dir, source, lv, epochs, years):
    """Facitår -> (unit_id -> bm_lon, bm_lat) ur jämförelsekällans epokfiler."""
    cache = {}
    out = {}
    for y in years:
        e = nearest(epochs, y)
        if e not in cache:
            f = results_dir / source / f"adm{lv}" / f"{e}.csv"
            cache[e] = pd.read_csv(f).set_index("unit_id")[["bm_lon", "bm_lat"]]
        out[y] = cache[e]
    return out


def geometric_series(cache_dir, bset, lv, isos):
    g = gpd.read_file(cache_dir / f"{bset}_adm{lv}.gpkg", engine="pyogrio")
    g = g.rename(columns={f"GID_{lv}": "unit_id"})
    g = g[g.unit_id.str[:3].isin(isos)]
    parts = []
    for iso, gi in g.groupby(g.unit_id.str[:3]):
        glon, glat = geometric_centroids(gi, iso)
        parts.append(pd.DataFrame({"unit_id": gi.unit_id.to_numpy(), "bm_lon": glon, "bm_lat": glat}))
    return pd.concat(parts, ignore_index=True).set_index("unit_id")


def coords(series, unit_ids):
    r = series.reindex(unit_ids)
    return r["bm_lon"].to_numpy(), r["bm_lat"].to_numpy()


def run_dataset(ds, cfg, truth_root, release_path, results_dir, cache_dir, bset,
                source_names, epochs, n_boot):
    pos_rows, chg_rows = [], []
    for lv in (1, 2):
        isos = facit_iso(cfg, truth_root, ds, lv)
        tru = {y: truth(truth_root, ds, y, lv, 10) for y in cfg["years"]}
        tru = {y: t[t.index.str[:3].isin(isos)] for y, t in tru.items()}
        series = {"light": light_series(release_path, lv, cfg["years"])}
        for s in source_names:
            series[s] = source_series(results_dir, s, lv, epochs, cfg["years"])
        geo = geometric_series(cache_dir, bset, lv, isos)
        names = ["light", *source_names, "geometric"]

        def est(name, y, ids):   # koordinater för en källa ett år
            s = geo if name == "geometric" else series[name][y]
            return coords(s, ids)

        # Läge per år
        for y in cfg["years"]:
            ids = tru[y].index
            tlon, tlat = tru[y].bm_lon.to_numpy(), tru[y].bm_lat.to_numpy()
            errs = {}
            for name in names:
                elon, elat = est(name, y, ids)
                errs[name] = great_circle_km(elon, elat, tlon, tlat)
            ok = np.isfinite(tlon)
            for v in errs.values():
                ok &= np.isfinite(v)
            kids = ids[ok]
            if len(kids) == 0:
                continue
            multi = kids.str[:3].nunique() >= 5
            blk = block_ids(kids, lv, multi)
            for name in names:
                e = errs[name][ok]
                _reseed()
                lo, hi = boot_ci(lambda idx, e=e: np.median(e[idx]), blk, n_boot)
                pos_rows.append({"dataset": ds, "level": lv, "year": y, "source": name,
                                 "n": len(kids), "err_med": float(np.median(e)),
                                 "err_med_lo": lo, "err_med_hi": hi})

        # Förändring per period
        for y1, y2 in cfg["pairs"]:
            ids = tru[y1].index.intersection(tru[y2].index)
            if len(ids) == 0:
                continue
            est_xy = {}
            okmask = np.ones(len(ids), bool)
            for name in names:
                l1 = np.column_stack(est(name, y1, ids))
                l2 = np.column_stack(est(name, y2, ids))
                est_xy[name] = (l1, l2)
                okmask &= np.isfinite(l1).all(1) & np.isfinite(l2).all(1)
            T1, T2 = tru[y1].loc[ids], tru[y2].loc[ids]
            okmask &= T1.bm_lon.notna().to_numpy() & T2.bm_lon.notna().to_numpy()
            kids = ids[okmask]
            if len(kids) == 0:
                continue
            multi = kids.str[:3].nunique() >= 5
            blk = block_ids(kids, lv, multi)
            t = np.column_stack(en(T1[okmask].bm_lon, T1[okmask].bm_lat,
                                   T2[okmask].bm_lon, T2[okmask].bm_lat))
            for name in names:
                l1, l2 = est_xy[name]
                e = np.column_stack(en(l1[okmask, 0], l1[okmask, 1], l2[okmask, 0], l2[okmask, 1]))
                m = metrics(t, e)
                _reseed()
                ci = bootstrap(t, e, blk, n_boot)
                chg_rows.append({"dataset": ds, "level": lv, "from": y1, "to": y2,
                                 "source": name, "n": len(kids), **m, **ci})
    return pd.DataFrame(pos_rows), pd.DataFrame(chg_rows)


def main(argv=None):
    p = argparse.ArgumentParser(description="B – jämförelsemått mot facit.")
    p.add_argument("--sources", nargs="+", default=["ghspop"])
    p.add_argument("--datasets", nargs="+", default=list(DS), choices=list(DS))
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--results", type=Path, default=Path("data/results/gadm41"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    from wpc.lights import available_years
    epochs = {s: available_years(s, Path("data/raw") / s, None) for s in a.sources}
    # alla källor förutsätts ha samma epoker här (GHS-POP). Ta första källans.
    ep = epochs[a.sources[0]]

    for ds in a.datasets:
        print(f"== {ds} ==")
        pos, chg = run_dataset(ds, DS[ds], a.truth, a.release, a.results, a.cache_dir,
                               a.boundary_set, a.sources, ep, a.n_boot)
        pos.to_csv(a.out / f"{ds}_benchmark_position.csv", index=False, float_format="%.4f")
        chg.to_csv(a.out / f"{ds}_benchmark_change.csv", index=False, float_format="%.4f")
        if len(pos):
            print(pos[["level", "year", "source", "n", "err_med", "err_med_lo", "err_med_hi"]].round(3).to_string(index=False))
        if len(chg):
            print(chg[["level", "from", "to", "source", "n", "err_med", "dir_right", "dir_right_lo", "dir_right_hi"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
