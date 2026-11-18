#!/usr/bin/env python3
"""
World Population Centroids v2 – jämför nattljusets centroider med
befolkningens (facit från tools/validation_centroids.py).

Två mått, per serie, nivå och facitår (endast enheter med complete = True
och centroid i båda):

1. Lägesfel: avståndet (storcirkel, km) mellan seriens och facits centroid
   samma år, för Barmore (bm) och Aboufadel–Austin (aa). Median, 90 %,
   medelvärde, samt medelfelet i km öster/norr (systematisk förskjutning).

2. Förändringsfel (--change år1 år2): facits förflyttning mellan två år
   jämförs med seriens. Felet är längden av skillnaden mellan
   förflyttningsvektorerna (km). Som jämförelse anges facits egen
   förflyttning (median) och felet om serien inte hade rört sig alls
   (= facits förflyttning). En serie som fångar verklig förändring har
   förändringsfel klart mindre än facits förflyttning.

En serie anges som namn=källa:år-år[+källa:år-år], som i compare_series.py.

Exempel
  python tools/validate_series.py --dataset census --years 2000 2010 2020 \\
      --change 2000 2010 --change 2010 2020 \\
      --series reg=dmsp-irqr-reg:1992-2013+viirs-dmsp-reg:2014-2025 \\
               spliced=spliced-reg:1992-2025 \\
      --results data/results_mask2km data/results --out docs/validation
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_series import find, parse_series, read_results  # noqa: E402

R_KM = 6371.0


def en_km(lon0, lat0, lon1, lat1):
    dlon = (np.asarray(lon1) - np.asarray(lon0) + 180.0) % 360.0 - 180.0
    return (R_KM * np.radians(dlon) * np.cos(np.radians(lat0)),
            R_KM * np.radians(np.asarray(lat1) - np.asarray(lat0)))


def gc_km(lon1, lat1, lon2, lat2):
    l1, p1, l2, p2 = map(np.radians, (lon1, lat1, lon2, lat2))
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin((l2 - l1) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


_RELEASE: dict = {}


def release_year(release: Path, level: int, year: int) -> pd.DataFrame | None:
    """Den publicerade datamängdens centroider (tools/build_release.py) för ett år."""
    if level not in _RELEASE:
        f = release / f"wpc_v2_gadm41_adm{level}.csv.gz"
        if not f.exists():
            f = f.with_suffix("")  # okomprimerad .csv
        _RELEASE[level] = read_results(f) if f.exists() else None
    d = _RELEASE[level]
    if d is None or year not in set(d.year):
        return None
    return d[d.year == year][["unit_id", "aa_lon", "aa_lat", "bm_lon", "bm_lat"]].set_index("unit_id")


def series_year(results, boundary, parts, level, year) -> pd.DataFrame | None:
    for src, years in parts:
        if year in years:
            f = find(results, boundary, src, level) / f"{year}.csv"
            if not f.exists():
                return None
            d = read_results(f)
            return d[["unit_id", "aa_lon", "aa_lat", "bm_lon", "bm_lat"]].set_index("unit_id")
    return None


def truth(root: Path, dataset: str, year: int, level: int) -> pd.DataFrame:
    d = read_results(root / dataset / str(year) / f"adm{level}.csv")
    d = d[d.complete.astype(str).str.lower() == "true"]
    return d.set_index("unit_id")


def main(argv=None):
    p = argparse.ArgumentParser(description="Validera nattljusets centroider mot befolkningens.")
    p.add_argument("--dataset", choices=["census", "geostat", "inegi", "ibge"], required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--change", type=int, nargs=2, action="append", default=[])
    p.add_argument("--series", nargs="+", default=[])
    p.add_argument("--release", type=Path, default=None,
                   help="Validera även den publicerade datamängden i denna mapp (serien 'release').")
    p.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--results", type=Path, nargs="+", default=[Path("data/results_mask2km"), Path("data/results")])
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--min-pop", type=float, default=0, help="Bara enheter med minst så många invånare.")
    p.add_argument("--min-points", type=int, default=10,
                   help="Bara enheter med minst så många facitpunkter (standard 10; utesluter "
                        "småenheter som Vatikanen och Monaco, där 1 km-rutor ger fel befolkning).")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args(argv)
    specs = [parse_series(s) for s in args.series]
    if args.release:
        specs.insert(0, ("release", None))
    if not specs:
        p.error("ange --series och/eller --release")

    pos, chg = [], []
    for level in args.levels:
        T = {y: truth(args.truth, args.dataset, y, level) for y in args.years}
        T = {y: t[(t["pop"] >= args.min_pop) & (t["n_points"] >= args.min_points)] for y, t in T.items()}
        for name, parts in specs:
            S = {y: (release_year(args.release, level, y) if parts is None
                     else series_year(args.results, args.boundary_set, parts, level, y)) for y in args.years}
            for y in args.years:
                if S[y] is None:
                    continue
                j = T[y].join(S[y], rsuffix="_s", how="inner").dropna(subset=["bm_lon_s", "aa_lon_s"])
                row = {"series": name, "level": level, "year": y, "n": len(j)}
                for m in ("bm", "aa"):
                    d = gc_km(j[f"{m}_lon"], j[f"{m}_lat"], j[f"{m}_lon_s"], j[f"{m}_lat_s"])
                    e, n = en_km(j[f"{m}_lon"], j[f"{m}_lat"], j[f"{m}_lon_s"], j[f"{m}_lat_s"])
                    row.update({f"{m}_median": np.median(d), f"{m}_p90": np.quantile(d, 0.9),
                                f"{m}_mean": np.mean(d), f"{m}_east_mean": np.mean(e),
                                f"{m}_north_mean": np.mean(n)})
                pos.append(row)
            for y1, y2 in args.change:
                if S.get(y1) is None or S.get(y2) is None:
                    continue
                j = (T[y1].join(T[y2], lsuffix="_t1", rsuffix="_t2", how="inner")
                     .join(S[y1].add_suffix("_s1"), how="inner").join(S[y2].add_suffix("_s2"), how="inner")
                     .dropna())
                row = {"series": name, "level": level, "from": y1, "to": y2, "n": len(j)}
                for m in ("bm", "aa"):
                    te, tn = en_km(j[f"{m}_lon_t1"], j[f"{m}_lat_t1"], j[f"{m}_lon_t2"], j[f"{m}_lat_t2"])
                    se, sn = en_km(j[f"{m}_lon_s1"], j[f"{m}_lat_s1"], j[f"{m}_lon_s2"], j[f"{m}_lat_s2"])
                    err = np.hypot(se - te, sn - tn)
                    tm = np.hypot(te, tn)
                    r = np.corrcoef(np.r_[te, tn], np.r_[se, sn])[0, 1] if len(j) > 2 else np.nan
                    row.update({f"{m}_true_move_median": np.median(tm), f"{m}_err_median": np.median(err),
                                f"{m}_err_p90": np.quantile(err, 0.9), f"{m}_corr": r})
                chg.append(row)

    pd.set_option("display.width", 220)
    if pos:
        P = pd.DataFrame(pos)
        print("\n== Lägesfel (km), Barmore")
        print(P[["series", "level", "year", "n", "bm_median", "bm_p90", "bm_east_mean", "bm_north_mean"]]
              .round(2).to_string(index=False))
    if chg:
        C = pd.DataFrame(chg)
        print("\n== Förändringsfel (km), Barmore")
        print(C[["series", "level", "from", "to", "n", "bm_true_move_median", "bm_err_median", "bm_err_p90",
                 "bm_corr"]].round(3).to_string(index=False))
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        if pos:
            pd.DataFrame(pos).to_csv(args.out / f"{args.dataset}_position.csv", index=False, float_format="%.4f")
        if chg:
            pd.DataFrame(chg).to_csv(args.out / f"{args.dataset}_change.csv", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
