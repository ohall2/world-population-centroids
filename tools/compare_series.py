#!/usr/bin/env python3
"""
World Population Centroids v2 – jämför centroidserier.

För varje serie och nivå beräknas avståndet (storcirkel, km, R = 6371 km)
mellan en enhets Barmore-centroid ett år och nästa år ("årshopp"), och för
varje årspar median, 90- och 99-percentil över enheterna. En skarv mellan
sensorer syns som ett årspar med ovanligt stora hopp.

En serie kan bestå av flera källor som fogas ihop, t.ex.
"dmsp-irqr:1992-2013+viirs-dmsp:2014-2025". Varje källa läses från
<results>/<gränser>/<källa>/adm<nivå>/<år>.csv; flera resultatmappar kan
anges med --results, och den första som har källan används.

Dessutom beräknas avståndet mellan två källor samma år (--same-year).

Kräver bara pandas och numpy.

Exempel
  python tools/compare_series.py --levels 0 1 \\
      --series ours=dmsp-irqr:1992-2013+viirs-dmsp:2014-2025 \\
               li=li2020-dmsp:1992-2013+li2020-viirs:2014-2021 \\
      --same-year dmsp-irqr,viirs-dmsp,2013 \\
      --results data/results_mask2km data/results --out docs/comparison
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

R_KM = 6371.0


def read_results(path) -> pd.DataFrame:
    """Läs en resultatfil. ID-kolumner läses som text utan NA-tolkning: GADM
    4.1 har en enhet med ID:t "NA", som pandas annars gör till ett saknat värde."""
    df = pd.read_csv(path, keep_default_na=False, na_values={})
    for c in df.columns:
        if c in ("unit_id", "source", "COUNTRY") or c.startswith(("GID_", "NAME_", "TYPE_", "ENGTYPE_")):
            continue
        if not pd.api.types.is_numeric_dtype(df[c]):
            try:
                df[c] = pd.to_numeric(df[c].replace("", np.nan))
            except (ValueError, TypeError):
                pass
    return df


def gc_km(lon1, lat1, lon2, lat2):
    l1, p1, l2, p2 = map(np.radians, (lon1, lat1, lon2, lat2))
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin((l2 - l1) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def find(results: list[Path], boundary: str, source: str, level: int) -> Path:
    for r in results:
        d = r / boundary / source / f"adm{level}"
        if d.exists():
            return d
    raise FileNotFoundError(f"{source} adm{level} finns inte i {', '.join(map(str, results))}")


def load(results, boundary, source, level, years) -> pd.DataFrame:
    d = find(results, boundary, source, level)
    parts = []
    for y in years:
        f = d / f"{y}.csv"
        if not f.exists():
            raise FileNotFoundError(f)
        parts.append(read_results(f)[["unit_id", "year", "bm_lon", "bm_lat"]].assign(src=source))
    return pd.concat(parts, ignore_index=True)


def parse_series(spec: str) -> tuple[str, list[tuple[str, list[int]]]]:
    name, _, rest = spec.partition("=")
    parts = []
    for p in rest.split("+"):
        src, _, yr = p.partition(":")
        a, _, b = yr.partition("-")
        parts.append((src, list(range(int(a), int(b or a) + 1))))
    return name, parts


def jumps(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["bm_lon"]).sort_values(["unit_id", "year"])
    g = df.groupby("unit_id")
    df["km"] = gc_km(g.bm_lon.shift(), g.bm_lat.shift(), df.bm_lon, df.bm_lat)
    df["from_src"] = g.src.shift()
    ok = g.year.diff() == 1
    return df[ok]


def summarize(x: pd.Series) -> dict:
    return {"n": int(x.size), "median": x.median(), "p90": x.quantile(0.9), "p99": x.quantile(0.99)}


def main(argv=None):
    p = argparse.ArgumentParser(description="Jämför centroidserier (årshopp och skarvar).")
    p.add_argument("--results", type=Path, nargs="+", default=[Path("data/results")])
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--series", nargs="+", default=[])
    p.add_argument("--same-year", nargs="*", default=[],
                   help="källa1,källa2,år – avstånd mellan källorna samma år")
    p.add_argument("--out", type=Path, default=None, help="Mapp för CSV-filer med resultaten.")
    args = p.parse_args(argv)

    rows, srows = [], []
    for level in args.levels:
        for spec in args.series:
            name, parts = parse_series(spec)
            df = pd.concat([load(args.results, args.boundary_set, s, level, ys) for s, ys in parts])
            j = jumps(df)
            for (y0, y1), x in j.groupby([j.year - 1, j.year]):
                rows.append({"series": name, "level": level, "from": int(y0), "to": int(y1),
                             "sources": f"{x.from_src.iloc[0]}->{x.src.iloc[0]}", **summarize(x.km)})
        for spec in args.same_year:
            a, b, y = spec.split(",")
            da = load(args.results, args.boundary_set, a, level, [int(y)]).set_index("unit_id")
            db = load(args.results, args.boundary_set, b, level, [int(y)]).set_index("unit_id")
            jn = da.join(db, lsuffix="_a", rsuffix="_b", how="inner").dropna(subset=["bm_lon_a", "bm_lon_b"])
            d = pd.Series(gc_km(jn.bm_lon_a, jn.bm_lat_a, jn.bm_lon_b, jn.bm_lat_b), index=jn.index)
            srows.append({"level": level, "a": a, "b": b, "year": int(y), **summarize(d)})

    pd.set_option("display.width", 200)
    if rows:
        t = pd.DataFrame(rows)
        for (name, level), x in t.groupby(["series", "level"]):
            print(f"\n== {name} adm{level}: årshopp (km)")
            print(x[["from", "to", "sources", "n", "median", "p90", "p99"]].round(2).to_string(index=False))
    if srows:
        s = pd.DataFrame(srows)
        print("\n== Samma år, avstånd mellan källor (km)")
        print(s.round(2).to_string(index=False))
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        if rows:
            pd.DataFrame(rows).to_csv(args.out / "year_jumps.csv", index=False, float_format="%.4f")
        if srows:
            pd.DataFrame(srows).to_csv(args.out / "same_year.csv", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
