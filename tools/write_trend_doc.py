#!/usr/bin/env python3
"""Skriver docs/extended/A2_trend.md ur A2:s CSV-filer (planen §4).

Alla tal hämtas ur docs/extended/<ds>_trend_change.csv. Ändra aldrig talen för
hand; kör om tools/trend_estimate.py och sedan detta skript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DS_NAMES = {"census": "USA", "inegi": "Mexiko", "geostat": "Europa"}
SERIES_NAMES = {"endpoint": "ändpunkter", "three_year": "treårsmedel", "trend": "trend"}
H2_THRESHOLD = 0.20


def f(x, nd=2) -> str:
    if pd.isna(x):
        return "–"
    return f"{x:.{nd}f}".replace(".", ",")


def pct(x) -> str:
    if pd.isna(x):
        return "–"
    return f"{100 * x:.0f} %".replace("-", "−")


def read(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size <= 1:
        return pd.DataFrame()
    return pd.read_csv(path)


def change_table(df: pd.DataFrame, lv: int) -> list[str]:
    d = df[df.level == lv]
    out = ["| Period | skattning | n | median \\|e−t\\| | rätt riktning |",
           "| --- | --- | --- | --- | --- |"]
    for (fr, to), g in d.groupby(["from", "to"]):
        g = g.set_index("series")
        for s in ("endpoint", "three_year", "trend"):
            if s not in g.index:
                continue
            r = g.loc[s]
            out.append(f"| {int(fr)}–{int(to)} | {SERIES_NAMES[s]} | {int(r['n'])} | "
                       f"{f(r['err_med'])} ({f(r['err_med_lo'])}–{f(r['err_med_hi'])}) | "
                       f"{f(r['dir_right'])} ({f(r['dir_right_lo'])}–{f(r['dir_right_hi'])}) |")
    return out


def reductions(all_df: pd.DataFrame) -> pd.DataFrame:
    """Nivå 2: minskning av median |e−t| mot ändpunkter för treårsmedel och trend."""
    rows = []
    d = all_df[all_df.level == 2]
    for (ds, fr, to), g in d.groupby(["dataset", "from", "to"]):
        g = g.set_index("series")
        ep = g.loc["endpoint", "err_med"]
        rows.append({"dataset": ds, "from": int(fr), "to": int(to), "endpoint": ep,
                     "three_year": g.loc["three_year", "err_med"] if "three_year" in g.index else float("nan"),
                     "trend": g.loc["trend", "err_med"] if "trend" in g.index else float("nan"),
                     "red_3yr": (ep - g.loc["three_year", "err_med"]) / ep if "three_year" in g.index else float("nan"),
                     "red_trend": (ep - g.loc["trend", "err_med"]) / ep if "trend" in g.index else float("nan")})
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv A2-resultatdokumentet ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    p.add_argument("--datasets", nargs="+", default=["census", "inegi", "geostat"])
    a = p.parse_args(argv)

    frames = {ds: read(a.dir / f"{ds}_trend_change.csv") for ds in a.datasets}
    frames = {ds: d for ds, d in frames.items() if not d.empty}
    all_df = pd.concat(frames.values(), ignore_index=True) if frames else pd.DataFrame()

    L = ["# A2 – Trendskattning (utforskande)", "",
         "Frågan (docs/ANALYSIS_PLAN.md §4): minskar felet om hela tidsserien "
         "används i stället för två ändpunkter?", "",
         "Förflyttningen skattas per enhet och koordinat ur en Theil–Sen-lutning "
         "inom varje sensorperiod (DMSP t.o.m. 2012, VIIRS fr.o.m. 2013), i ett "
         "lokalt plan kring enhetens ljuscentroid 2013. Minst fyra år krävs i "
         "varje del. Jämförs med ändpunkter (nuvarande, ettårsfönster) och "
         "treårsmedel, på samma enheter. Ändpunkts- och treårsskattningarna och "
         "måtten är de frysta i `tools/change_detection.py` och "
         "`tools/robustness_checks.py`.", "",
         "Alla tal är skrivna ur `docs/extended/<facit>_trend_change.csv` med "
         "`tools/write_trend_doc.py`. Lägesändringsmått: median great circle-fel "
         "mot facit |e−t| (km) och andel rätt riktning, med 95 %-blockbootstrap.", ""]

    for ds, df in frames.items():
        L += [f"## {DS_NAMES.get(ds, ds)}", "", "**Nivå 1.**", ""]
        L += change_table(df, 1) + ["", "**Nivå 2.**", ""]
        L += change_table(df, 2) + [""]

    red = reductions(all_df)
    L += ["## Hypotes H2", "",
          "H2: trendskattningen minskar median |e−t| på nivå 2 med minst 20 % mot "
          "ändpunkter i minst fyra av de fem perioderna (treårsmedel gav 6–15 %).", "",
          "| Period | ändpunkter | treårsmedel (minskning) | trend (minskning) |",
          "| --- | --- | --- | --- |"]
    n_ge = 0
    for _, r in red.iterrows():
        if r["red_trend"] >= H2_THRESHOLD:
            n_ge += 1
        L += [f"| {DS_NAMES.get(r['dataset'], r['dataset'])} {r['from']}–{r['to']} | "
              f"{f(r['endpoint'])} | {f(r['three_year'])} ({pct(r['red_3yr'])}) | "
              f"{f(r['trend'])} ({pct(r['red_trend'])}) |"]
    n_tot = len(red)
    verdict = "stöds" if n_ge >= 4 else "stöds inte"
    L += ["",
          f"Trenden minskar felet med minst {pct(H2_THRESHOLD)} i {n_ge} av {n_tot} "
          f"perioder. **H2 {verdict}.** Treårsmedlets minskningar "
          f"({pct(red['red_3yr'].min())}–{pct(red['red_3yr'].max())}) ligger i linje "
          f"med det tidigare resultatet.", ""]

    out = a.dir / "A2_trend.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}  (H2 {verdict}: {n_ge}/{n_tot} perioder ≥ {pct(H2_THRESHOLD)})")


if __name__ == "__main__":
    main()
