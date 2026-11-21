#!/usr/bin/env python3
"""Skriver docs/extended/B_benchmarks.md ur B:s CSV-filer (planen §6).

Alla tal hämtas ur docs/extended/<facit>_benchmark_*.csv. Ändra aldrig talen för
hand; kör om tools/compare_sources.py, tools/compare_benchmarks.py och sedan
detta skript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DS = {"census": "USA", "inegi": "Mexiko", "ibge": "Brasilien", "geostat": "Europa"}
SRC = {"light": "ljus (utgåvan)", "ghspop": "GHS-POP", "ghsbuilts": "GHS-BUILT-S",
       "geometric": "geometrisk"}
ORDER = ["light", "ghspop", "ghsbuilts", "geometric"]


def f(x, nd=2) -> str:
    if pd.isna(x):
        return "–"
    return f"{x:.{nd}f}".replace(".", ",").replace("-", "−")


def read(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size <= 1:
        return pd.DataFrame()
    return pd.read_csv(path)


def position_table(pos: pd.DataFrame) -> list[str]:
    head = " | ".join(SRC[s] for s in ORDER)
    out = [f"| Nivå | År | {head} |", "| " + " | ".join(["---"] * (2 + len(ORDER))) + " |"]
    for (lv, y), g in pos.groupby(["level", "year"]):
        g = g.set_index("source")
        n = int(g["n"].iloc[0])
        cells = [f"{f(g.loc[s,'err_med'])} ({f(g.loc[s,'err_med_lo'])}–{f(g.loc[s,'err_med_hi'])})"
                 if s in g.index else "–" for s in ORDER]
        out.append(f"| {lv} | {y} (n={n}) | " + " | ".join(cells) + " |")
    return out


def change_table(chg: pd.DataFrame) -> list[str]:
    out = ["| Nivå | Period | källa | median \\|e−t\\| | rätt riktning |",
           "| --- | --- | --- | --- | --- |"]
    for (lv, fr, to), g in chg.groupby(["level", "from", "to"]):
        g = g.set_index("source")
        for s in ORDER:
            if s not in g.index:
                continue
            r = g.loc[s]
            out.append(f"| {lv} | {int(fr)}–{int(to)} | {SRC[s]} | {f(r['err_med'])} | "
                       f"{f(r['dir_right'])} ({f(r['dir_right_lo'])}–{f(r['dir_right_hi'])}) |")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv B-resultatdokumentet ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    p.add_argument("--datasets", nargs="+", default=["census", "inegi", "ibge", "geostat"])
    a = p.parse_args(argv)

    L = ["# B – Jämförelsemått (utforskande)", "",
         "Frågan (docs/ANALYSIS_PLAN.md §6): hur står sig utgåvans ljuscentroid mot "
         "det en användare kan få i dag? Varje jämförelsekälla körs genom samma "
         "centroidkod som ljuset och jämförs mot facit med samma läges- och "
         "förändringsmått, på nivå 1 och 2, bara för facitländerna. Geometrisk "
         "mittpunkt är med som jämförelse.", "",
         "**Denna körning:** ljus (utgåvan), **GHS-POP** och **GHS-BUILT-S** (JRC "
         "GHSL, R2023A, 30\", epokerna 2000/2010/2020; facitår knutet till närmaste "
         "epok). WorldPop läggs till i en senare körning.", "",
         "**Förbehåll.** GHS-POP bygger på folkräkning och GHS-BUILT-S på bebyggd "
         "yta, inte på nattljus; de är ändå inte oberoende av facit (GHS-POP "
         "fördelar folkräkningens folk), och förändringen mellan epoker är delvis "
         "interpolerad. Jämförelsen visar vad en användare kan få, inte ett "
         "oberoende test av nattljuset. WorldPop (senare körning) använder nattljus "
         "som förklaringsvariabel och är inte oberoende av ljuset; den redovisas "
         "skilt från GHS-POP och GHS-BUILT-S.", "",
         "Alla tal är skrivna ur `docs/extended/<facit>_benchmark_*.csv` med "
         "`tools/write_comparison_doc.py`. Läge: median great circle-fel mot facit "
         "(km); förändring: de frysta måtten i `tools/change_detection.py`; båda med "
         "95 %-blockbootstrap. Geometrisk mittpunkt är statisk, så dess "
         "förändringsskattning är alltid noll (rätt riktning 0).", ""]

    for ds in a.datasets:
        pos = read(a.dir / f"{ds}_benchmark_position.csv")
        chg = read(a.dir / f"{ds}_benchmark_change.csv")
        if pos.empty:
            continue
        L += [f"## {DS.get(ds, ds)}", "", "**Läge: median great circle-fel mot facit (km).**", ""]
        L += position_table(pos) + [""]
        if not chg.empty:
            L += ["**Förändring.**", ""]
            L += change_table(chg) + [""]

    out = a.dir / "B_benchmarks.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}")


if __name__ == "__main__":
    main()
