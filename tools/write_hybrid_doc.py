#!/usr/bin/env python3
"""Skriver docs/extended/A1_hybrid.md ur A1:s CSV-filer (planen §3).

Alla tal hämtas ur docs/extended/<ds>_hybrid_*.csv. Ändra aldrig talen för hand;
kör om tools/hybrid_estimate.py och sedan detta skript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DS_NAMES = {"census": "USA (delstater)", "inegi": "Mexiko (delstater)",
            "ibge": "Brasilien (delstater)", "geostat": "Europa (regioner)"}
EST_NAMES = {"light": "ljus", "hybrid_light": "hybrid-ljus", "hybrid_geom": "hybrid-geom"}


def f(x, nd=2) -> str:
    """Svensk decimalkomma."""
    if pd.isna(x):
        return "–"
    return f"{x:.{nd}f}".replace(".", ",")


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size <= 1:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def position_table(pos: pd.DataFrame) -> list[str]:
    out = ["| År | n | ljus | hybrid-ljus | hybrid-geom |", "| --- | --- | --- | --- | --- |"]
    for y, g in pos.groupby("year"):
        g = g.set_index("estimate")
        n = int(g["n"].iloc[0])
        cells = [f"{f(g.loc[e,'err_med'])} ({f(g.loc[e,'err_med_lo'])}–{f(g.loc[e,'err_med_hi'])})"
                 if e in g.index else "–" for e in ("light", "hybrid_light", "hybrid_geom")]
        out.append(f"| {int(y)} | {n} | " + " | ".join(cells) + " |")
    return out


def pairs_table(pair: pd.DataFrame) -> list[str]:
    out = ["| År | ljus − hybrid-ljus (H1a) | hybrid-geom − hybrid-ljus (H1c) |",
           "| --- | --- | --- |"]
    for y, g in pair.groupby("year"):
        g = g.set_index("pair")
        cells = []
        for key in ("light_minus_hybrid_light", "hybrid_geom_minus_hybrid_light"):
            r = g.loc[key]
            cells.append(f"{f(r['diff_med'])} ({f(r['diff_lo'])}–{f(r['diff_hi'])})")
        out.append(f"| {int(y)} | " + " | ".join(cells) + " |")
    return out


def change_table(chg: pd.DataFrame) -> list[str]:
    out = ["| Period | skattning | n | median \\|e−t\\| | ratio | ρ | lutning | rätt riktning |",
           "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for (fr, to), g in chg.groupby(["from", "to"]):
        g = g.set_index("series")
        for e in ("light", "hybrid_light", "hybrid_geom"):
            if e not in g.index:
                continue
            r = g.loc[e]
            out.append(f"| {int(fr)}–{int(to)} | {EST_NAMES[e]} | {int(r['n'])} | "
                       f"{f(r['err_med'])} | {f(r['ratio'])} | {f(r['rho'])} | {f(r['slope'])} | "
                       f"{f(r['dir_right'])} ({f(r['dir_right_lo'])}–{f(r['dir_right_hi'])}) |")
    return out


def verdicts(pair: pd.DataFrame, chg: pd.DataFrame) -> tuple[str, str, str]:
    """H1a, H1b, H1c som 'stöds (alla år/perioder)', 'blandat' eller '–'."""
    def judge(mask_series):
        vals = list(mask_series)
        if not vals:
            return "–"
        return "stöds" if all(vals) else ("stöds ej" if not any(vals) else "blandat")

    h1a = judge(pair[pair.pair == "light_minus_hybrid_light"].diff_lo > 0)
    h1c = judge(pair[pair.pair == "hybrid_geom_minus_hybrid_light"].diff_lo > 0)
    if len(chg):
        rows = []
        for _, g in chg.groupby(["from", "to"]):
            g = g.set_index("series")
            if "light" in g.index and "hybrid_light" in g.index:
                rows.append(g.loc["hybrid_light", "dir_right_lo"] > g.loc["light", "dir_right_hi"])
        h1b = judge(pd.Series(rows))
    else:
        h1b = "–"
    return h1a, h1b, h1c


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv A1-resultatdokumentet ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    p.add_argument("--datasets", nargs="+", default=["census", "inegi", "ibge", "geostat"])
    a = p.parse_args(argv)

    L = []
    L += ["# A1 – Hybridskattning (utforskande)", "",
          "Frågan (docs/ANALYSIS_PLAN.md §3): tillför ljuset något när "
          "befolkningstal per nivå 2-enhet är kända?", "",
          "Tre skattningar av en nivå 1-enhets centroid jämförs mot facit:",
          "*ljus* (utgåvans nivå 1-ljuscentroid), *hybrid-ljus* (befolkningsviktat "
          "3D-medel av nivå 2-enheternas ljuscentroider) och *hybrid-geom* (samma "
          "vikter, men nivå 2-enheternas geometriska ytcentroider). Hybriden "
          "använder alltså läget inom varje nivå 2-enhet från ljuset respektive "
          "geometrin, och känd befolkning per nivå 2-enhet som vikt.", "",
          "**Förbehåll.** Hybriden använder facits egna befolkningstal per nivå "
          "2-enhet. Den prövar hur mycket fel som återstår från läget inom "
          "enheterna, inte en oberoende skattning.", "",
          "Alla tal är skrivna ur `docs/extended/<facit>_hybrid_*.csv` med "
          "`tools/write_hybrid_doc.py`. Lägesmått: median great circle-fel mot "
          "facit (km) med 95 %-blockbootstrap. Parvis skillnad: per enhet, median "
          "med intervall (positiv = hybrid-ljus bättre). Förändringsmåtten är de "
          "frysta i `tools/change_detection.py`. Nivå 1-urval: enheter där nivå "
          "2-enheter med ljuscentroid täcker ≥ 99 % av facitbefolkningen.", ""]

    verdict_rows = []
    for ds in a.datasets:
        pos = read_csv(a.dir / f"{ds}_hybrid_position.csv")
        pair = read_csv(a.dir / f"{ds}_hybrid_position_pairs.csv")
        chg = read_csv(a.dir / f"{ds}_hybrid_change.csv")
        cov = read_csv(a.dir / f"{ds}_hybrid_coverage.csv")
        if pos.empty:
            continue
        L += [f"## {DS_NAMES.get(ds, ds)}", ""]
        L += ["**Läge: medianfel mot facit (km), med intervall.**", ""]
        L += position_table(pos) + [""]
        if not pair.empty:
            L += ["**Parvis felskillnad per enhet (km), med intervall.**", ""]
            L += pairs_table(pair) + [""]
        if not chg.empty:
            L += ["**Förändring (nivå 1, frysta mått).**", ""]
            L += change_table(chg) + [""]
        if not cov.empty:
            n_units = cov.adm1.nunique()
            drop = cov[(~cov.used) & (cov.n_adm2 > 0)]        # äkta 99 %-gallring
            none2 = cov[cov.n_adm2 == 0]                       # saknar nivå 2-indelning
            parts = []
            if len(drop):
                items = ", ".join(f"{r.adm1} {int(r.year)} ({f(r.cover)})" for _, r in drop.iterrows())
                parts.append(f"{len(drop)} enheter föll på 99 %-regeln (täckning inom "
                             f"parentes): {items}")
            else:
                parts.append("alla enheter med nivå 2-indelning klarade 99 %-regeln")
            if len(none2):
                isos = ", ".join(sorted(set(none2.adm1.str[:3])))
                parts.append(f"{none2.adm1.nunique()} nivå 1-enheter saknar nivå "
                             f"2-indelning och kan inte få en hybrid ({isos})")
            L += [f"**Urval.** Av {n_units} nivå 1-enheter med facit: " + "; ".join(parts) + ".", ""]
        h1a, h1b, h1c = verdicts(pair, chg)
        verdict_rows.append((DS_NAMES.get(ds, ds), h1a, h1b, h1c))

    L += ["## Hypoteser", "",
          "H1a: hybrid-ljus har mindre lägesfel än ljus. "
          "H1b: hybrid-ljus ger rätt riktning oftare än ljus (förändring). "
          "H1c: hybrid-ljus har mindre lägesfel än hybrid-geom.", "",
          "Stöds = intervallet pekar åt hypotesens håll för alla år/perioder "
          "(H1a/H1c: parvis skillnad > 0; H1b: hybrid-ljusets riktningsintervall "
          "över ljusets). \"–\" = ingen förändringsperiod.", "",
          "| Facit | H1a | H1b | H1c |", "| --- | --- | --- | --- |"]
    for name, h1a, h1b, h1c in verdict_rows:
        L += [f"| {name} | {h1a} | {h1b} | {h1c} |"]
    L += ["",
          "**Not om förändringen.** För förändringsmåtten (tabellerna ovan) är "
          "hybrid-geom lika bra som eller bättre än hybrid-ljus i alla facit. "
          "Ljuset bidrar alltså inte till förändringen när befolkningstalen per "
          "nivå 2-enhet är kända; det som bär förändringen är omviktningen mellan "
          "enheterna, inte var ljuset ligger inom dem.", ""]

    out = a.dir / "A1_hybrid.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}")


if __name__ == "__main__":
    main()
