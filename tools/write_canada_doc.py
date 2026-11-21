#!/usr/bin/env python3
"""Skriver docs/extended/C_canada.md ur C:s CSV-filer (planen §7).

Alla tal hämtas ur docs/extended/C_canada_*.csv. Ändra aldrig talen för hand;
kör om tools/canada_confirm.py och sedan detta skript. Tolkningstabellen i §7
återges ordagrant; vilken rad utfallet hamnar på räknas ur huvudmåttets CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

SRC = {"light": "ljus (utgåvan)", "ghspop": "GHS-POP", "ghsbuilts": "GHS-BUILT-S",
       "worldpop": "WorldPop", "geometric": "geometrisk"}
ORDER = ["light", "ghspop", "ghsbuilts", "worldpop", "geometric"]

# §7-tabellen, ordagrant ur docs/ANALYSIS_PLAN.md.
TABLE = [
    ("Punktskattning ≥ 0,65 och intervallets nedre gräns > 0,50",
     "VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade perioderna (0,51–0,66)."),
    ("Punktskattning 0,50–0,65 och intervallbredd ≤ 0,15",
     "Resultatet från DMSP-tiden gäller också VIIRS."),
    ("Intervallbredd > 0,15 eller färre än 150 enheter",
     "Inkonklusivt. Ett andra facit inom VIIRS-tiden läggs till enligt samma plan, valt innan dess data öppnas."),
]


def f(x, nd=2):
    if pd.isna(x):
        return "–"
    return f"{x:.{nd}f}".replace(".", ",").replace("-", "−")


def read(d, name):
    p = d / f"C_canada_{name}.csv"
    return pd.read_csv(p) if p.exists() and p.stat().st_size > 1 else pd.DataFrame()


def outcome_row(m) -> int:
    """Vilken rad i §7-tabellen (1–3) utfallet hamnar på."""
    dr, lo, hi, n = m.dir_right, m.dir_right_lo, m.dir_right_hi, m.n
    width = m.ci_width
    if dr >= 0.65 and lo > 0.50:
        return 1
    if 0.50 <= dr <= 0.65 and width <= 0.15:
        return 2
    return 3   # intervallbredd > 0,15 eller färre än 150 enheter


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv C-resultatdokumentet ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    d = a.dir
    m = read(d, "main").iloc[0]
    change = read(d, "change")
    pos = read(d, "position")
    hyb = read(d, "hybrid")
    pair = read(d, "benchmark_pairs")
    row = outcome_row(m)

    L = ["# C – Kanada, bekräftande (planen §7)", "",
         "Frågan (docs/ANALYSIS_PLAN.md §7): håller resultatet om nattljusets "
         "centroider mäter förflyttning också för en förändring helt inom "
         "VIIRS-tiden (2016–2021)? Facit: Statistics Canada Geographic Attribute "
         "File 2016 och 2021, dissemination area-punkter viktade med DA-befolkning "
         "(se docs/DATA.md). Serie: utgåvan v2.0.0 (VIIRS 2016 och 2021, GADM 4.1 "
         "CAN).", "",
         "Måtten körs med de frysta funktionerna via drivrutinen "
         "`tools/canada_confirm.py` (de frysta skripten ändras inte). Att "
         "drivrutinen räknar som de frysta skripten är prövat mot census nivå 2 "
         "2010–2020 (`tests/test_canada_confirm.py`): n=2383, rätt riktning "
         "0,5820, median |e−t| 0,7947 återges exakt. Urvalet är som för övriga "
         "facit: complete och minst tio DA-punkter båda åren. Alla tal är skrivna "
         "ur `docs/extended/C_canada_*.csv` med `tools/write_canada_doc.py`.", "",
         "## Huvudmått: andel rätt riktning, nivå 2, 2016–2021", "",
         f"Blockbootstrap med block = provins ({int(m.n_blocks)} provinser), "
         f"{int(m.n)} census divisions, och permutationstest.", "",
         "| Mått | Värde |", "| --- | --- |",
         f"| Andel rätt riktning | **{f(m.dir_right)}** |",
         f"| 95 %-intervall (blockbootstrap) | {f(m.dir_right_lo)} – {f(m.dir_right_hi)} |",
         f"| Intervallbredd | {f(m.ci_width)} |",
         f"| Antal enheter (nivå 2) | {int(m.n)} |",
         f"| Permutationstest, nollfördelningens medel | {f(m.perm_null_mean)} ({f(m.perm_null_lo)}–{f(m.perm_null_hi)}) |",
         f"| Permutationstest, p (andel ≥ observerad) | {f(m.perm_p_upper, 4)} |", "",
         "### Tolkning, bestämd i förväg (§7, ordagrant)", "",
         "| Utfall | Slutsats |", "| --- | --- |"]
    for i, (utfall, slutsats) in enumerate(TABLE, 1):
        mark = " **←**" if i == row else ""
        L.append(f"| {utfall}{mark} | {slutsats} |")
    concl = TABLE[row - 1][1]
    L += ["",
          f"**Utfallet hamnar på rad {row}.** Punktskattningen är {f(m.dir_right)} "
          f"och intervallet {f(m.dir_right_lo)}–{f(m.dir_right_hi)} (bredd "
          f"{f(m.ci_width)}); antalet enheter är {int(m.n)}. "
          + ("Intervallbredden överstiger 0,15, så utfallet är: "
             if row == 3 else "Slutsatsen enligt tabellen är: ")
          + f"*{concl}*", "",
          "Riktningen ligger ändå över slumpen (permutationstestets nollmedel "
          f"{f(m.perm_null_mean)}, p = {f(m.perm_p_upper, 4)}), men intervallet är "
          "för brett för en bestämd slutsats enligt tabellen.", ""]

    # Bimått
    L += ["## Bimått", "", "**Förändring 2016–2021 (median |e−t| mot median |t|).**", "",
          "| Nivå | n | median \\|t\\| | median \\|e−t\\| | kvot | rätt riktning |",
          "| --- | --- | --- | --- | --- | --- |"]
    for _, r in change.iterrows():
        L.append(f"| {int(r.level)} | {int(r.n)} | {f(r.true_med)} | {f(r.err_med)} | {f(r.ratio)} | "
                 f"{f(r.dir_right)} ({f(r.dir_right_lo)}–{f(r.dir_right_hi)}) |")
    L += ["", "**Lägesfel 2021 (median great circle-fel mot facit, km).**", "",
          "| Nivå | källa | n | median fel | 95 %-intervall |", "| --- | --- | --- | --- | --- |"]
    for lv in (1, 2):
        for s in ORDER:
            r = pos[(pos.level == lv) & (pos.source == s)]
            if len(r):
                r = r.iloc[0]
                L.append(f"| {lv} | {SRC[s]} | {int(r.n)} | {f(r.err_med)} | {f(r.err_med_lo)}–{f(r.err_med_hi)} |")
    L += ["",
          "Ljuset och utgåvan använder exakt 2016 och 2021. Jämförelsekällornas "
          "epoker är de närmaste: GHS-POP och GHS-BUILT-S 2015 och 2020, WorldPop "
          "2016 och 2020 – alltså perioderna 2015–2020 respektive 2016–2020 mot "
          "facits 2016–2021 (Avvikelser §6/§7).", ""]
    if len(pair):
        L += ["**Parvis mot ljuset, per enhet (läge 2021).**", "",
              "| Nivå | källa | n | andel närmare facit än ljuset | median skillnad ljus−källa (km) |",
              "| --- | --- | --- | --- | --- |"]
        for lv in (1, 2):
            for s in ORDER:
                r = pair[(pair.level == lv) & (pair.source == s)]
                if len(r):
                    r = r.iloc[0]
                    L.append(f"| {lv} | {SRC[s]} | {int(r.n)} | {f(r.closer_share)} "
                             f"({f(r.closer_lo)}–{f(r.closer_hi)}) | {f(r.diff_med)} "
                             f"({f(r.diff_lo)}–{f(r.diff_hi)}) |")
        L.append("")
    if len(hyb):
        L += ["**Hybridskattning på provinsnivå (A1, beskrivande, 13 provinser).**", "",
              "| År | skattning | median lägesfel (km) |", "| --- | --- | --- |"]
        name = {"light": "ljus", "hybrid_light": "hybrid-ljus", "hybrid_geom": "hybrid-geom"}
        for _, r in hyb.iterrows():
            L.append(f"| {int(r.year)} | {name.get(r.estimate, r.estimate)} | {f(r.err_med)} |")
        L.append("")
    L += ["**A3 (detektionsmodellen) – redovisas inte.** Den frysta A3-modellens "
          "brusmått σ är definierat för perioder inom DMSP och över skarven, men "
          "inte för en period helt inom VIIRS (2016–2021). A3-modellen klarade "
          "dessutom inte sin egen prövning (H3b, H3c; se `A3_detection.md`). A3-"
          "bimåttet redovisas därför inte utan en uttrycklig σ-definition för "
          "VIIRS-perioden.", ""]

    L += ["## Förbehåll", "",
          "- Omkring 0,9 % av Kanadas befolkning ligger i DA-representativa punkter "
          "som faller utanför GADM-polygonerna (kustnära punkter).",
          "- WorldPop använder nattljus som förklaringsvariabel och är inte "
          "oberoende av det utgåvan mäter; redovisas tillsammans med de andra här "
          "men ska läsas med det i åtanke. GHS-POP och GHS-BUILT-S bygger på "
          "folkräkning respektive bebyggelse.",
          "- Jämförelsekällornas epoker (2015/2016–2020) skiljer sig från facits "
          "2016–2021; ljuset använder exakt 2016 och 2021.", ""]

    out = d / "C_canada.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}  (utfall: rad {row})")


if __name__ == "__main__":
    main()
