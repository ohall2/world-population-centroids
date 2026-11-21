#!/usr/bin/env python3
"""Skriver docs/extended/C_japan.md ur C:s japanska CSV-filer (planen §7).

Japan är det andra facit inom VIIRS-tiden (rad 3-uppföljningen). Huvudmåttet
redovisas först, ensamt mot §7-tolkningstabellen (ordagrant); rad 4 är tillagd i
efterhand (Avvikelser §7) och märks så. Därefter samma bimått som för Kanada men
utan A3, och Kanadas resultat bredvid. Alla tal hämtas ur docs/extended/C_japan_*.csv
och C_canada_main.csv. Ändra aldrig talen för hand.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tools.write_canada_doc import TABLE as TABLE3, f, SRC, ORDER

# Rad 4, tillagd i efterhand efter Japans resultat (Avvikelser §7), ordagrant.
ROW4 = ("Punktskattning < 0,50 och intervallbredd ≤ 0,15 *(tillagd i efterhand)*",
        "VIIRS mäter inte riktningen i detta facit. Ligger intervallets övre gräns "
        "under 0,50 redovisas utfallet som signifikant under slumpen, och orsaken som "
        "outredd.")
TABLE = TABLE3 + [ROW4]


def outcome_row(m) -> int:
    dr, lo, width = m.dir_right, m.dir_right_lo, m.ci_width
    if dr >= 0.65 and lo > 0.50:
        return 1
    if 0.50 <= dr <= 0.65 and width <= 0.15:
        return 2
    if dr < 0.50 and width <= 0.15:
        return 4
    return 3     # intervallbredd > 0,15 eller färre än 150 enheter


def read(d: Path, tag: str, name: str) -> pd.DataFrame:
    p = d / f"C_{tag}_{name}.csv"
    return pd.read_csv(p) if p.exists() and p.stat().st_size > 1 else pd.DataFrame()


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv C_japan.md ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    d = a.dir
    m = read(d, "japan", "main").iloc[0]
    change = read(d, "japan", "change")
    pos = read(d, "japan", "position")
    hyb = read(d, "japan", "hybrid")
    pair = read(d, "japan", "benchmark_pairs")
    row = outcome_row(m)
    can = read(d, "canada", "main")
    can_m = can.iloc[0] if len(can) else None

    L = ["# C (andra facit) – Japan, bekräftande (planen §7)", "",
         "Japan är det andra facit inom VIIRS-tiden, valt efter att Kanada hamnat "
         "på rad 3 (inkonklusivt). Frågan (§7): mäter nattljusets centroider "
         "förflyttning för en förändring helt inom VIIRS-tiden (2015–2020)? Facit: "
         "e-Stats rutnätsstatistik ur folkräkningen (räknade tal), 500 m-mesh, "
         "rutans mittpunkt viktad med rutans befolkning (se docs/DATA.md; den "
         "nationella totalen är certifierad mot folkräkningen båda åren). Serie: "
         "utgåvan v2.0.0 (VIIRS 2015 och 2020, GADM 4.1 JPN).", "",
         "Måtten körs med de frysta funktionerna via samma drivrutin som Kanada "
         "(`tools/canada_confirm.py`, facit som parameter); de frysta skripten "
         "ändras inte. Att drivrutinen räknar som dem är prövat mot census nivå 2 "
         "2010–2020 (`tests/test_canada_confirm.py`). Urvalet är complete och minst "
         "tio rutpunkter båda åren. Block = nivå 1-enhet (prefektur). Alla tal är "
         "skrivna ur `docs/extended/C_japan_*.csv` med `tools/write_japan_doc.py`.", "",
         "## Huvudmått: andel rätt riktning, nivå 2, 2015–2020", "",
         f"Blockbootstrap med block = prefektur ({int(m.n_blocks)} block), "
         f"{int(m.n)} kommuner, och permutationstest.", "",
         "| Mått | Värde |", "| --- | --- |",
         f"| Andel rätt riktning | **{f(m.dir_right)}** |",
         f"| 95 %-intervall (blockbootstrap) | {f(m.dir_right_lo)} – {f(m.dir_right_hi)} |",
         f"| Intervallbredd | {f(m.ci_width)} |",
         f"| Antal enheter (nivå 2) | {int(m.n)} |",
         f"| Permutationstest, nollfördelningens medel | {f(m.perm_null_mean)} ({f(m.perm_null_lo)}–{f(m.perm_null_hi)}) |",
         f"| Permutationstest, nedre svans p (andel ≤ observerad) | {f(m.perm_p_lower, 4)} |", "",
         "### Tolkning, bestämd i förväg (§7, ordagrant; rad 4 tillagd i efterhand)", "",
         "| Utfall | Slutsats |", "| --- | --- |"]
    for i, (utfall, slutsats) in enumerate(TABLE, 1):
        mark = " **←**" if i == row else ""
        L.append(f"| {utfall}{mark} | {slutsats} |")
    concl = TABLE[row - 1][1]
    L += ["",
          f"**Utfallet hamnar på rad {row}.** Punktskattningen är {f(m.dir_right)}, "
          f"intervallet {f(m.dir_right_lo)}–{f(m.dir_right_hi)} (bredd {f(m.ci_width)}), "
          f"n = {int(m.n)}. Slutsatsen enligt tabellen: *{concl}*", ""]
    if row == 4:
        L += [f"Intervallets övre gräns ({f(m.dir_right_hi)}) ligger under 0,50, och "
              f"permutationstestets nedre svans ger p = {f(m.perm_p_lower, 4)} (andel av "
              f"nollfördelningen ≤ den observerade, nollmedel {f(m.perm_null_mean)}). "
              "Utfallet redovisas därför som **signifikant under slumpen**. Orsaken är "
              "**outredd** – inga analyser görs här för att förklara den.", ""]

    # Kanada bredvid + samlad slutsats för VIIRS-frågan
    if can_m is not None:
        L += ["### Samlad slutsats för VIIRS-frågan", "",
              f"Kanada (2016–2021) blev **inkonklusivt** (rätt riktning "
              f"{f(can_m.dir_right)}, intervall {f(can_m.dir_right_lo)}–{f(can_m.dir_right_hi)}, "
              "rad 3). Japan (2015–2020) ligger **under slumpen** "
              f"({f(m.dir_right)}, rad 4). **Inget av de två facit ger stöd för att "
              "VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade serierna "
              "(0,51–0,66).** Stoppregeln (§7) gäller: inget tredje facit läggs till.", ""]

    # Bimått (samma som Kanada, utan A3)
    L += ["## Bimått", "", "**Förändring 2015–2020 (median |e−t| mot median |t|).**", "",
          "| Nivå | n | median \\|t\\| | median \\|e−t\\| | kvot | rätt riktning |",
          "| --- | --- | --- | --- | --- | --- |"]
    for _, r in change.iterrows():
        L.append(f"| {int(r.level)} | {int(r.n)} | {f(r.true_med)} | {f(r.err_med)} | {f(r.ratio)} | "
                 f"{f(r.dir_right)} ({f(r.dir_right_lo)}–{f(r.dir_right_hi)}) |")
    L += ["", "**Lägesfel 2020 (median great circle-fel mot facit, km).**", "",
          "| Nivå | källa | n | median fel | 95 %-intervall |", "| --- | --- | --- | --- | --- |"]
    for lv in (1, 2):
        for s in ORDER:
            r = pos[(pos.level == lv) & (pos.source == s)]
            if len(r):
                r = r.iloc[0]
                L.append(f"| {lv} | {SRC[s]} | {int(r.n)} | {f(r.err_med)} | {f(r.err_med_lo)}–{f(r.err_med_hi)} |")
    L += ["",
          "Ljuset/utgåvan och alla jämförelsekällor använder exakt facitåren 2015 "
          "och 2020 (GHS-POP/GHS-BUILT-S epok 2015 och 2020, WorldPop 2015 och "
          "2020) – ingen epokförskjutning.", ""]
    if len(pair):
        L += ["**Parvis mot ljuset, per enhet (läge 2020).**", "",
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
        L += ["**Hybridskattning på prefekturnivå (A1, beskrivande).**", "",
              "| År | skattning | median lägesfel (km) |", "| --- | --- | --- |"]
        name = {"light": "ljus", "hybrid_light": "hybrid-ljus", "hybrid_geom": "hybrid-geom"}
        for _, r in hyb.iterrows():
            L.append(f"| {int(r.year)} | {name.get(r.estimate, r.estimate)} | {f(r.err_med)} |")
        L.append("")
    L += ["*A3 (detektionsmodellen) redovisas inte, som för Kanada: σ är inte "
          "definierat för en period helt inom VIIRS och modellen föll på sin egen "
          "prövning (se `A3_detection.md`).*", ""]

    L += ["## Förbehåll", "",
          "- **2020 var ett pandemiår**, vilket kan ha påverkat folkräkningens läge "
          "och jämförbarheten 2015→2020.",
          "- WorldPop använder nattljus som förklaringsvariabel och är inte oberoende "
          "av det utgåvan mäter; GHS-POP och GHS-BUILT-S bygger på folkräkning "
          "respektive bebyggelse.",
          "- Omkring 0,5 % av befolkningen ligger i rutor vars mittpunkt faller "
          "utanför GADM-polygonerna (kustnära rutor).", ""]

    out = d / "C_japan.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}  (Japan: rad {row})")


if __name__ == "__main__":
    main()
