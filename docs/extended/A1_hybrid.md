# A1 – Hybridskattning (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §3): tillför ljuset något när befolkningstal per nivå 2-enhet är kända?

Tre skattningar av en nivå 1-enhets centroid jämförs mot facit:
*ljus* (utgåvans nivå 1-ljuscentroid), *hybrid-ljus* (befolkningsviktat 3D-medel av nivå 2-enheternas ljuscentroider) och *hybrid-geom* (samma vikter, men nivå 2-enheternas geometriska ytcentroider). Hybriden använder alltså läget inom varje nivå 2-enhet från ljuset respektive geometrin, och känd befolkning per nivå 2-enhet som vikt.

**Förbehåll.** Hybriden använder facits egna befolkningstal per nivå 2-enhet. Den prövar hur mycket fel som återstår från läget inom enheterna, inte en oberoende skattning.

Alla tal är skrivna ur `docs/extended/<facit>_hybrid_*.csv` med `tools/write_hybrid_doc.py`. Lägesmått: median great circle-fel mot facit (km) med 95 %-blockbootstrap. Parvis skillnad: per enhet, median med intervall (positiv = hybrid-ljus bättre). Förändringsmåtten är de frysta i `tools/change_detection.py`. Nivå 1-urval: enheter där nivå 2-enheter med ljuscentroid täcker ≥ 99 % av facitbefolkningen.

## USA (delstater)

**Läge: medianfel mot facit (km), med intervall.**

| År | n | ljus | hybrid-ljus | hybrid-geom |
| --- | --- | --- | --- | --- |
| 2000 | 50 | 13,15 (10,23–17,62) | 0,52 (0,35–0,66) | 3,02 (2,05–4,97) |
| 2010 | 50 | 10,53 (7,13–17,57) | 0,35 (0,25–0,59) | 2,66 (1,89–4,84) |
| 2020 | 50 | 12,14 (8,21–16,58) | 0,38 (0,28–0,60) | 2,57 (1,81–4,67) |

**Parvis felskillnad per enhet (km), med intervall.**

| År | ljus − hybrid-ljus (H1a) | hybrid-geom − hybrid-ljus (H1c) |
| --- | --- | --- |
| 2000 | 12,24 (8,12–17,41) | 2,29 (1,63–4,12) |
| 2010 | 10,33 (6,36–16,82) | 2,36 (1,63–4,48) |
| 2020 | 11,09 (7,87–15,68) | 2,35 (1,55–4,15) |

**Förändring (nivå 1, frysta mått).**

| Period | skattning | n | median \|e−t\| | ratio | ρ | lutning | rätt riktning |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2000–2010 | ljus | 50 | 5,50 | 2,01 | 0,10 | 0,72 | 0,50 (0,36–0,64) |
| 2000–2010 | hybrid-ljus | 50 | 0,30 | 0,11 | 0,98 | 0,98 | 0,96 (0,90–1,00) |
| 2000–2010 | hybrid-geom | 50 | 0,31 | 0,11 | 0,99 | 1,04 | 0,98 (0,94–1,00) |
| 2000–2020 | ljus | 50 | 6,56 | 1,69 | 0,16 | 0,31 | 0,58 (0,44–0,72) |
| 2000–2020 | hybrid-ljus | 50 | 0,36 | 0,09 | 0,99 | 0,98 | 0,98 (0,94–1,00) |
| 2000–2020 | hybrid-geom | 50 | 0,44 | 0,11 | 0,99 | 1,02 | 1,00 (1,00–1,00) |
| 2010–2020 | ljus | 50 | 5,67 | 2,18 | 0,02 | 0,02 | 0,54 (0,40–0,68) |
| 2010–2020 | hybrid-ljus | 50 | 0,20 | 0,08 | 0,97 | 0,98 | 0,98 (0,94–1,00) |
| 2010–2020 | hybrid-geom | 50 | 0,18 | 0,07 | 0,99 | 0,99 | 0,98 (0,94–1,00) |

**Urval.** Av 50 nivå 1-enheter med facit: alla enheter med nivå 2-indelning klarade 99 %-regeln.

## Mexiko (delstater)

**Läge: medianfel mot facit (km), med intervall.**

| År | n | ljus | hybrid-ljus | hybrid-geom |
| --- | --- | --- | --- | --- |
| 2010 | 32 | 8,88 (3,81–11,54) | 0,87 (0,41–1,02) | 2,86 (1,65–6,04) |
| 2020 | 32 | 7,38 (3,95–10,57) | 0,87 (0,55–1,20) | 2,74 (1,44–6,25) |

**Parvis felskillnad per enhet (km), med intervall.**

| År | ljus − hybrid-ljus (H1a) | hybrid-geom − hybrid-ljus (H1c) |
| --- | --- | --- |
| 2010 | 7,86 (3,41–10,84) | 1,51 (0,76–5,51) |
| 2020 | 5,63 (2,88–9,27) | 1,43 (0,62–5,18) |

**Förändring (nivå 1, frysta mått).**

| Period | skattning | n | median \|e−t\| | ratio | ρ | lutning | rätt riktning |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2010–2020 | ljus | 32 | 4,74 | 2,28 | 0,42 | 0,38 | 0,72 (0,56–0,88) |
| 2010–2020 | hybrid-ljus | 32 | 0,42 | 0,20 | 0,96 | 0,98 | 0,97 (0,91–1,00) |
| 2010–2020 | hybrid-geom | 32 | 0,25 | 0,12 | 0,95 | 0,93 | 1,00 (1,00–1,00) |

**Urval.** Av 32 nivå 1-enheter med facit: alla enheter med nivå 2-indelning klarade 99 %-regeln.

## Brasilien (delstater)

**Läge: medianfel mot facit (km), med intervall.**

| År | n | ljus | hybrid-ljus | hybrid-geom |
| --- | --- | --- | --- | --- |
| 2022 | 27 | 16,52 (11,97–19,63) | 0,44 (0,23–0,62) | 3,12 (2,21–6,10) |

**Parvis felskillnad per enhet (km), med intervall.**

| År | ljus − hybrid-ljus (H1a) | hybrid-geom − hybrid-ljus (H1c) |
| --- | --- | --- |
| 2022 | 15,89 (11,90–19,52) | 2,81 (1,79–5,64) |

**Urval.** Av 27 nivå 1-enheter med facit: alla enheter med nivå 2-indelning klarade 99 %-regeln.

## Europa (regioner)

**Läge: medianfel mot facit (km), med intervall.**

| År | n | ljus | hybrid-ljus | hybrid-geom |
| --- | --- | --- | --- | --- |
| 2011 | 341 | 2,97 (2,55–3,56) | 0,47 (0,34–0,63) | 1,57 (1,21–2,45) |
| 2021 | 419 | 2,47 (2,06–3,18) | 0,48 (0,34–0,62) | 1,66 (1,29–2,39) |

**Parvis felskillnad per enhet (km), med intervall.**

| År | ljus − hybrid-ljus (H1a) | hybrid-geom − hybrid-ljus (H1c) |
| --- | --- | --- |
| 2011 | 2,37 (1,75–2,85) | 1,00 (0,69–1,74) |
| 2021 | 1,73 (1,39–2,31) | 1,09 (0,82–1,83) |

**Förändring (nivå 1, frysta mått).**

| Period | skattning | n | median \|e−t\| | ratio | ρ | lutning | rätt riktning |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2011–2021 | ljus | 322 | 1,96 | 4,60 | -0,07 | -0,13 | 0,49 (0,43–0,54) |
| 2011–2021 | hybrid-ljus | 322 | 0,28 | 0,65 | 0,70 | 0,80 | 0,84 (0,79–0,88) |
| 2011–2021 | hybrid-geom | 322 | 0,09 | 0,22 | 0,92 | 0,85 | 0,92 (0,86–0,97) |

**Urval.** Av 451 nivå 1-enheter med facit: 3 enheter föll på 99 %-regeln (täckning inom parentes): EST.12_1 2011 (0,93), ROU.13_1 2011 (0,92), ROU.28_1 2011 (0,97); 13 nivå 1-enheter saknar nivå 2-indelning och kan inte få en hybrid (ALA, CYP, LIE, XAD).

## Hypoteser

H1a: hybrid-ljus har mindre lägesfel än ljus. H1b: hybrid-ljus ger rätt riktning oftare än ljus (förändring). H1c: hybrid-ljus har mindre lägesfel än hybrid-geom.

Stöds = intervallet pekar åt hypotesens håll för alla år/perioder (H1a/H1c: parvis skillnad > 0; H1b: hybrid-ljusets riktningsintervall över ljusets). "–" = ingen förändringsperiod.

| Facit | H1a | H1b | H1c |
| --- | --- | --- | --- |
| USA (delstater) | stöds | stöds | stöds |
| Mexiko (delstater) | stöds | stöds | stöds |
| Brasilien (delstater) | stöds | – | stöds |
| Europa (regioner) | stöds | stöds | stöds |

