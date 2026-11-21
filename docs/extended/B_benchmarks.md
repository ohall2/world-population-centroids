# B – Jämförelsemått (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §6): hur står sig utgåvans ljuscentroid mot det en användare kan få i dag? Varje jämförelsekälla körs genom samma centroidkod som ljuset och jämförs mot facit med samma läges- och förändringsmått, på nivå 1 och 2, bara för facitländerna. Geometrisk mittpunkt är med som jämförelse.

**Denna körning:** ljus (utgåvan), **GHS-POP** och **GHS-BUILT-S** (JRC GHSL, R2023A, 30", epokerna 2000/2010/2020; facitår knutet till närmaste epok). WorldPop läggs till i en senare körning.

**Förbehåll.** GHS-POP bygger på folkräkning och GHS-BUILT-S på bebyggd yta, inte på nattljus; de är ändå inte oberoende av facit (GHS-POP fördelar folkräkningens folk), och förändringen mellan epoker är delvis interpolerad. Jämförelsen visar vad en användare kan få, inte ett oberoende test av nattljuset. WorldPop (senare körning) använder nattljus som förklaringsvariabel och är inte oberoende av ljuset; den redovisas skilt från GHS-POP och GHS-BUILT-S.

Alla tal är skrivna ur `docs/extended/<facit>_benchmark_*.csv` med `tools/write_comparison_doc.py`. Läge: median great circle-fel mot facit (km); förändring: de frysta måtten i `tools/change_detection.py`; båda med 95 %-blockbootstrap. Geometrisk mittpunkt är statisk, så dess förändringsskattning är alltid noll (rätt riktning 0).

## USA

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2000 (n=50) | 13,15 (10,23–17,62) | 0,31 (0,23–0,46) | 16,03 (11,99–21,94) | 62,23 (35,53–96,08) |
| 1 | 2010 (n=50) | 10,53 (6,92–16,93) | 0,18 (0,11–0,29) | 16,35 (12,27–22,31) | 61,20 (35,60–97,40) |
| 1 | 2020 (n=50) | 12,14 (8,19–16,58) | 2,04 (1,48–2,46) | 17,18 (13,80–20,12) | 63,39 (37,55–101,90) |
| 2 | 2000 (n=2441) | 1,51 (1,41–1,63) | 0,40 (0,37–0,44) | 1,22 (1,11–1,35) | 5,40 (4,72–6,25) |
| 2 | 2010 (n=2400) | 1,52 (1,43–1,62) | 0,09 (0,08–0,10) | 1,19 (1,08–1,34) | 5,50 (4,83–6,42) |
| 2 | 2020 (n=2436) | 1,46 (1,37–1,60) | 0,29 (0,26–0,32) | 1,27 (1,14–1,41) | 5,58 (4,94–6,55) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2000–2010 | ljus (utgåvan) | 5,50 | 0,50 (0,36–0,64) |
| 1 | 2000–2010 | GHS-POP | 0,18 | 0,98 (0,94–1,00) |
| 1 | 2000–2010 | GHS-BUILT-S | 2,11 | 0,78 (0,66–0,88) |
| 1 | 2000–2010 | geometrisk | 2,74 | 0,00 (0,00–0,00) |
| 1 | 2000–2020 | ljus (utgåvan) | 6,56 | 0,58 (0,44–0,72) |
| 1 | 2000–2020 | GHS-POP | 1,77 | 0,88 (0,78–0,96) |
| 1 | 2000–2020 | GHS-BUILT-S | 3,01 | 0,76 (0,62–0,86) |
| 1 | 2000–2020 | geometrisk | 3,88 | 0,00 (0,00–0,00) |
| 1 | 2010–2020 | ljus (utgåvan) | 5,67 | 0,54 (0,40–0,68) |
| 1 | 2010–2020 | GHS-POP | 1,86 | 0,74 (0,62–0,86) |
| 1 | 2010–2020 | GHS-BUILT-S | 1,89 | 0,70 (0,58–0,82) |
| 1 | 2010–2020 | geometrisk | 2,60 | 0,00 (0,00–0,00) |
| 2 | 2000–2010 | ljus (utgåvan) | 0,74 | 0,61 (0,59–0,64) |
| 2 | 2000–2010 | GHS-POP | 0,37 | 0,57 (0,54–0,59) |
| 2 | 2000–2010 | GHS-BUILT-S | 0,33 | 0,69 (0,67–0,71) |
| 2 | 2000–2010 | geometrisk | 0,37 | 0,00 (0,00–0,00) |
| 2 | 2000–2020 | ljus (utgåvan) | 1,05 | 0,66 (0,64–0,69) |
| 2 | 2000–2020 | GHS-POP | 0,53 | 0,59 (0,57–0,60) |
| 2 | 2000–2020 | GHS-BUILT-S | 0,48 | 0,70 (0,68–0,72) |
| 2 | 2000–2020 | geometrisk | 0,54 | 0,00 (0,00–0,00) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,79 | 0,58 (0,56–0,60) |
| 2 | 2010–2020 | GHS-POP | 0,26 | 0,59 (0,57–0,62) |
| 2 | 2010–2020 | GHS-BUILT-S | 0,25 | 0,65 (0,63–0,68) |
| 2 | 2010–2020 | geometrisk | 0,26 | 0,00 (0,00–0,00) |

## Mexiko

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2010 (n=32) | 8,88 (4,25–11,54) | 0,78 (0,44–1,12) | 8,10 (4,85–10,41) | 40,77 (26,57–48,31) |
| 1 | 2020 (n=32) | 7,38 (3,95–9,84) | 3,02 (1,68–4,34) | 8,64 (5,82–12,35) | 39,83 (27,00–48,58) |
| 2 | 2010 (n=1950) | 1,12 (0,98–1,27) | 0,35 (0,33–0,38) | 0,85 (0,73–1,03) | 3,19 (2,70–4,31) |
| 2 | 2020 (n=1947) | 1,03 (0,91–1,23) | 0,44 (0,40–0,51) | 0,87 (0,73–1,15) | 3,28 (2,74–4,43) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2010–2020 | ljus (utgåvan) | 4,74 | 0,72 (0,56–0,88) |
| 1 | 2010–2020 | GHS-POP | 3,22 | 0,41 (0,22–0,59) |
| 1 | 2010–2020 | GHS-BUILT-S | 3,16 | 0,62 (0,47–0,78) |
| 1 | 2010–2020 | geometrisk | 2,07 | 0,00 (0,00–0,00) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,85 | 0,56 (0,53–0,58) |
| 2 | 2010–2020 | GHS-POP | 0,31 | 0,51 (0,47–0,54) |
| 2 | 2010–2020 | GHS-BUILT-S | 0,49 | 0,50 (0,46–0,53) |
| 2 | 2010–2020 | geometrisk | 0,25 | 0,00 (0,00–0,00) |

## Brasilien

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2022 (n=27) | 16,52 (11,97–19,63) | 2,51 (0,94–4,60) | 22,23 (15,53–26,97) | 98,29 (75,69–141,35) |
| 2 | 2022 (n=5567) | 1,11 (0,99–1,32) | 0,34 (0,30–0,43) | 0,86 (0,78–1,02) | 3,82 (3,31–4,80) |

## Europa

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2011 (n=348) | 2,97 (2,54–3,62) | 0,20 (0,15–0,25) | 2,37 (1,81–3,09) | 8,67 (6,77–12,61) |
| 1 | 2021 (n=331) | 2,42 (2,05–3,25) | 0,44 (0,36–0,54) | 2,61 (2,22–3,50) | 9,08 (7,04–13,30) |
| 2 | 2011 (n=6266) | 0,95 (0,89–1,01) | 0,22 (0,21–0,24) | 0,46 (0,42–0,52) | 1,79 (1,64–1,99) |
| 2 | 2021 (n=6247) | 0,71 (0,66–0,76) | 0,27 (0,25–0,28) | 0,47 (0,43–0,53) | 1,81 (1,65–2,00) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2011–2021 | ljus (utgåvan) | 1,97 | 0,49 (0,43–0,53) |
| 1 | 2011–2021 | GHS-POP | 0,37 | 0,71 (0,63–0,80) |
| 1 | 2011–2021 | GHS-BUILT-S | 0,41 | 0,56 (0,51–0,62) |
| 1 | 2011–2021 | geometrisk | 0,42 | 0,00 (0,00–0,00) |
| 2 | 2011–2021 | ljus (utgåvan) | 0,64 | 0,51 (0,49–0,52) |
| 2 | 2011–2021 | GHS-POP | 0,13 | 0,62 (0,60–0,64) |
| 2 | 2011–2021 | GHS-BUILT-S | 0,14 | 0,59 (0,58–0,60) |
| 2 | 2011–2021 | geometrisk | 0,13 | 0,00 (0,00–0,00) |

