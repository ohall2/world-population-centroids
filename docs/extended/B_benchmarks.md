# B – Jämförelsemått (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §6): hur står sig utgåvans ljuscentroid mot det en användare kan få i dag? Varje jämförelsekälla körs genom samma centroidkod som ljuset och jämförs mot facit med samma läges- och förändringsmått, på nivå 1 och 2, bara för facitländerna. Geometrisk mittpunkt är med som jämförelse.

**Denna körning:** ljus (utgåvan) och **GHS-POP** (JRC GHSL, R2023A, 30", epokerna 2000/2010/2020; facitår knutet till närmaste epok). WorldPop och GHS-BUILT-S läggs till i senare körningar.

**Förbehåll.** GHS-POP bygger på folkräkning och bebyggelse, inte på nattljus; den är inte oberoende av facit, och förändringen mellan epoker är delvis interpolerad. Jämförelsen visar vad en användare kan få, inte ett oberoende test av nattljuset. WorldPop (senare) använder nattljus som förklaringsvariabel och redovisas skilt från GHS-POP och GHS-BUILT-S.

Alla tal är skrivna ur `docs/extended/<facit>_benchmark_*.csv` med `tools/write_comparison_doc.py`. Läge: median great circle-fel mot facit (km); förändring: de frysta måtten i `tools/change_detection.py`; båda med 95 %-blockbootstrap. Geometrisk mittpunkt är statisk, så dess förändringsskattning är alltid noll (rätt riktning 0).

## USA

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | geometrisk |
| --- | --- | --- | --- | --- |
| 1 | 2000 (n=50) | 13,15 (10,23–17,62) | 0,31 (0,23–0,46) | 62,23 (36,52–92,43) |
| 1 | 2010 (n=50) | 10,53 (7,13–17,84) | 0,18 (0,12–0,30) | 61,20 (35,49–96,49) |
| 1 | 2020 (n=50) | 12,14 (8,27–15,94) | 2,04 (1,52–2,46) | 63,39 (37,14–101,90) |
| 2 | 2000 (n=2441) | 1,51 (1,41–1,62) | 0,40 (0,37–0,44) | 5,40 (4,72–6,31) |
| 2 | 2010 (n=2400) | 1,52 (1,44–1,62) | 0,09 (0,08–0,10) | 5,50 (4,77–6,42) |
| 2 | 2020 (n=2436) | 1,46 (1,37–1,59) | 0,29 (0,26–0,32) | 5,58 (4,90–6,56) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2000–2010 | ljus (utgåvan) | 5,50 | 0,50 (0,36–0,64) |
| 1 | 2000–2010 | GHS-POP | 0,18 | 0,98 (0,92–1,00) |
| 1 | 2000–2010 | geometrisk | 2,74 | 0,00 (0,00–0,00) |
| 1 | 2000–2020 | ljus (utgåvan) | 6,56 | 0,58 (0,44–0,72) |
| 1 | 2000–2020 | GHS-POP | 1,77 | 0,88 (0,78–0,96) |
| 1 | 2000–2020 | geometrisk | 3,88 | 0,00 (0,00–0,00) |
| 1 | 2010–2020 | ljus (utgåvan) | 5,67 | 0,54 (0,40–0,68) |
| 1 | 2010–2020 | GHS-POP | 1,86 | 0,74 (0,62–0,84) |
| 1 | 2010–2020 | geometrisk | 2,60 | 0,00 (0,00–0,00) |
| 2 | 2000–2010 | ljus (utgåvan) | 0,74 | 0,61 (0,59–0,64) |
| 2 | 2000–2010 | GHS-POP | 0,37 | 0,57 (0,54–0,59) |
| 2 | 2000–2010 | geometrisk | 0,37 | 0,00 (0,00–0,00) |
| 2 | 2000–2020 | ljus (utgåvan) | 1,05 | 0,66 (0,64–0,68) |
| 2 | 2000–2020 | GHS-POP | 0,53 | 0,59 (0,57–0,61) |
| 2 | 2000–2020 | geometrisk | 0,54 | 0,00 (0,00–0,00) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,79 | 0,58 (0,56–0,60) |
| 2 | 2010–2020 | GHS-POP | 0,26 | 0,59 (0,57–0,62) |
| 2 | 2010–2020 | geometrisk | 0,26 | 0,00 (0,00–0,00) |

## Mexiko

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | geometrisk |
| --- | --- | --- | --- | --- |
| 1 | 2010 (n=32) | 8,88 (4,25–12,06) | 0,78 (0,47–1,12) | 40,77 (26,53–48,31) |
| 1 | 2020 (n=32) | 7,38 (3,62–10,57) | 3,02 (1,68–4,34) | 39,83 (27,30–48,58) |
| 2 | 2010 (n=1950) | 1,12 (0,97–1,26) | 0,35 (0,33–0,38) | 3,19 (2,72–4,43) |
| 2 | 2020 (n=1947) | 1,03 (0,91–1,28) | 0,44 (0,40–0,51) | 3,28 (2,73–4,45) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2010–2020 | ljus (utgåvan) | 4,74 | 0,72 (0,56–0,88) |
| 1 | 2010–2020 | GHS-POP | 3,22 | 0,41 (0,25–0,59) |
| 1 | 2010–2020 | geometrisk | 2,07 | 0,00 (0,00–0,00) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,85 | 0,56 (0,53–0,59) |
| 2 | 2010–2020 | GHS-POP | 0,31 | 0,51 (0,48–0,54) |
| 2 | 2010–2020 | geometrisk | 0,25 | 0,00 (0,00–0,00) |

## Brasilien

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | geometrisk |
| --- | --- | --- | --- | --- |
| 1 | 2022 (n=27) | 16,52 (11,97–19,63) | 2,51 (0,94–4,60) | 98,29 (75,69–141,35) |
| 2 | 2022 (n=5567) | 1,11 (0,99–1,34) | 0,34 (0,30–0,43) | 3,82 (3,30–4,69) |

## Europa

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | geometrisk |
| --- | --- | --- | --- | --- |
| 1 | 2011 (n=348) | 2,97 (2,53–3,62) | 0,20 (0,16–0,25) | 8,67 (6,77–12,67) |
| 1 | 2021 (n=331) | 2,42 (2,04–3,27) | 0,44 (0,36–0,54) | 9,08 (6,85–13,02) |
| 2 | 2011 (n=6266) | 0,95 (0,89–1,01) | 0,22 (0,21–0,24) | 1,79 (1,64–1,99) |
| 2 | 2021 (n=6247) | 0,71 (0,66–0,76) | 0,27 (0,25–0,28) | 1,81 (1,65–2,01) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2011–2021 | ljus (utgåvan) | 1,97 | 0,49 (0,43–0,53) |
| 1 | 2011–2021 | GHS-POP | 0,37 | 0,71 (0,64–0,79) |
| 1 | 2011–2021 | geometrisk | 0,42 | 0,00 (0,00–0,00) |
| 2 | 2011–2021 | ljus (utgåvan) | 0,64 | 0,51 (0,49–0,52) |
| 2 | 2011–2021 | GHS-POP | 0,13 | 0,62 (0,60–0,64) |
| 2 | 2011–2021 | geometrisk | 0,13 | 0,00 (0,00–0,00) |

