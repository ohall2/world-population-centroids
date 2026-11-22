# B – Jämförelsemått (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §6): hur står sig utgåvans ljuscentroid mot det en användare kan få i dag? Varje jämförelsekälla körs genom samma centroidkod som ljuset och jämförs mot facit med samma läges- och förändringsmått, på nivå 1 och 2, bara för facitländerna. Geometrisk mittpunkt är med som jämförelse.

**Källor:** ljus (utgåvan), **GHS-POP** och **GHS-BUILT-S** (JRC GHSL, R2023A, 30"), och **WorldPop** (Global unconstrained, 1 km). Alla med epokerna/åren 2000/2010/2011/2020; facitår knutet till närmaste år.

**Förbehåll.** GHS-POP bygger på folkräkning och GHS-BUILT-S på bebyggd yta, inte på nattljus; de är ändå inte oberoende av facit (GHS-POP fördelar folkräkningens folk), och förändringen mellan epoker är delvis interpolerad. WorldPop använder nattljus som förklaringsvariabel och är **inte oberoende av ljuset**; den redovisas därför i en egen sektion nedan, skild från de andra. Jämförelsen visar vad en användare kan få, inte ett oberoende test av nattljuset.

Alla tal är skrivna ur `docs/extended/<facit>_benchmark_*.csv` med `tools/write_comparison_doc.py`. Läge: median great circle-fel mot facit (km); förändring: de frysta måtten i `tools/change_detection.py`; den parvisa jämförelsen är per enhet mot ljuset (andel enheter där källan ligger närmare facit, och median av skillnaden ljus−källa, positiv = källan närmare). Allt med 95 %-blockbootstrap. Geometrisk mittpunkt är statisk, så dess förändringsskattning är alltid noll.

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

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2000 | GHS-POP | 1,00 (1,00–1,00) | 12,48 (9,79–17,53) |
| 1 | 2000 | GHS-BUILT-S | 0,46 (0,32–0,60) | −0,49 (−2,80–2,02) |
| 1 | 2000 | geometrisk | 0,04 (0,00–0,10) | −53,58 (−73,21–−22,68) |
| 1 | 2010 | GHS-POP | 0,98 (0,94–1,00) | 10,07 (6,67–16,84) |
| 1 | 2010 | GHS-BUILT-S | 0,40 (0,26–0,54) | −2,66 (−4,62–1,36) |
| 1 | 2010 | geometrisk | 0,06 (0,00–0,14) | −53,48 (−74,88–−24,39) |
| 1 | 2020 | GHS-POP | 0,92 (0,84–0,98) | 8,74 (5,78–13,71) |
| 1 | 2020 | GHS-BUILT-S | 0,46 (0,32–0,60) | −1,72 (−7,20–2,11) |
| 1 | 2020 | geometrisk | 0,06 (0,00–0,14) | −54,57 (−80,95–−24,55) |
| 2 | 2000 | GHS-POP | 0,89 (0,87–0,91) | 1,01 (0,93–1,11) |
| 2 | 2000 | GHS-BUILT-S | 0,61 (0,58–0,64) | 0,26 (0,20–0,34) |
| 2 | 2000 | geometrisk | 0,11 (0,09–0,13) | −3,55 (−4,38–−2,97) |
| 2 | 2010 | GHS-POP | 0,99 (0,98–0,99) | 1,38 (1,28–1,49) |
| 2 | 2010 | GHS-BUILT-S | 0,60 (0,56–0,64) | 0,26 (0,16–0,36) |
| 2 | 2010 | geometrisk | 0,12 (0,10–0,13) | −3,70 (−4,59–−3,12) |
| 2 | 2020 | GHS-POP | 0,94 (0,93–0,95) | 1,10 (1,02–1,18) |
| 2 | 2020 | GHS-BUILT-S | 0,57 (0,53–0,61) | 0,20 (0,12–0,29) |
| 2 | 2020 | geometrisk | 0,11 (0,09–0,13) | −3,95 (−4,68–−3,23) |

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

**Bortfall (den gemensamma mängden kräver värde från alla källor).**

| Nivå | År | facitenheter | använda | bortfall | källor med saknade värden |
| --- | --- | --- | --- | --- | --- |
| 1 | 2000 | 50 | 50 | 0 | – |
| 1 | 2010 | 50 | 50 | 0 | – |
| 1 | 2020 | 50 | 50 | 0 | – |
| 2 | 2000 | 2441 | 2441 | 0 | – |
| 2 | 2010 | 2400 | 2400 | 0 | – |
| 2 | 2020 | 2436 | 2436 | 0 | – |

## Mexiko

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2010 (n=32) | 8,88 (4,25–11,54) | 0,78 (0,44–1,12) | 8,10 (4,85–10,41) | 40,77 (26,57–48,31) |
| 1 | 2020 (n=32) | 7,38 (3,95–9,84) | 3,02 (1,68–4,34) | 8,64 (5,82–12,35) | 39,83 (27,00–48,58) |
| 2 | 2010 (n=1950) | 1,12 (0,98–1,27) | 0,35 (0,33–0,38) | 0,85 (0,73–1,03) | 3,19 (2,70–4,31) |
| 2 | 2020 (n=1947) | 1,03 (0,91–1,23) | 0,44 (0,40–0,51) | 0,87 (0,73–1,15) | 3,28 (2,74–4,43) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2010 | GHS-POP | 0,97 (0,91–1,00) | 7,52 (3,77–10,77) |
| 1 | 2010 | GHS-BUILT-S | 0,59 (0,41–0,75) | 0,39 (−1,78–4,61) |
| 1 | 2010 | geometrisk | 0,06 (0,00–0,16) | −26,28 (−36,46–−12,25) |
| 1 | 2020 | GHS-POP | 0,72 (0,56–0,88) | 3,89 (0,48–7,04) |
| 1 | 2020 | GHS-BUILT-S | 0,47 (0,31–0,66) | −0,10 (−2,27–6,26) |
| 1 | 2020 | geometrisk | 0,06 (0,00–0,16) | −26,27 (−33,65–−10,88) |
| 2 | 2010 | GHS-POP | 0,84 (0,81–0,86) | 0,64 (0,53–0,80) |
| 2 | 2010 | GHS-BUILT-S | 0,62 (0,58–0,66) | 0,20 (0,13–0,30) |
| 2 | 2010 | geometrisk | 0,19 (0,14–0,24) | −1,87 (−2,62–−1,44) |
| 2 | 2020 | GHS-POP | 0,80 (0,78–0,83) | 0,51 (0,44–0,64) |
| 2 | 2020 | GHS-BUILT-S | 0,60 (0,55–0,65) | 0,13 (0,07–0,19) |
| 2 | 2020 | geometrisk | 0,15 (0,11–0,18) | −1,94 (−2,79–−1,54) |

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

**Bortfall (den gemensamma mängden kräver värde från alla källor).**

| Nivå | År | facitenheter | använda | bortfall | källor med saknade värden |
| --- | --- | --- | --- | --- | --- |
| 1 | 2010 | 32 | 32 | 0 | – |
| 1 | 2020 | 32 | 32 | 0 | – |
| 2 | 2010 | 1952 | 1950 | 2 | ljus (utgåvan) (2) |
| 2 | 2020 | 1947 | 1947 | 0 | – |

## Brasilien

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2022 (n=27) | 16,52 (11,97–19,63) | 2,51 (0,94–4,60) | 22,23 (15,53–26,97) | 98,29 (75,69–141,35) |
| 2 | 2022 (n=5567) | 1,11 (0,99–1,32) | 0,34 (0,30–0,43) | 0,86 (0,78–1,02) | 3,82 (3,31–4,80) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2022 | GHS-POP | 0,96 (0,89–1,00) | 13,64 (9,28–17,16) |
| 1 | 2022 | GHS-BUILT-S | 0,30 (0,15–0,48) | −3,13 (−8,04–−0,13) |
| 1 | 2022 | geometrisk | 0,00 (0,00–0,00) | −74,15 (−126,39–−61,47) |
| 2 | 2022 | GHS-POP | 0,88 (0,87–0,90) | 0,70 (0,62–0,82) |
| 2 | 2022 | GHS-BUILT-S | 0,62 (0,60–0,65) | 0,20 (0,15–0,27) |
| 2 | 2022 | geometrisk | 0,12 (0,10–0,13) | −2,48 (−3,12–−2,07) |

**Bortfall (den gemensamma mängden kräver värde från alla källor).**

| Nivå | År | facitenheter | använda | bortfall | källor med saknade värden |
| --- | --- | --- | --- | --- | --- |
| 1 | 2022 | 27 | 27 | 0 | – |
| 2 | 2022 | 5567 | 5567 | 0 | – |

## Europa

**Läge: median great circle-fel mot facit (km).**

| Nivå | År | ljus (utgåvan) | GHS-POP | GHS-BUILT-S | geometrisk |
| --- | --- | --- | --- | --- | --- |
| 1 | 2011 (n=348) | 2,97 (2,54–3,62) | 0,20 (0,15–0,25) | 2,37 (1,81–3,09) | 8,67 (6,77–12,61) |
| 1 | 2021 (n=331) | 2,42 (2,05–3,25) | 0,44 (0,36–0,54) | 2,61 (2,22–3,50) | 9,08 (7,04–13,30) |
| 2 | 2011 (n=6266) | 0,95 (0,89–1,01) | 0,22 (0,21–0,24) | 0,46 (0,42–0,52) | 1,79 (1,64–1,99) |
| 2 | 2021 (n=6247) | 0,71 (0,66–0,76) | 0,27 (0,25–0,28) | 0,47 (0,43–0,53) | 1,81 (1,65–2,00) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2011 | GHS-POP | 0,97 (0,95–0,99) | 2,68 (2,30–3,34) |
| 1 | 2011 | GHS-BUILT-S | 0,61 (0,53–0,67) | 0,43 (0,17–0,69) |
| 1 | 2011 | geometrisk | 0,12 (0,08–0,17) | −4,98 (−8,00–−3,59) |
| 1 | 2021 | GHS-POP | 0,93 (0,91–0,96) | 1,97 (1,55–2,52) |
| 1 | 2021 | GHS-BUILT-S | 0,51 (0,44–0,59) | 0,02 (−0,33–0,26) |
| 1 | 2021 | geometrisk | 0,09 (0,06–0,14) | −5,80 (−8,90–−4,54) |
| 2 | 2011 | GHS-POP | 0,88 (0,87–0,90) | 0,62 (0,56–0,67) |
| 2 | 2011 | GHS-BUILT-S | 0,76 (0,73–0,78) | 0,38 (0,34–0,43) |
| 2 | 2011 | geometrisk | 0,27 (0,25–0,30) | −0,78 (−0,93–−0,67) |
| 2 | 2021 | GHS-POP | 0,81 (0,80–0,83) | 0,36 (0,33–0,41) |
| 2 | 2021 | GHS-BUILT-S | 0,68 (0,66–0,69) | 0,18 (0,17–0,20) |
| 2 | 2021 | geometrisk | 0,18 (0,17–0,20) | −0,99 (−1,13–−0,88) |

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

**Bortfall (den gemensamma mängden kräver värde från alla källor).**

| Nivå | År | facitenheter | använda | bortfall | källor med saknade värden |
| --- | --- | --- | --- | --- | --- |
| 1 | 2011 | 348 | 348 | 0 | – |
| 1 | 2021 | 331 | 331 | 0 | – |
| 2 | 2011 | 6298 | 6266 | 32 | ljus (utgåvan) (32) |
| 2 | 2021 | 6254 | 6247 | 7 | ljus (utgåvan) (7) |

## WorldPop (redovisas skilt: beroende av nattljus)

WorldPop använder nattljus som en förklaringsvariabel i fördelningen av befolkningen och är därför inte oberoende av det utgåvan mäter. Talen nedan ska läsas som en jämförelse av vad en användare kan få, inte som ett oberoende test.

### USA

**Läge (km).**

| Nivå | År | ljus (utgåvan) | WorldPop |
| --- | --- | --- | --- |
| 1 | 2000 (n=50) | 13,15 (10,23–17,62) | 0,34 (0,26–0,42) |
| 1 | 2010 (n=50) | 10,53 (6,92–16,93) | 0,16 (0,09–0,28) |
| 1 | 2020 (n=50) | 12,14 (8,19–16,58) | 1,91 (1,26–2,54) |
| 2 | 2000 (n=2441) | 1,51 (1,41–1,63) | 0,40 (0,37–0,45) |
| 2 | 2010 (n=2400) | 1,52 (1,43–1,62) | 0,08 (0,07–0,09) |
| 2 | 2020 (n=2436) | 1,46 (1,37–1,60) | 0,28 (0,26–0,31) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2000 | WorldPop | 1,00 (1,00–1,00) | 12,69 (9,78–17,49) |
| 1 | 2010 | WorldPop | 0,98 (0,94–1,00) | 10,08 (6,68–16,87) |
| 1 | 2020 | WorldPop | 0,96 (0,88–1,00) | 9,02 (5,74–14,09) |
| 2 | 2000 | WorldPop | 0,89 (0,87–0,91) | 1,01 (0,94–1,10) |
| 2 | 2010 | WorldPop | 0,98 (0,98–0,99) | 1,38 (1,29–1,50) |
| 2 | 2020 | WorldPop | 0,94 (0,93–0,95) | 1,10 (1,03–1,19) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2000–2010 | ljus (utgåvan) | 5,50 | 0,50 (0,36–0,64) |
| 1 | 2000–2010 | WorldPop | 0,16 | 1,00 (1,00–1,00) |
| 1 | 2000–2020 | ljus (utgåvan) | 6,56 | 0,58 (0,44–0,72) |
| 1 | 2000–2020 | WorldPop | 1,75 | 0,92 (0,84–0,98) |
| 1 | 2010–2020 | ljus (utgåvan) | 5,67 | 0,54 (0,40–0,68) |
| 1 | 2010–2020 | WorldPop | 1,77 | 0,74 (0,62–0,86) |
| 2 | 2000–2010 | ljus (utgåvan) | 0,74 | 0,61 (0,59–0,64) |
| 2 | 2000–2010 | WorldPop | 0,37 | 0,56 (0,54–0,58) |
| 2 | 2000–2020 | ljus (utgåvan) | 1,05 | 0,66 (0,64–0,69) |
| 2 | 2000–2020 | WorldPop | 0,54 | 0,57 (0,55–0,59) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,79 | 0,58 (0,56–0,60) |
| 2 | 2010–2020 | WorldPop | 0,26 | 0,54 (0,51–0,56) |

### Mexiko

**Läge (km).**

| Nivå | År | ljus (utgåvan) | WorldPop |
| --- | --- | --- | --- |
| 1 | 2010 (n=32) | 8,88 (4,25–11,54) | 0,74 (0,42–1,09) |
| 1 | 2020 (n=32) | 7,38 (3,95–9,84) | 2,05 (1,46–2,53) |
| 2 | 2010 (n=1950) | 1,12 (0,98–1,27) | 0,44 (0,41–0,49) |
| 2 | 2020 (n=1947) | 1,03 (0,91–1,23) | 0,53 (0,47–0,60) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2010 | WorldPop | 0,97 (0,91–1,00) | 7,54 (3,80–10,83) |
| 1 | 2020 | WorldPop | 0,81 (0,66–0,94) | 4,34 (1,06–8,26) |
| 2 | 2010 | WorldPop | 0,80 (0,77–0,83) | 0,56 (0,45–0,70) |
| 2 | 2020 | WorldPop | 0,76 (0,74–0,80) | 0,45 (0,35–0,58) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2010–2020 | ljus (utgåvan) | 4,74 | 0,72 (0,56–0,88) |
| 1 | 2010–2020 | WorldPop | 1,88 | 0,88 (0,75–0,97) |
| 2 | 2010–2020 | ljus (utgåvan) | 0,85 | 0,56 (0,53–0,58) |
| 2 | 2010–2020 | WorldPop | 0,32 | 0,52 (0,50–0,54) |

### Brasilien

**Läge (km).**

| Nivå | År | ljus (utgåvan) | WorldPop |
| --- | --- | --- | --- |
| 1 | 2022 (n=27) | 16,52 (11,97–19,63) | 2,50 (1,57–4,61) |
| 2 | 2022 (n=5567) | 1,11 (0,99–1,32) | 0,38 (0,32–0,47) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2022 | WorldPop | 0,96 (0,89–1,00) | 12,89 (8,30–17,16) |
| 2 | 2022 | WorldPop | 0,86 (0,85–0,88) | 0,67 (0,59–0,76) |

### Europa

**Läge (km).**

| Nivå | År | ljus (utgåvan) | WorldPop |
| --- | --- | --- | --- |
| 1 | 2011 (n=348) | 2,97 (2,54–3,62) | 0,26 (0,21–0,32) |
| 1 | 2021 (n=331) | 2,42 (2,05–3,25) | 0,53 (0,42–0,69) |
| 2 | 2011 (n=6266) | 0,95 (0,89–1,01) | 0,43 (0,39–0,46) |
| 2 | 2021 (n=6247) | 0,71 (0,66–0,76) | 0,44 (0,41–0,47) |

**Parvis mot ljuset, per enhet.**

| Nivå | År | källa | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | 2011 | WorldPop | 0,95 (0,91–0,97) | 2,71 (2,21–3,26) |
| 1 | 2021 | WorldPop | 0,89 (0,86–0,93) | 1,78 (1,42–2,58) |
| 2 | 2011 | WorldPop | 0,73 (0,71–0,76) | 0,38 (0,32–0,45) |
| 2 | 2021 | WorldPop | 0,66 (0,63–0,69) | 0,19 (0,15–0,24) |

**Förändring.**

| Nivå | Period | källa | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 1 | 2011–2021 | ljus (utgåvan) | 1,97 | 0,49 (0,43–0,53) |
| 1 | 2011–2021 | WorldPop | 0,41 | 0,73 (0,65–0,81) |
| 2 | 2011–2021 | ljus (utgåvan) | 0,64 | 0,51 (0,49–0,52) |
| 2 | 2011–2021 | WorldPop | 0,22 | 0,53 (0,51–0,55) |

