# A2 – Trendskattning (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §4): minskar felet om hela tidsserien används i stället för två ändpunkter?

Förflyttningen skattas per enhet och koordinat ur en Theil–Sen-lutning inom varje sensorperiod (DMSP t.o.m. 2012, VIIRS fr.o.m. 2013), i ett lokalt plan kring enhetens ljuscentroid 2013. Minst fyra år krävs i varje del. Jämförs med ändpunkter (nuvarande, ettårsfönster) och treårsmedel, på samma enheter. Ändpunkts- och treårsskattningarna och måtten är de frysta i `tools/change_detection.py` och `tools/robustness_checks.py`.

Alla tal är skrivna ur `docs/extended/<facit>_trend_change.csv` med `tools/write_trend_doc.py`. Lägesändringsmått: median great circle-fel mot facit |e−t| (km) och andel rätt riktning, med 95 %-blockbootstrap.

## USA

**Nivå 1.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2000–2010 | ändpunkter | 50 | 5,50 (4,50–8,35) | 0,50 (0,36–0,64) |
| 2000–2010 | treårsmedel | 50 | 5,17 (2,88–6,19) | 0,48 (0,36–0,62) |
| 2000–2010 | trend | 50 | 3,83 (2,94–5,05) | 0,62 (0,50–0,76) |
| 2000–2020 | ändpunkter | 50 | 6,56 (4,32–8,46) | 0,58 (0,44–0,72) |
| 2000–2020 | treårsmedel | 50 | 4,71 (3,54–8,84) | 0,64 (0,50–0,76) |
| 2000–2020 | trend | 50 | 4,70 (3,88–5,45) | 0,70 (0,58–0,82) |
| 2010–2020 | ändpunkter | 50 | 5,67 (4,49–7,32) | 0,54 (0,40–0,68) |
| 2010–2020 | treårsmedel | 50 | 4,31 (3,37–6,12) | 0,56 (0,42–0,70) |
| 2010–2020 | trend | 50 | 4,40 (2,88–6,44) | 0,68 (0,54–0,80) |

**Nivå 2.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2000–2010 | ändpunkter | 2381 | 0,74 (0,68–0,83) | 0,61 (0,59–0,64) |
| 2000–2010 | treårsmedel | 2381 | 0,63 (0,58–0,69) | 0,63 (0,61–0,66) |
| 2000–2010 | trend | 2381 | 0,57 (0,51–0,66) | 0,63 (0,60–0,65) |
| 2000–2020 | ändpunkter | 2373 | 1,05 (0,96–1,13) | 0,66 (0,64–0,68) |
| 2000–2020 | treårsmedel | 2373 | 0,93 (0,87–1,01) | 0,67 (0,65–0,69) |
| 2000–2020 | trend | 2373 | 0,87 (0,80–0,97) | 0,67 (0,64–0,69) |
| 2010–2020 | ändpunkter | 2383 | 0,79 (0,73–0,89) | 0,58 (0,56–0,60) |
| 2010–2020 | treårsmedel | 2383 | 0,71 (0,65–0,79) | 0,58 (0,56–0,60) |
| 2010–2020 | trend | 2383 | 0,73 (0,66–0,80) | 0,58 (0,56–0,61) |

## Mexiko

**Nivå 1.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2010–2020 | ändpunkter | 32 | 4,74 (2,10–7,55) | 0,72 (0,56–0,88) |
| 2010–2020 | treårsmedel | 32 | 5,21 (2,16–7,32) | 0,72 (0,56–0,88) |
| 2010–2020 | trend | 32 | 4,56 (3,85–6,68) | 0,62 (0,44–0,78) |

**Nivå 2.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2010–2020 | ändpunkter | 1903 | 0,85 (0,71–1,00) | 0,56 (0,54–0,59) |
| 2010–2020 | treårsmedel | 1903 | 0,76 (0,64–0,92) | 0,55 (0,52–0,58) |
| 2010–2020 | trend | 1903 | 0,74 (0,61–0,90) | 0,51 (0,48–0,53) |

## Europa

**Nivå 1.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2011–2021 | ändpunkter | 329 | 1,97 (1,57–2,52) | 0,49 (0,43–0,53) |
| 2011–2021 | treårsmedel | 329 | 1,65 (1,35–2,10) | 0,48 (0,42–0,52) |
| 2011–2021 | trend | 329 | 1,97 (1,47–2,51) | 0,41 (0,34–0,47) |

**Nivå 2.**

| Period | skattning | n | median \|e−t\| | rätt riktning |
| --- | --- | --- | --- | --- |
| 2011–2021 | ändpunkter | 5997 | 0,62 (0,59–0,66) | 0,51 (0,49–0,52) |
| 2011–2021 | treårsmedel | 5997 | 0,59 (0,55–0,63) | 0,51 (0,49–0,52) |
| 2011–2021 | trend | 5997 | 0,65 (0,60–0,71) | 0,46 (0,45–0,47) |

## Hypotes H2

H2: trendskattningen minskar median |e−t| på nivå 2 med minst 20 % mot ändpunkter i minst fyra av de fem perioderna (treårsmedel gav 6–15 %).

| Period | ändpunkter | treårsmedel (minskning) | trend (minskning) |
| --- | --- | --- | --- |
| USA 2000–2010 | 0,74 | 0,63 (15 %) | 0,57 (22 %) |
| USA 2000–2020 | 1,05 | 0,93 (12 %) | 0,87 (17 %) |
| USA 2010–2020 | 0,79 | 0,71 (10 %) | 0,73 (8 %) |
| Europa 2011–2021 | 0,62 | 0,59 (5 %) | 0,65 (−4 %) |
| Mexiko 2010–2020 | 0,85 | 0,76 (11 %) | 0,74 (13 %) |

Trenden minskar felet med minst 20 % i 1 av 5 perioder. **H2 stöds inte.** Treårsmedlets minskningar (5 %–15 %) ligger i linje med det tidigare resultatet.

