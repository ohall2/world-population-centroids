# C – Kanada, bekräftande (planen §7)

Frågan (docs/ANALYSIS_PLAN.md §7): håller resultatet om nattljusets centroider mäter förflyttning också för en förändring helt inom VIIRS-tiden (2016–2021)? Facit: Statistics Canada Geographic Attribute File 2016 och 2021, dissemination area-punkter viktade med DA-befolkning (se docs/DATA.md). Serie: utgåvan v2.0.0 (VIIRS 2016 och 2021, GADM 4.1 CAN).

Måtten körs med de frysta funktionerna via drivrutinen `tools/canada_confirm.py` (de frysta skripten ändras inte). Att drivrutinen räknar som de frysta skripten är prövat mot census nivå 2 2010–2020 (`tests/test_canada_confirm.py`): n=2383, rätt riktning 0,5820, median |e−t| 0,7947 återges exakt. Urvalet är som för övriga facit: complete och minst tio DA-punkter båda åren. Alla tal är skrivna ur `docs/extended/C_canada_*.csv` med `tools/write_canada_doc.py`.

## Huvudmått: andel rätt riktning, nivå 2, 2016–2021

Blockbootstrap med block = provins (13 provinser), 287 census divisions, och permutationstest.

| Mått | Värde |
| --- | --- |
| Andel rätt riktning | **0,59** |
| 95 %-intervall (blockbootstrap) | 0,51 – 0,68 |
| Intervallbredd | 0,17 |
| Antal enheter (nivå 2) | 287 |
| Permutationstest, nollfördelningens medel | 0,50 (0,45–0,56) |
| Permutationstest, p (andel ≥ observerad) | 0,0020 |

### Tolkning, bestämd i förväg (§7, ordagrant)

| Utfall | Slutsats |
| --- | --- |
| Punktskattning ≥ 0,65 och intervallets nedre gräns > 0,50 | VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade perioderna (0,51–0,66). |
| Punktskattning 0,50–0,65 och intervallbredd ≤ 0,15 | Resultatet från DMSP-tiden gäller också VIIRS. |
| Intervallbredd > 0,15 eller färre än 150 enheter **←** | Inkonklusivt. Ett andra facit inom VIIRS-tiden läggs till enligt samma plan, valt innan dess data öppnas. |

**Utfallet hamnar på rad 3.** Punktskattningen är 0,59 och intervallet 0,51–0,68 (bredd 0,17); antalet enheter är 287. Intervallbredden överstiger 0,15, så utfallet är: *Inkonklusivt. Ett andra facit inom VIIRS-tiden läggs till enligt samma plan, valt innan dess data öppnas.*

Riktningen ligger ändå över slumpen (permutationstestets nollmedel 0,50, p = 0,0020), men intervallet är för brett för en bestämd slutsats enligt tabellen.

## Bimått

**Förändring 2016–2021 (median |e−t| mot median |t|).**

| Nivå | n | median \|t\| | median \|e−t\| | kvot | rätt riktning |
| --- | --- | --- | --- | --- | --- |
| 1 | 13 | 2,10 | 8,34 | 3,98 | 0,62 (0,31–0,85) |
| 2 | 287 | 0,33 | 1,04 | 3,18 | 0,59 (0,51–0,68) |

**Lägesfel 2021 (median great circle-fel mot facit, km).**

| Nivå | källa | n | median fel | 95 %-intervall |
| --- | --- | --- | --- | --- |
| 1 | ljus (utgåvan) | 13 | 35,77 | 17,87–64,28 |
| 1 | GHS-POP | 13 | 1,85 | 1,52–5,66 |
| 1 | GHS-BUILT-S | 13 | 26,38 | 15,99–46,41 |
| 1 | WorldPop | 13 | 3,53 | 1,91–3,96 |
| 1 | geometrisk | 13 | 292,06 | 272,45–563,44 |
| 2 | ljus (utgåvan) | 288 | 2,57 | 1,84–4,38 |
| 2 | GHS-POP | 288 | 0,61 | 0,47–0,91 |
| 2 | GHS-BUILT-S | 288 | 2,37 | 1,63–4,72 |
| 2 | WorldPop | 288 | 0,58 | 0,47–0,91 |
| 2 | geometrisk | 288 | 15,56 | 11,43–28,56 |

Ljuset och utgåvan använder exakt 2016 och 2021. Jämförelsekällornas epoker är de närmaste: GHS-POP och GHS-BUILT-S 2015 och 2020, WorldPop 2016 och 2020 – alltså perioderna 2015–2020 respektive 2016–2020 mot facits 2016–2021 (Avvikelser §6/§7).

**Parvis mot ljuset, per enhet (läge 2021).**

| Nivå | källa | n | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | GHS-POP | 13 | 1,00 (1,00–1,00) | 32,72 (16,14–45,98) |
| 1 | GHS-BUILT-S | 13 | 0,69 (0,38–0,92) | 1,88 (−0,58–7,52) |
| 1 | WorldPop | 13 | 1,00 (1,00–1,00) | 32,11 (17,24–49,79) |
| 1 | geometrisk | 13 | 0,08 (0,00–0,23) | −273,23 (−472,78–−216,37) |
| 2 | GHS-POP | 288 | 0,93 (0,90–0,96) | 1,76 (1,26–3,09) |
| 2 | GHS-BUILT-S | 288 | 0,57 (0,49–0,61) | 0,20 (−0,13–0,42) |
| 2 | WorldPop | 288 | 0,94 (0,90–0,96) | 1,72 (1,25–3,13) |
| 2 | geometrisk | 288 | 0,08 (0,04–0,15) | −11,50 (−21,54–−8,68) |

**Hybridskattning på provinsnivå (A1, beskrivande, 13 provinser).**

| År | skattning | median lägesfel (km) |
| --- | --- | --- |
| 2016 | ljus | 42,34 |
| 2016 | hybrid-ljus | 1,51 |
| 2016 | hybrid-geom | 14,75 |
| 2021 | ljus | 35,77 |
| 2021 | hybrid-ljus | 0,99 |
| 2021 | hybrid-geom | 15,01 |

**A3 (detektionsmodellen) – redovisas inte.** Den frysta A3-modellens brusmått σ är definierat för perioder inom DMSP och över skarven, men inte för en period helt inom VIIRS (2016–2021). A3-modellen klarade dessutom inte sin egen prövning (H3b, H3c; se `A3_detection.md`). A3-bimåttet redovisas därför inte utan en uttrycklig σ-definition för VIIRS-perioden.

## Förbehåll

- Omkring 0,9 % av Kanadas befolkning ligger i DA-representativa punkter som faller utanför GADM-polygonerna (kustnära punkter).
- WorldPop använder nattljus som förklaringsvariabel och är inte oberoende av det utgåvan mäter; redovisas tillsammans med de andra här men ska läsas med det i åtanke. GHS-POP och GHS-BUILT-S bygger på folkräkning respektive bebyggelse.
- Jämförelsekällornas epoker (2015/2016–2020) skiljer sig från facits 2016–2021; ljuset använder exakt 2016 och 2021.

