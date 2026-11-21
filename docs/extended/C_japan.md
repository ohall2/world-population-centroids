# C (andra facit) – Japan, bekräftande (planen §7)

Japan är det andra facit inom VIIRS-tiden, valt efter att Kanada hamnat på rad 3 (inkonklusivt). Frågan (§7): mäter nattljusets centroider förflyttning för en förändring helt inom VIIRS-tiden (2015–2020)? Facit: e-Stats rutnätsstatistik ur folkräkningen (räknade tal), 500 m-mesh, rutans mittpunkt viktad med rutans befolkning (se docs/DATA.md; den nationella totalen är certifierad mot folkräkningen båda åren). Serie: utgåvan v2.0.0 (VIIRS 2015 och 2020, GADM 4.1 JPN).

Måtten körs med de frysta funktionerna via samma drivrutin som Kanada (`tools/canada_confirm.py`, facit som parameter); de frysta skripten ändras inte. Att drivrutinen räknar som dem är prövat mot census nivå 2 2010–2020 (`tests/test_canada_confirm.py`). Urvalet är complete och minst tio rutpunkter båda åren. Block = nivå 1-enhet (prefektur). Alla tal är skrivna ur `docs/extended/C_japan_*.csv` med `tools/write_japan_doc.py`.

## Huvudmått: andel rätt riktning, nivå 2, 2015–2020

Blockbootstrap med block = prefektur (47 block), 1792 kommuner, och permutationstest.

| Mått | Värde |
| --- | --- |
| Andel rätt riktning | **0,46** |
| 95 %-intervall (blockbootstrap) | 0,44 – 0,49 |
| Intervallbredd | 0,05 |
| Antal enheter (nivå 2) | 1792 |
| Permutationstest, nollfördelningens medel | 0,49 (0,47–0,52) |
| Permutationstest, nedre svans p (andel ≤ observerad) | 0,0030 |

### Tolkning, bestämd i förväg (§7, ordagrant; rad 4 tillagd i efterhand)

| Utfall | Slutsats |
| --- | --- |
| Punktskattning ≥ 0,65 och intervallets nedre gräns > 0,50 | VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade perioderna (0,51–0,66). |
| Punktskattning 0,50–0,65 och intervallbredd ≤ 0,15 | Resultatet från DMSP-tiden gäller också VIIRS. |
| Intervallbredd > 0,15 eller färre än 150 enheter | Inkonklusivt. Ett andra facit inom VIIRS-tiden läggs till enligt samma plan, valt innan dess data öppnas. |
| Punktskattning < 0,50 och intervallbredd ≤ 0,15 *(tillagd i efterhand)* **←** | VIIRS mäter inte riktningen i detta facit. Ligger intervallets övre gräns under 0,50 redovisas utfallet som signifikant under slumpen, och orsaken som outredd. |

**Utfallet hamnar på rad 4.** Punktskattningen är 0,46, intervallet 0,44–0,49 (bredd 0,05), n = 1792. Slutsatsen enligt tabellen: *VIIRS mäter inte riktningen i detta facit. Ligger intervallets övre gräns under 0,50 redovisas utfallet som signifikant under slumpen, och orsaken som outredd.*

Intervallets övre gräns (0,49) ligger under 0,50, och permutationstestets nedre svans ger p = 0,0030 (andel av nollfördelningen ≤ den observerade, nollmedel 0,49). Utfallet redovisas därför som **signifikant under slumpen**. Orsaken är **outredd** – inga analyser görs här för att förklara den.

### Samlad slutsats för VIIRS-frågan

Kanada (2016–2021) blev **inkonklusivt** (rätt riktning 0,59, intervall 0,51–0,68, rad 3). Japan (2015–2020) ligger **under slumpen** (0,46, rad 4). **Inget av de två facit ger stöd för att VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade serierna (0,51–0,66).** Stoppregeln (§7) gäller: inget tredje facit läggs till.

## Bimått

**Förändring 2015–2020 (median |e−t| mot median |t|).**

| Nivå | n | median \|t\| | median \|e−t\| | kvot | rätt riktning |
| --- | --- | --- | --- | --- | --- |
| 1 | 46 | 0,36 | 0,78 | 2,15 | 0,59 (0,46–0,72) |
| 2 | 1792 | 0,06 | 0,21 | 3,81 | 0,46 (0,44–0,49) |

**Lägesfel 2020 (median great circle-fel mot facit, km).**

| Nivå | källa | n | median fel | 95 %-intervall |
| --- | --- | --- | --- | --- |
| 1 | ljus (utgåvan) | 47 | 2,29 | 1,65–2,58 |
| 1 | GHS-POP | 47 | 0,61 | 0,52–0,89 |
| 1 | GHS-BUILT-S | 47 | 3,66 | 2,90–4,56 |
| 1 | WorldPop | 47 | 0,77 | 0,57–0,99 |
| 1 | geometrisk | 47 | 17,96 | 13,22–20,85 |
| 2 | ljus (utgåvan) | 1797 | 0,67 | 0,57–0,77 |
| 2 | GHS-POP | 1797 | 0,30 | 0,26–0,34 |
| 2 | GHS-BUILT-S | 1797 | 0,47 | 0,40–0,57 |
| 2 | WorldPop | 1797 | 0,65 | 0,54–0,81 |
| 2 | geometrisk | 1797 | 2,44 | 1,85–3,21 |

Ljuset/utgåvan och alla jämförelsekällor använder exakt facitåren 2015 och 2020 (GHS-POP/GHS-BUILT-S epok 2015 och 2020, WorldPop 2015 och 2020) – ingen epokförskjutning.

**Parvis mot ljuset, per enhet (läge 2020).**

| Nivå | källa | n | andel närmare facit än ljuset | median skillnad ljus−källa (km) |
| --- | --- | --- | --- | --- |
| 1 | GHS-POP | 47 | 0,94 (0,85–1,00) | 1,35 (0,97–1,96) |
| 1 | GHS-BUILT-S | 47 | 0,26 (0,13–0,38) | −1,37 (−2,37–−0,49) |
| 1 | WorldPop | 47 | 0,94 (0,85–1,00) | 1,47 (1,00–1,92) |
| 1 | geometrisk | 47 | 0,00 (0,00–0,00) | −15,07 (−18,67–−11,28) |
| 2 | GHS-POP | 1797 | 0,81 (0,79–0,83) | 0,32 (0,27–0,37) |
| 2 | GHS-BUILT-S | 1797 | 0,68 (0,65–0,72) | 0,17 (0,13–0,20) |
| 2 | WorldPop | 1797 | 0,55 (0,50–0,60) | 0,04 (−0,00–0,06) |
| 2 | geometrisk | 1797 | 0,14 (0,10–0,17) | −1,53 (−2,21–−1,05) |

**Hybridskattning på prefekturnivå (A1, beskrivande).**

| År | skattning | median lägesfel (km) |
| --- | --- | --- |
| 2015 | ljus | 1,79 |
| 2015 | hybrid-ljus | 0,27 |
| 2015 | hybrid-geom | 2,67 |
| 2020 | ljus | 2,29 |
| 2020 | hybrid-ljus | 0,28 |
| 2020 | hybrid-geom | 2,76 |

*A3 (detektionsmodellen) redovisas inte, som för Kanada: σ är inte definierat för en period helt inom VIIRS och modellen föll på sin egen prövning (se `A3_detection.md`).*

## Förbehåll

- **2020 var ett pandemiår**, vilket kan ha påverkat folkräkningens läge och jämförbarheten 2015→2020.
- WorldPop använder nattljus som förklaringsvariabel och är inte oberoende av det utgåvan mäter; GHS-POP och GHS-BUILT-S bygger på folkräkning respektive bebyggelse.
- Omkring 0,5 % av befolkningen ligger i rutor vars mittpunkt faller utanför GADM-polygonerna (kustnära rutor).

