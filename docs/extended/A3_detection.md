# A3 – Detektionsmodell (utforskande)

Frågan (docs/ANALYSIS_PLAN.md §5): går det att ange per enhet hur stor en förflyttning måste vara för att synas?

Brusmåttet σ per enhet tas ur utgåvan: `noise_dmsp_km` inom DMSP, √(`noise_dmsp_km`² + `noise_viirs_km`²) över skarven. e är ändpunktsskattningen (utgåvans ljuscentroider y1→y2); |t| den verkliga förflyttningen i facit. Allt på nivå 2. Måtten och facitläsningen är de frysta i `tools/change_detection.py`; logistisk regression och AUC görs i numpy. Alla tal är skrivna ur `docs/extended/detection_*.csv` med `tools/write_detection_doc.py`.

## Steg 1: σ mot fel (Spearmans ρ, nivå 2)

| Period | n | ρ(σ, \|e−t\|) |
| --- | --- | --- |
| USA 2000–2010 | 2381 | 0,42 (0,35–0,48) |
| USA 2010–2020 | 2383 | 0,57 (0,53–0,61) |
| USA 2000–2020 | 2373 | 0,49 (0,45–0,54) |
| Mexiko 2010–2020 | 1915 | 0,65 (0,59–0,69) |
| Europa 2011–2021 | 6161 | 0,65 (0,62–0,67) |

**H3a** (ρ > 0 i alla facit): stöds. Större brus hör ihop med större fel i alla facit och perioder.

## Steg 2: logistisk modell, anpassad på USA counties 2000–2020

Rätt riktning ~ β0 + β1 · prediktor.

| Modell | prediktor | n | β0 | β1 |
| --- | --- | --- | --- | --- |
| snr | log(\|t\|/σ) | 2373 | 0,699 | 0,379 |
| abs_t | log\|t\| | 2373 | 0,926 | 0,356 |

## Steg 3: prövning (AUC och Brier)

| Prövningsfacit | modell | n | basandel | AUC | Brier | Brier (konstant) |
| --- | --- | --- | --- | --- | --- | --- |
| Mexiko 2010–2020 | log(\|t\|/σ) | 1915 | 0,56 | 0,511 | 0,256 | 0,247 |
| Mexiko 2010–2020 | log\|t\| | 1915 | 0,56 | 0,533 | 0,252 | 0,247 |
| USA 2010–2020 | log(\|t\|/σ) | 2383 | 0,58 | 0,579 | 0,239 | 0,243 |
| USA 2010–2020 | log\|t\| | 2383 | 0,58 | 0,566 | 0,241 | 0,243 |

**H3b** (σ-modellen högre AUC än \|t\|-modellen på Mexiko): stöds inte. σ-modellens AUC 0,511 mot 0,533 för \|t\|-modellen.

### Kalibrering på Mexiko (σ-modellen, tio lika stora klasser)

| Klass | n | förutsagd | observerad | skillnad |
| --- | --- | --- | --- | --- |
| 1 | 192 | 0,43 | 0,49 | 0,06 |
| 2 | 192 | 0,51 | 0,52 | 0,01 |
| 3 | 192 | 0,55 | 0,58 | 0,03 |
| 4 | 192 | 0,58 | 0,57 | 0,00 |
| 5 | 192 | 0,60 | 0,64 | 0,03 |
| 6 | 191 | 0,62 | 0,54 | 0,08 |
| 7 | 191 | 0,65 | 0,61 | 0,04 |
| 8 | 191 | 0,67 | 0,53 | 0,14 |
| 9 | 191 | 0,71 | 0,56 | 0,15 |
| 10 | 191 | 0,78 | 0,54 | 0,24 |

**H3c** (förutsagd och observerad skiljer ≤ 0,05 i minst åtta av tio klasser på Mexiko): stöds inte. 5 av tio klasser inom 0,05.

## Steg 4: minsta upptäckbara förflyttning (km)

Det |t| där σ-modellen (anpassad på USA 2000–2020) ger 0,8 sannolikhet för rätt riktning, per enhet. Fördelning över enheterna i varje facit.

| Period | n | σ median | min. flytt median | q25 | q75 |
| --- | --- | --- | --- | --- | --- |
| USA 2000–2010 | 2381 | 0,44 | 2,72 | 1,61 | 4,30 |
| USA 2010–2020 | 2383 | 0,53 | 3,22 | 2,06 | 5,00 |
| USA 2000–2020 | 2373 | 0,53 | 3,22 | 2,06 | 5,00 |
| Mexiko 2010–2020 | 1915 | 0,50 | 3,04 | 1,55 | 5,61 |
| Europa 2011–2021 | 6161 | 0,46 | 2,84 | 1,40 | 4,86 |

## Hypoteser

| Hypotes | Utfall |
| --- | --- |
| H3a: ρ > 0 i alla facit | stöds |
| H3b: σ-modellen högre AUC på Mexiko | stöds inte |
| H3c: kalibrering ≤ 0,05 i ≥ 8/10 klasser på Mexiko | stöds inte |

