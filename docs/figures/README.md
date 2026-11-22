# Figurer till artikeln

Alla figurer ritas från filer i repot. PDF för manuset, PNG för förhandsvisning.
Skripten kräver matplotlib (`requirements.txt`) och körs från repots rot, t.ex.
`python tools/figures/fig6_sources.py`.

| Figur | Innehåll | Skript | Indata |
| --- | --- | --- | --- |
| Fig. 1 | Lägesfel per enhet efter nivå och facit (USA och Puerto Rico var för sig; Kanada 2021 och Japan 2020 tillagda) | tools/figures/fig1_location_error.py | docs/figures/data/unit_errors.csv.gz |
| Fig. 2 | Ljus per invånare per region och landsfelets riktning (MEX, BRA, SWE, NOR, FIN); pilar 3 × verklig längd | tools/figures/fig2_map.py | docs/check_identity/id_*.json, GADM nivå 1 (data/derived/check/adm1_simpl.gpkg) |
| Fig. 3 | Landsfelets uppdelning (A, Bb, Bw) och rutnätets betydelse för Mexiko | tools/figures/fig3_national_bias.py | docs/check_identity/id_*.json |
| Fig. 4 | Mätt mot sann förflyttning per kommun (nivå 2): USA, Mexiko, Europa, Kanada (±4 km) och Japan (±1 km) | tools/figures/fig4_measured_vs_true.py | docs/figures/data/unit_shifts.csv.gz |
| Fig. 5 | Andel rätt riktning efter den sanna förflyttningens storlek, fem facit | tools/figures/fig5_direction_by_shift.py | docs/change_detection/*_bins.csv (USA, Mexiko, Europa); docs/figures/data/unit_shifts.csv.gz (Kanada, Japan) |
| Fig. 6 | Jämförelsekällor på nivå 2, senaste folkräkningsår: (a) medianfel per källa och facit, (b) andel enheter där källan ligger närmare facit än ljuset | tools/figures/fig6_sources.py | docs/extended/*_benchmark_{position,pairs}.csv, docs/extended/C_{canada,japan}_{position,benchmark_pairs}.csv |
| Fig. 7 | Förändring i översikt: (a) andel rätt riktning nivå 2 per facit och period, ordnat efter sensorperiod; (b) nivå 1, ljus och hybridskattningarna | tools/figures/fig7_change_overview.py | docs/change_detection/*_change_detection.csv, docs/extended/*_hybrid_change.csv, docs/extended/C_{canada,japan}_change.csv |
| Fig. S1 | Årshopp i centroiden 2010–2017, fem serier, gemensamt urval (tidigare Fig. 6) | tools/figures/figS1_seam.py | docs/check_seam/*.csv.gz, utgåvans adm1-fil |

Mellanfilerna i docs/figures/data skapas av tools/figures/unit_errors.py och
tools/figures/unit_shifts.py från data/validation, data/validation_canada,
data/validation_japan och data/release/v2.0.0 (körs från repots rot). De återskapar
medianerna i docs/validation och docs/change_detection, och för Kanada och Japan
antal enheter, medianfel och andel rätt riktning i docs/extended/C_*.csv (prövas i
tests/test_figures_common.py). Mellanfilen är avrundad till 0,0001 km, så enstaka
medianer kan skilja i fjärde decimalen.

Gemensamma färger, markörer och hjälpfunktioner ligger i tools/figures/common.py.
USA (#2a78d6), Mexiko (#eb6834) och Europa (#1baf7a) har samma färger som tidigare.
Brasilien (#882255), Kanada (#ddaa33) och Japan (#3b1f6e) är nya; alla par av de
sex färgerna skiljer sig med ΔE ≥ 8,8 (OKLab × 100) under simulerad protanopi,
deuteranopi och tritanopi, och varje facit har dessutom en egen markör.

Att notera:

- Fig. 5: Kanada och Japan saknar klassfiler i docs/change_detection (fryst). Deras
  klasser räknas i skriptet ur mellanfilen med samma klassgränser
  (`common.direction_bins`), vilket för USA återger de frysta klasserna (test).
- Fig. 7b: hybridskattningarnas förändring finns inte för Kanada och Japan
  (tools/canada_confirm.py redovisar hybriden bara för läge), så där visas bara
  ljuset. Nivå 1-värdena för USA, Mexiko och Europa är tagna ur A1-filerna och gäller
  A1:s urval (för Europa 322 regioner mot 329 i Table 3).
- Fig. 7a: intervallen är blockbootstrap som i Table 3; Fig. 5 har binomiala
  Wilson-intervall.

Figurtexterna står i manuset, som inte ligger i repot.
