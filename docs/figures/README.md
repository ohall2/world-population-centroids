# Figurer till artikeln

Alla figurer ritas från filer i repot. PDF för manuset, PNG för förhandsvisning.

| Figur | Innehåll | Skript | Indata |
| --- | --- | --- | --- |
| Fig. 1 | Lägesfel per enhet efter nivå och facit (USA och Puerto Rico var för sig) | tools/figures/fig1_location_error.py | docs/figures/data/unit_errors.csv.gz |
| Fig. 2 | Ljus per invånare per region och landsfelets riktning (MEX, BRA, SWE, NOR, FIN); pilar 3 × verklig längd | tools/figures/fig2_map.py | docs/check_identity/id_*.json, GADM nivå 1 (data/derived/check/adm1_simpl.gpkg) |
| Fig. 3 | Landsfelets uppdelning (A, Bb, Bw) och rutnätets betydelse för Mexiko | tools/figures/fig3_national_bias.py | docs/check_identity/id_*.json |
| Fig. 4 | Mätt mot sann förflyttning per kommun (nivå 2) | tools/figures/fig4_measured_vs_true.py | docs/figures/data/unit_shifts.csv.gz |
| Fig. 5 | Andel rätt riktning efter den sanna förflyttningens storlek | tools/figures/fig5_direction_by_shift.py | docs/change_detection/*_bins.csv |
| Fig. 6 | Årshopp i centroiden 2010–2017, fem serier, gemensamt urval | tools/figures/fig6_seam.py | docs/check_seam/*.csv.gz, utgåvans adm1-fil |

Mellanfilerna i docs/figures/data skapas av tools/figures/unit_errors.py och
tools/figures/unit_shifts.py från data/validation och data/release/v2.0.0
(körs från repots rot). De återskapar exakt medianerna i docs/validation och
docs/change_detection.
