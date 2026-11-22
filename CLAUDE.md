# World Population Centroids v2

Befolkningscentroider ur nattljus (DMSP 1992–2012, VIIRS 2013–2025), prövade mot
folkräkningar. Översikt: README.md. Metod: docs/METHODS.md. Aktuellt arbete:
docs/ANALYSIS_PLAN.md och docs/review/granskning-2026-11-20.md.

## Miljö

    source ~/wpc-venv/bin/activate
    python -m pytest -q tests

Alla kommandon körs från repots rot. `data/` ingår inte i git.

## Regler

- Följ docs/ANALYSIS_PLAN.md. Går något inte att genomföra: stanna och fråga.
  Ändringar i planen förs in som en rad under *Avvikelser*, aldrig tyst.
- En analys per session och en commit per analys.
- Tester först: syntetiska data med känt svar, innan analysskriptet skrivs.
- Ändra inte `data/release/v2.0.0/`, `data/validation/`, `docs/validation/`,
  `docs/change_detection/` eller `docs/robustness/`. Nya resultat läggs i
  `docs/extended/`.
- Kanadas facit får inte hämtas eller öppnas innan skripten för A1–A3 och B är
  committade (planen, avsnitt 7).
- Tal i dokument och manus skrivs ur CSV-filerna med skript, inte för hand.
- Hitta aldrig på referenser, datakällor eller versioner. Kontrollera i källan
  och ange vad som inte är kontrollerat.
- Repot på GitHub är offentligt. `docs/manuscript/` och `docs/review/` ligger utanför
  git (se `.gitignore`). Committa aldrig manus, bilaga, granskningar, läsexemplar
  eller Zenodo-länken med nyckel, och kontrollera synligheten före varje push.
- Dokumentation och kodkommentarer på svenska, som i resten av repot. Manuset
  och utgåvans README är på engelska.
