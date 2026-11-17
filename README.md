# World Population Centroids – v2 (utvecklingsversion)

Befolkningscentroider för världens administrativa enheter (nivå 0–2), beräknade från satellitobserverat nattljus.

- **v1** (Hall m.fl. 2019, *Scientific Data* 6:235, [doi:10.1038/s41597-019-0250-z](https://doi.org/10.1038/s41597-019-0250-z)) täcker 1992–2013 med DMSP-OLS och GADM 2.8. Originalfilerna ligger orörda i [`legacy/v1/`](legacy/v1/).
- **v2** räknas om från grunden med en öppen, testad pipeline och förlängs med VIIRS. Arbetet pågår på grenen `v2`.

## Status (17 november 2026)

Huvudserien är klar som utvecklingsversion: GADM 4.1 nivå 0–2, 1992–2025. Den är
validerad mot US Census (2000, 2010, 2020), Eurostat GEOSTAT (2011, 2021) och
INEGI:s ortsregister för Mexiko (2010, 2020). Nivå 0 beräknas utan ljus norr om
65° N, som i v1; nivå 1–2 med allt ljus.
Resonemang och resultat steg för steg: [`docs/HARMONIZATION.md`](docs/HARMONIZATION.md).

| Del | Status |
| --- | --- |
| DMSP-OLS 1992–2013, IRQR-kalibrering som i v1 | Klar |
| Samregistrering av DMSP-kompositerna mot VIIRS | Klar |
| VIIRS VNL v2 2013–2025, fackelmask 2 km | Klar |
| Skarvning per enhet (DMSP flyttad till VIIRS 2013) = huvudserien | Klar |
| Validering mot folkräkningar (USA, Europa) | Klar |
| Jämförelser: CCNL, Li m.fl. 2020, Chiovelli m.fl. 2026, VIIRS harmoniserad till DMSP | Klara |
| Oberoende validering i ett tredje område (Mexiko) | Klar |
| Maskning norr om 65° N på nivå 0, som i v1 | Klar |
| Publicering (Zenodo-DOI, licens, dataartikel) | Planerad |

## Snabbstart

Alla kommandon körs från projektmappen. Data hamnar i `data/` (ingår inte i git).

```bash
python -m venv ~/wpc-venv && source ~/wpc-venv/bin/activate
pip install -r requirements.txt

# 1. Indata (se docs/DATA.md)
python tools/download_ee.py --product dmsp --project <earth-engine-projekt>
python tools/download_ee.py --product viirs --collection NOAA/VIIRS/DNB/ANNUAL_V21 --bands average_masked --project <projekt>
python tools/download_ee.py --product viirs --bands average_masked --project <projekt>
python tools/download_zenodo.py --product flares-viirs
python tools/prepare_flares.py
python tools/prepare_flares_viirs.py --radii 2
#    GADM 4.1: ladda ner gadm_410-levels.zip till data/raw/gadm41/ och packa upp den (unzip)

# 2. Samregistrering av DMSP (kräver VDMSP_2013.tif från tools/harmonize_viirs.py)
python tools/harmonize_viirs.py --fit-only && python tools/harmonize_viirs.py --years 2013
python tools/coregister_dmsp.py

# 3. Centroider och huvudserien
python -m wpc.run --source dmsp-irqr-reg --flares data/processed/flares_viirs_r2km.gpkg --out data/results_mask2km
python -m wpc.run --source viirs
python tools/splice_series.py --direction backward --viirs viirs --viirs-results data/results --name spliced-back
python -m wpc.run --source dmsp-irqr-reg --max-lat 65 --flares data/processed/flares_viirs_r2km.gpkg --out data/results_m65
python -m wpc.run --source viirs --max-lat 65 --out data/results_m65
python tools/splice_series.py --direction backward --viirs viirs --results data/results_m65 --name spliced-back
python tools/build_release.py             # -> data/release/v<version>/ (nivå 0 ur results_m65, 1–2 ur results_mask2km)

# 4. Validering
python tools/download_validation.py
python tools/validation_centroids.py --dataset census --years 2000 2010 2020
python tools/validation_centroids.py --dataset geostat --years 2011 2021
python tools/validation_centroids.py --dataset inegi --years 2010 2020
python tools/validate_series.py --dataset census --years 2000 2010 2020 --change 2000 2010 --change 2010 2020 \
    --release data/release/v2.0.0.dev0 --out docs/validation   # likadant för geostat och inegi
python tools/build_release.py             # igen, så att metadata får valideringen

# 5. Tester
python -m pytest -q tests
```

Resultaten skrivs till `data/results*/<gränser>/<källa>/adm<N>/<år>.csv`, med en JSON-fil per år som beskriver körningen. Färdiga år hoppas över, så en avbruten körning kan startas om med samma kommando.

## Metod i korthet

Se [`docs/METHODS.md`](docs/METHODS.md) för detaljer och för skillnaderna mot v1.

1. DMSP-kompositerna kalibreras per satellit med koefficienterna i v1 ([`config/calibration_v1.csv`](config/calibration_v1.csv)), och medelvärdet tas per år. Gasfacklor maskas.
2. För varje enhet viktas varje pixel med **ljus × andel av pixeln i enheten × pixelns yta**.
3. Två centroider beräknas: Aboufadel–Austin (viktat medel av 3D-vektorer på sfären) och Barmore (iterativt viktat medel i azimutal ekvidistant projektion).
4. Varje centroid flaggas: 0 = i egen enhet, 1 = utanför alla enheter, 2 = i en annan enhet, -1 = inget ljus.

## Utdata (kolumner)

| Kolumn | Innehåll |
| --- | --- |
| `source`, `level`, `year` | Ljuskälla, nivå, år |
| `unit_id` + GADM-fält | Enhetens ID och namn |
| `aa_lon`, `aa_lat` | Centroid enligt Aboufadel–Austin (grader) |
| `bm_lon`, `bm_lat`, `bm_iter` | Centroid enligt Barmore, antal iterationer |
| `aa_bm_m` | Avstånd mellan metoderna (m) |
| `sum_light` | Σ ljus × andel |
| `sum_weight` | Σ ljus × andel × yta (km²) |
| `n_cells`, `n_lit` | Pixlar som enheten berör, varav med ljus |
| `area_km2` | Enhetens yta enligt pixelandelarna |
| `aa_flag`, `bm_flag` | Flaggor enligt ovan |

## Citera

Hall, O., Bustos, M. F. A., Olén, N. B. & Niedomysl, T. (2019). Population centroids of the world administrative units from nighttime lights 1992–2013. *Scientific Data* 6, 235.
