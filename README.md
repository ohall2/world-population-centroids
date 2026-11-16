# World Population Centroids – v2 (under arbete)

Befolkningscentroider för världens administrativa enheter (nivå 0–2), beräknade från satellitobserverat nattljus.

- **v1** (Hall m.fl. 2019, *Scientific Data* 6:235, [doi:10.1038/s41597-019-0250-z](https://doi.org/10.1038/s41597-019-0250-z)) täcker 1992–2013 med DMSP-OLS och GADM 2.8. Originalfilerna ligger orörda i [`legacy/v1/`](legacy/v1/).
- **v2** räknas om från grunden med en öppen, testad pipeline och förlängs med VIIRS. Arbetet pågår på grenen `v2`.

## Status

| Del | Status |
| --- | --- |
| DMSP-OLS 1992–2013, IRQR-kalibrering som i v1 | Kod klar och testad |
| Gränser GADM 4.1 (v2) och GADM 2.8 (jämförelse med v1) | Klart |
| CCNL 1992–2013 (Zhao m.fl. 2022) | Kod klar och testad, körning pågår |
| VIIRS 2012– | Nedladdning och källa klara, körning planerad |
| Validering mot SCB, Eurostat GISCO, US Census | Planerad |

## Snabbstart

Alla kommandon körs från projektmappen. Data hamnar i `data/` (ingår inte i git).

```bash
python -m venv ~/wpc-venv && source ~/wpc-venv/bin/activate
pip install -r requirements.txt

# 1. Indata (se docs/DATA.md)
python tools/download_ee.py --product dmsp --project <earth-engine-projekt>
python tools/download_ee.py --product ccnl --project <earth-engine-projekt>
python tools/prepare_flares.py
#    GADM 4.1: ladda ner gadm_410-levels.zip till data/raw/gadm41/ och packa upp den (unzip)

# 2. Centroider
python -m wpc.run --source dmsp-irqr      # alla år 1992–2013, adm0–2
python -m wpc.run --source ccnl

# 3. Tester
python -m pytest -q tests
```

Resultaten skrivs till `data/results/<gränser>/<ljuskälla>/adm<N>/<år>.csv`, med en JSON-fil per år som beskriver körningen. Färdiga år hoppas över, så en avbruten körning kan startas om med samma kommando.

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
