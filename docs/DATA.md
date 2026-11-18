# Indata

Alla filer läggs under `data/`, som inte ingår i git. Sökvägarna nedan är de som skripten använder som standard.

| Data | Sökväg | Hur den hämtas |
| --- | --- | --- |
| DMSP-OLS v4 Stable Lights, 35 kompositer | `data/raw/dmsp/DMSP_<satellit><år>_stable_lights.tif` + `.json` | `python tools/download_ee.py --product dmsp --project <earth-engine-projekt>` (kräver Earth Engine-konto och `earthengine authenticate`) |
| CCNL v1 (Zhao m.fl. 2022), 30″, en fil per år 1992–2013 | `data/raw/ccnl/CCNL_DMSP_<år>_V1.tif` + `.json` | `python tools/download_zenodo.py --product ccnl` (Zenodo doi:10.5281/zenodo.6644980, 22 filer, 7,3 GB; md5 kontrolleras mot Zenodo). Versionen i Earth Engine (`BNU/FGS/CCNL/v1`) är omsamplad till 0,008983° (1 km vid ekvatorn, 40 076 × 15 585 pixlar) och används inte. |
| VIIRS VNL v2, årskompositer 2013–2025 | `data/raw/viirs/VIIRS_<id>_average_masked.tif` + `.json` | `python tools/download_ee.py --product viirs --collection NOAA/VIIRS/DNB/ANNUAL_V21 --bands average_masked --project <projekt>` (2013–2021) och samma utan `--collection` (`ANNUAL_V22`, 2022–2025). Bara bandet `average_masked` används. |
| Fackelpolygoner (Elvidge m.fl. 2009) | `data/raw/flares/` → `data/processed/flares_dmsp.gpkg` | `python tools/prepare_flares.py` |
| Fackelkatalog MYC25Q4 (VIIRS Nightfire 2012–2025) | `data/raw/flares_viirs/annual_flare_site_records.csv` m.fl. | `python tools/download_zenodo.py --product flares-viirs` (Zenodo doi:10.5281/zenodo.23091822) |
| EOG:s årliga fackelkataloger 2012–2024 (Excel) | `data/raw/flares_viirs/legacy/*.xlsx` + `SHA256SUMS` | https://eogdata.mines.edu/global_flare_data/ (nio filer, länkade från https://eogdata.mines.edu/products/vnf/global_gas_flare.html) |
| Fackelmasker för VIIRS | `data/processed/flares_viirs_r<r>km.gpkg`, `flares_viirs_sites.csv`, `flares_viirs_report.json` | `python tools/prepare_flares_viirs.py --radii 0 1 2 3 5` |
| GADM 4.1, nivåer | `data/raw/gadm41/gadm_410-levels.gpkg` | Manuell nedladdning från https://geodata.ucdavis.edu/gadm/gadm4.1/gadm_410-levels.zip (2,5 GB), uppackad med `unzip` (5,0 GB). Zip-filen fungerar också, men är mycket långsam att läsa. |
| GADM 2.8, geometrier (för jämförelse med v1) | `data/raw/gadm28/gadm28_adm{0,1,2}.{shp,shx,prj,cpg}` | https://gml.noaa.gov/aftp/pub/basu/Borders/GADM/ |
| GADM 2.8, adm0-tabell | `data/raw/gadm28_levels.shp.zip` | Zenodo, doi:10.5281/zenodo.5103739 (md5 090efe30d3cb99f04a4814966832da6e) |
| Li m.fl. 2020 (jämförelse) | `data/raw/li2020/` | `python tools/download_ee.py --product li2020-dmsp --project <projekt>` och `--product li2020-viirs` |
| Chiovelli m.fl. 2026 (jämförelse) | `data/raw/chiovelli/<variant>/` | `python tools/download_zenodo.py --product chiovelli` (Zenodo doi:10.5281/zenodo.20558766) och `python tools/prepare_chiovelli.py` |
| Facit: US Census, befolkning och befolkningsmittpunkt per block group 2000, 2010, 2020 | `data/raw/validation/census/` | `python tools/download_validation.py --only census` |
| Facit: Eurostat GEOSTAT 1 km-rutnät 2011 (V2.0.1) och 2021 (V3) | `data/raw/validation/geostat/` | `python tools/download_validation.py --only geostat` |
| Facit: INEGI ortsregister (ITER) 2010 och 2020 | `data/raw/validation/inegi/` | `python tools/download_validation.py --only inegi` |
| Facit: IBGE Grade Estatística, Censo 2022, rutor om 1 km² (BR1KM_20251002.zip, 416 MB, sha256 944579f2…870a) | `data/raw/validation/ibge/2022/` | `curl -C - -o data/raw/validation/ibge/2022/BR1KM_20251002.zip https://geoftp.ibge.gov.br/recortes_para_fins_estatisticos/grade_estatistica/censo_2022/grade_1km/BR1KM_20251002.zip` (203 002 064 personer i 1 501 555 rutor) |
| PCNL (Li, S. m.fl. 2023), jämförelse | `data/raw/pcnl/PCNL<år>.tif` | `python tools/download_zenodo.py --record 7612389 --out data/raw/pcnl --pattern 'PCNL201[0-7].tif'` |
| GADM 2.8, färdig | `data/processed/gadm28.gpkg` | `python tools/prepare_gadm28.py` (använder `legacy/v1/` för attributen) |

## Licenser

- DMSP-OLS (NOAA/EOG): public domain.
- CCNL (Zhao m.fl. 2022, Scientific Data 9:424, doi:10.1038/s41597-022-01540-x; data doi:10.5281/zenodo.6644980): CC BY 4.0.
- VIIRS VNL v2 (Elvidge m.fl. 2021, EOG): CC BY 4.0.
- MYC25Q4 (Zhizhin m.fl. 2026): CC BY-NC 4.0, dvs. bara icke-kommersiellt bruk. Masken som byggs av katalogen omfattas.
- EOG:s årliga fackelkataloger: ingen licenstext på nedladdningssidan; EOG ska anges som källa.
- GADM: fritt för akademiskt och annat icke-kommersiellt bruk; vidarespridning och kommersiellt bruk kräver tillstånd.

## Om nedladdningen

`tools/download_ee.py` hämtar varje band i källans eget rutnät (inga pixlar
resamplas) och skriver bredvid varje GeoTIFF en JSON-fil med samling, bild-ID,
år, rutnät, Earth Engines datatyp och SHA-256. Pipelinen hittar årsprodukternas
filer via JSON-filernas år. `--dry-run` visar år, band, rutnät, datatyp och
okomprimerad storlek utan att hämta något. Varje anrop har en tidsgräns på två
minuter och görs om vid fel; en avbruten nedladdning startas om med samma
kommando.
