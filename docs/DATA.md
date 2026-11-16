# Indata

Alla filer läggs under `data/`, som inte ingår i git. Sökvägarna nedan är de som skripten använder som standard.

| Data | Sökväg | Hur den hämtas |
| --- | --- | --- |
| DMSP-OLS v4 Stable Lights, 35 kompositer | `data/raw/dmsp/DMSP_<satellit><år>_stable_lights.tif` + `.json` | `python tools/download_ee.py --product dmsp --project <earth-engine-projekt>` (kräver Earth Engine-konto och `earthengine authenticate`) |
| CCNL v1 (Zhao m.fl. 2022), en fil per år 1992–2013 | `data/raw/ccnl/CCNL_<id>_b1.tif` + `.json` | `python tools/download_ee.py --product ccnl --project <earth-engine-projekt>` (Earth Engine `BNU/FGS/CCNL/v1`) |
| VIIRS VNL v2, årskompositer | `data/raw/viirs/VIIRS_<id>_{average_masked,cf_cvg}.tif` + `.json` | `python tools/download_ee.py --product viirs --project <earth-engine-projekt>` (Earth Engine `NOAA/VIIRS/DNB/ANNUAL_V22`) |
| Fackelpolygoner (Elvidge m.fl. 2009) | `data/raw/flares/` → `data/processed/flares_dmsp.gpkg` | `python tools/prepare_flares.py` |
| GADM 4.1, nivåer | `data/raw/gadm41/gadm_410-levels.gpkg` | Manuell nedladdning från https://geodata.ucdavis.edu/gadm/gadm4.1/gadm_410-levels.zip (2,5 GB), uppackad med `unzip` (5,0 GB). Zip-filen fungerar också, men är mycket långsam att läsa. |
| GADM 2.8, geometrier (för jämförelse med v1) | `data/raw/gadm28/gadm28_adm{0,1,2}.{shp,shx,prj,cpg}` | https://gml.noaa.gov/aftp/pub/basu/Borders/GADM/ |
| GADM 2.8, adm0-tabell | `data/raw/gadm28_levels.shp.zip` | Zenodo, doi:10.5281/zenodo.5103739 (md5 090efe30d3cb99f04a4814966832da6e) |
| GADM 2.8, färdig | `data/processed/gadm28.gpkg` | `python tools/prepare_gadm28.py` (använder `legacy/v1/` för attributen) |

## Licenser

- DMSP-OLS (NOAA/EOG): public domain.
- CCNL (Zhao m.fl. 2022, doi:10.5194/essd-14-517-2022): CC BY 4.0.
- VIIRS VNL v2 (Elvidge m.fl. 2021, EOG): CC BY 4.0.
- GADM: fritt för akademiskt och annat icke-kommersiellt bruk; vidarespridning och kommersiellt bruk kräver tillstånd.

## Om nedladdningen

`tools/download_ee.py` hämtar varje band i källans eget rutnät (inga pixlar
resamplas) och skriver bredvid varje GeoTIFF en JSON-fil med samling, bild-ID,
år, rutnät, Earth Engines datatyp och SHA-256. Pipelinen hittar årsprodukternas
filer via JSON-filernas år. `--dry-run` visar år, band, rutnät, datatyp och
okomprimerad storlek utan att hämta något. Varje anrop har en tidsgräns på två
minuter och görs om vid fel; en avbruten nedladdning startas om med samma
kommando.
