# Indata

Alla filer läggs under `data/`, som inte ingår i git. Sökvägarna nedan är de som skripten använder som standard.

| Data | Sökväg | Hur den hämtas |
| --- | --- | --- |
| DMSP-OLS v4 Stable Lights, 35 kompositer | `data/raw/dmsp/DMSP_<satellit><år>_stable_lights.tif` + `.json` | `python tools/download_ee.py --project <earth-engine-projekt> --years all` (kräver Earth Engine-konto och `earthengine authenticate`) |
| Fackelpolygoner (Elvidge m.fl. 2009) | `data/raw/flares/` → `data/processed/flares_dmsp.gpkg` | `python tools/prepare_flares.py` |
| GADM 4.1, nivåer | `data/raw/gadm41/gadm_410-levels.gpkg` | Manuell nedladdning från https://geodata.ucdavis.edu/gadm/gadm4.1/gadm_410-levels.zip (2,5 GB), uppackad med `unzip` (5,0 GB). Zip-filen fungerar också, men är mycket långsam att läsa. |
| GADM 2.8, geometrier (för jämförelse med v1) | `data/raw/gadm28/gadm28_adm{0,1,2}.{shp,shx,prj,cpg}` | https://gml.noaa.gov/aftp/pub/basu/Borders/GADM/ |
| GADM 2.8, adm0-tabell | `data/raw/gadm28_levels.shp.zip` | Zenodo, doi:10.5281/zenodo.5103739 (md5 090efe30d3cb99f04a4814966832da6e) |
| GADM 2.8, färdig | `data/processed/gadm28.gpkg` | `python tools/prepare_gadm28.py` (använder `legacy/v1/` för attributen) |

## Licenser

- DMSP-OLS (NOAA/EOG): public domain.
- GADM: fritt för akademiskt och annat icke-kommersiellt bruk; vidarespridning och kommersiellt bruk kräver tillstånd.
