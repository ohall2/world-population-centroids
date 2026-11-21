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
| GHS-POP R2023A, befolkning per cell, WGS84 30″ (epoker 2000, 2010, 2020), jämförelse (analys B) | `data/raw/ghspop/GHS_POP_E<epok>_GLOBE_R2023A_4326_30ss_V1_0.tif` + `.json` | `python tools/download_ghsl.py --product pop --epochs 2000 2010 2020` (JRC GHSL, https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/ ; tre zip om ca 447/454/460 MB; GeoTIFF:en packas upp och SHA-256 skrivs i JSON) |
| GHS-BUILT-S R2023A, bebyggd yta (m²) per cell, WGS84 30″ (epoker 2000, 2010, 2020), jämförelse (analys B) | `data/raw/ghsbuilts/GHS_BUILT_S_E<epok>_GLOBE_R2023A_4326_30ss_V1_0.tif` + `.json` | `python tools/download_ghsl.py --product built_s --epochs 2000 2010 2020` (JRC GHSL, https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/ ; tre zip om ca 171/175/178 MB; GeoTIFF:en packas upp och SHA-256 skrivs i JSON) |
| WorldPop Global 2000–2020, unconstrained, 1 km, befolkning per cell (åren 2000, 2010, 2011, 2020), jämförelse (analys B) | `data/raw/worldpop/ppp_<år>_1km_Aggregated.tif` + `.json` | `python tools/download_worldpop.py --years 2000 2010 2011 2020` (WorldPop, https://data.worldpop.org/GIS/Population/Global_2000_2020/ ; fyra GeoTIFF om ca 1,1/1,1/1,1/0,87 GB; SHA-256 skrivs i JSON) |
| GADM 2.8, färdig | `data/processed/gadm28.gpkg` | `python tools/prepare_gadm28.py` (använder `legacy/v1/` för attributen) |

## Licenser

- DMSP-OLS (NOAA/EOG): public domain.
- CCNL (Zhao m.fl. 2022, Scientific Data 9:424, doi:10.1038/s41597-022-01540-x; data doi:10.5281/zenodo.6644980): CC BY 4.0.
- VIIRS VNL v2 (Elvidge m.fl. 2021, EOG): CC BY 4.0.
- MYC25Q4 (Zhizhin m.fl. 2026): CC BY-NC 4.0, dvs. bara icke-kommersiellt bruk. Masken som byggs av katalogen omfattas.
- EOG:s årliga fackelkataloger: ingen licenstext på nedladdningssidan; EOG ska anges som källa.
- GADM: fritt för akademiskt och annat icke-kommersiellt bruk; vidarespridning och kommersiellt bruk kräver tillstånd.
- GHS-POP och GHS-BUILT-S (JRC GHSL, release R2023A; referens Pesaresi m.fl. 2024, *Int. J. Digital Earth* 17(1)): CC BY 4.0 (bekräftad på GHSL:s FAQ-sida). Dataseten hämtas med `tools/download_ghsl.py`. Datasetens exakta DOI är inte kontrollerade.
- WorldPop (www.worldpop.org, Global 2000–2020, unconstrained, 1 km): CC BY 4.0. WorldPop använder nattljus som en förklaringsvariabel och är därför inte oberoende av det utgåvan mäter; redovisas skilt i analys B. Licensvillkoren är inte kontrollerade i detalj mot WorldPops sida.

## Kanadas facit (Statistics Canada Geographic Attribute File)

Statistics Canada Geographic Attribute File (GAF), katalog 92-151-X, folkräkning
2016 och 2021. Filen har en rad per dissemination block (DB) med befolknings- och
bostadstal per DB och, på dissemination area-nivå (DA), den representativa
punkten `DArplat`/`DArplong` (latitud/longitud i decimalgrader). Ingen
DB-representativ punkt finns i filen något år.

**DA-punktens definition (båda åren, oförändrad metod).** DA-punkten är ett
befolkningsviktat medelcentrum av DA:ns dissemination blocks representativa
punkter, viktat med blockens befolkning. Faller den viktade punkten utanför DA:n
eller i vatten flyttas den med minsta befolkningsviktade kvadratavstånd så att
den ligger inom DA:n. Metoden är densamma 2016 och 2021; varje års fil väger med
det årets folkräkningsbefolkning (vilket specifikt år anges inte ordagrant i
referensguiden, men följer av att varje fil är en produkt för respektive
folkräkning). Källa: StatCans referensguide 92-151-G (2016 och 2021) och
folkräkningsordlistan; sidorna på statcan.gc.ca kunde inte hämtas direkt (403)
utan lästes via sökträffar på samma sidor.

**Facit som byggs här.** DB-raderna läggs ihop till DA: en punkt per DA
(`DArplat`/`DArplong`) med DA-befolkning = summan av dess blocks befolkning.
DA-punkterna viktas med DA-befolkningen in i GADM-enheterna (nivå 1 = provins,
nivå 2 = census division), som för övriga facit. Samma nivå och punkttyp båda
åren (planens §7, reservregeln: blockpunkter saknas i filen båda åren).

Licens: Open Government Licence – Canada (`ca-ogl-lgo`).

## Japans facit (e-Stat, folkräkningens rutnätsstatistik)

Statistics Bureau of Japan / e-Stat, rutnätsstatistik (Grid Square Statistics) ur
folkräkningen (räknade tal, ej modellerade), 2015 och 2020. Rutstorlek: 500 m
(halvmesh) om den finns med räknade tal i samma form båda åren, annars 1 km
(3:e ordningens mesh). Samma JIS-meshsystem båda åren. Rutans mittpunkt används
som punkt (härledd ur meshkoden `KEY_CODE`, `tools/validation_centroids.py:
mesh_centroids`), rutans befolkning som vikt. Läsaren hittar befolkningskolumnen
via den japanska etiketten (人口…総数) i CSV:ns andra rubrikrad (cp932).

**Dataset (statsId) och hämtväg.** 500 m-mesh (halvmesh) används båda åren:
**2015 = `T000847`**, **2020 = `T001101`** (e-Stats definitionsfiler
`https://www.e-stat.go.jp/help/data-definition-information/shuroku/T000847.pdf`
och `.../T001101.pdf`). Hämtas per prefektur (47 filer per år) med
`tools/download_estat_jp.py`, som anropar e-Stats statistiska GIS:
`https://www.e-stat.go.jp/gis/statmap-search/data?statsId=<ID>&code=<01..47>&downloadType=2&aggregateUnit=S`
(zip med en cp932-CSV per prefektur). En ruta som ligger på en prefekturgräns
delas mellan prefekturfilerna (varje fil bär sin del av rutans befolkning); därför
summerar läsaren rader med samma rutkod. Då certifieras nationella totalen exakt
mot folkräkningen: 2015 = 127 094 745 och 2020 = 126 146 099 (kontrollerat).

**Sekretess (markering och hantering, lika båda åren).** Vi använder **total
befolkning** (人口（総数）). Enligt e-Stats definitionsfil `T000880` döljs inte
den totala befolkningen per ruta i sekretessbehandlingen (bara köns- och
åldersfördelade tal kan döljas i små rutor); totalbefolkningens kolumn har alltså
inga hemlighållna värden, och den nationella summan bevaras. e-Stat har fälten
`HTKSYORI` (秘匿処理), `HTKSAKI` (秘匿先) och `GASSAN` (合算) för sekretess och
sammanslagning av de fördelade talen; de rör inte totalen. Läsaren låter ändå en
ruta utgå om totalbefolkningen inte är ett tal, och regeln är densamma 2015 och
2020. Källa: definitionsfilen T000880
(`https://www.e-stat.go.jp/help/data-definition-information/shuroku/T000880.pdf`).

Facit läggs i `data/validation_japan/` (inte den frysta `data/validation/`), på
samma sätt som Kanada. Licens: Government of Japan Standard Terms of Use v2.0,
förenlig med CC BY 4.0; källa anges som "Source: e-Stat".

## Om nedladdningen

`tools/download_ee.py` hämtar varje band i källans eget rutnät (inga pixlar
resamplas) och skriver bredvid varje GeoTIFF en JSON-fil med samling, bild-ID,
år, rutnät, Earth Engines datatyp och SHA-256. Pipelinen hittar årsprodukternas
filer via JSON-filernas år. `--dry-run` visar år, band, rutnät, datatyp och
okomprimerad storlek utan att hämta något. Varje anrop har en tidsgräns på två
minuter och görs om vid fel; en avbruten nedladdning startas om med samma
kommando.
