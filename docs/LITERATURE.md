# Litteratur till artikeln: genomgång och kontroll

Påbörjad 2026-11-17. Syftet är att avgöra vad som redan är känt, så att
artikeln bara påstår det som är nytt. Varje post anger vad källan visar,
var det har kontrollerats och hur den berör våra påståenden. **Status**:
*kontrollerad* = uppgiften lästs i källan (abstract, fulltext eller
förlagets sida); *att kontrollera i fulltext* = uppgiften kommer från
sammanfattning, referat eller andrahandskälla.

Artikelns tre kärnpåståenden (se diskussionen 2026-11-17):

- **P1.** Skarven DMSP–VIIRS är också ett lägesproblem som radiometrisk
  harmonisering inte löser; bakåtskarvning förankrad i VIIRS ger bäst lägen.
- **P2.** Lägesfelet hos ljus som mått på befolkning bestäms av hur ljus per
  invånare samvarierar med läget: x_ljus − x_bef = Cov_P(x, r) / E_P(r).
- **P3.** Ljuset har ett brusgolv: förändring över tio år går inte att mäta på
  nivå 1–2 i de facit vi har.

## 1. Lägesfel i DMSP-OLS

| Källa | Vad den visar | Status | Betydelse för oss |
| --- | --- | --- | --- |
| Tuttle, Anderson, Sutton, Elvidge & Baugh (2013), "It used to be dark here: geolocation calibration of the DMSP OLS". Tidskrift, volym och DOI att fastställa (citeras som *Photogrammetric Engineering & Remote Sensing* 79(3):287–297 i andrahandskällor). | GPS-placerade lampor i Colorado och New Mexico 2009–2011 (F16, F18): systematisk förskjutning i medel 2,9 km, nästan rakt norrut (bäring −0,05°), "about one smooth pixel". Rekommenderar att data flyttas en smooth pixel söderut. Gäller överflygningsdata, inte årskompositerna; 28 punkter. | Kontrollerad (inskickat manus, P. Suttons webbsida). Tidskrift att kontrollera. | Riktningen (norr) stämmer med vår samregistrering (alla 35 kompositer ligger norrut om VIIRS, dy −0,14 … −1,38 px). Lägesfel i DMSP är alltså känt; vårt bidrag är mätningen per komposit mot VIIRS och följden för en serie av lägen. |
| Stathakis (2016), "Intercalibration of DMSP/OLS by parallel regressions", *IEEE Geoscience and Remote Sensing Letters*. | Citeras för att F162009, F182010 och F182011 är förskjutna en halv pixel. | Att kontrollera i fulltext. | Se nästa rad. |
| Geng, Xue, Yan, Xie & Huo (2020), "Shift estimation for DMSP-OLS night-time light rasters in 2009–2011", *J. Remote Sensing & GIS* 9(2):273, doi:10.35248/2469-4134.20.9.273. | F162009, F182010, F182011 saknar TFW-filer och är förskjutna av okänd anledning. FFT-korrelation mot medel av F162008 och F182012: x −0,64 / −1,27 / +0,27 px, y +0,55 / ±0,09 / +0,36 px (positivt = nedåt). | Kontrollerad (förlagets webbsida). Tidskriften är inte en etablerad fjärranalystidskrift. | Vår samregistrering mot VIIRS ger, räknat relativt grannkompositerna, x −0,43 / −1,06 / +0,49 px och y +0,69 / +0,45 / +0,61 px (relativt medlet av F162008 och F182012, som hos Geng m.fl.): samma tecken i x för alla tre och i y för F162009 och F182011; för F182010 anger Geng m.fl. motstridiga tecken i tabell och text. Oberoende stöd för att samregistreringen fungerar. Earth Engine anger hörnet (−180; 75) just för dessa tre kompositer (se HARMONIZATION.md §7.1). |
| Ash & Mazur (2020), "Identifying and correcting signal shift in DMSP-OLS data", *Remote Sensing* 12(14):2219, doi:10.3390/rs12142219. | Systematisk förskjutning av DN-lägen under en satellits livstid (F16 2004–2009), upptäckt i en studie av Syrien; föreslår korrigering med förskjutning av rastret. | Kontrollerad (abstract, MDPI). | Stämmer med vår iakttagelse att F15 2002–2008 och F16 2005–2009 avviker mest och regionalt (avsnitt 7). |
| Kontroll av vår egen inläsning (2026-11-17) | Earth Engines `crs_transform` för DMSP: hörn (−180,00417; 75,00417) för 32 av 35 kompositer, samma som EOG:s rutnät; (−180; 75) för F162009, F182010, F182011. | Kontrollerad i data/raw/dmsp/*.json. | Tolkningen i HARMONIZATION.md §7 punkt 2 (att den genomsnittliga förskjutningen beror på Earth Engines rutnät) är fel och har rättats. |

## 2. Harmoniserade serier och hur de utvärderas

| Källa | Utvärdering | Status | Betydelse |
| --- | --- | --- | --- |
| Li, Zhou, Zhao & Zhao (2020), "A harmonized global nighttime light dataset 1992–2018", *Scientific Data* 7:168, doi:10.1038/s41597-020-0510-y. | Summor och tända pixlar över tid; "consistent temporal trends". Rekommenderar DN > 7. | Kontrollerad (abstract, Nature). | Utvärderas med summor, inte med lägen eller fördelning. Vår kontroll kör både utan och med DN > 7. |
| Zhao, Cao, Chen & Cui (2022), CCNL 1992–2013, *Scientific Data* 9:424, doi:10.1038/s41597-022-01540-x. | Mellanårskonsistens, stadskärnor, blooming. | Kontrollerad (abstract). | Bara DMSP; ingen sensorskarv. |
| Chen m.fl. (2024), "A global annual simulated VIIRS nighttime light dataset from 1992 to 2023", *Scientific Data*, doi:10.1038/s41597-024-04228-6. | Pixel- och regionnivå, R² mot VIIRS. Konstaterar tydliga diskontinuiteter mellan DMSP och VIIRS i rådata. | Kontrollerad (förlagets sida). | Pixelöverensstämmelse, inte lägen. |
| Chiovelli, Michalopoulos, Papaioannou & Regan (2026), "Illuminating the Global South", *The Economic Journal*, doi:10.1093/ej/ueaf134 (advance article jan. 2026). Data: Zenodo doi:10.5281/zenodo.20558766. | Korrigering för mättnad, blooming och sensorbyte; validering mot lokala utvecklingsmått (utbildning, förmögenhet, elektrifiering) i låginkomstländer. Ingen synlig diskontinuitet 2012–2013 enligt författarnas sammanfattning (VoxDev). | Kontrollerad (CEP/LSE-sidan, arbetspapper nov. 2025). Utgivningsår: advance article 2026; zenodo-metadatan säger 2025 – välj 2026 med doi. | Validering mot utvecklingsmått, inte lägen. I vår jämförelse har den fortfarande en skarv 1,5 gånger ett DMSP-årshopp. |
| Ma m.fl.? (2026), "A temporally consistent global 500 m-resolution monthly VIIRS-like nighttime light dataset (1992–2024)", *ESSD* 18:3449, doi:10.5194/essd-18-3449-2026. Författare att fastställa. | R² mot radianskalibrerad DMSP på stadsnivå; nämner att hjälpdata har "spatial misalignments" mot DMSP. | Att kontrollera i fulltext. | Visar att lägesfel uppmärksammas som indata-problem, inte som utvärderingsmått. |
| Gibson, Olivia, Boe-Gibson & Li (2021), "Which night lights data should we use in economics, and where?", *J. Development Economics* 149:102602, doi:10.1016/j.jdeveco.2020.102602. | DMSP är ett dåligt mått på BNP utanför städer; skillnaden mot VIIRS är störst på låg nivå i hierarkin och i glesa områden. | Kontrollerad (abstract). | Stöder att VIIRS bör vara förankringen (P1). |

**Användning av lägesmått över skarven.** Flera studier beräknar ljusets
tyngdpunkt eller standardavvikelseellipser över 1992–2020 med DMSP och VIIRS
ihopfogade, t.ex. Dong m.fl. (2020, *J. Geo-information Science*
22(5):1161–1174, doi:10.12082/dqxxkx.2020.190711; tre kinesiska
stadsregioner, VIIRS omräknad till DMSP-lik med potensfunktion), Xu & Xu
(YRD 1992–2018, samma tidskrift) och Wang m.fl. (2025, *National Remote
Sensing Bulletin*, doi:10.11834/jrs.20254567; tyngdpunkt och
standardavvikelseellips 30 år). Status: kontrollerad (abstract). Det visar att
lägesmått över skarven används i praktiken.

**Sökt men inte funnet (2026-11-17).** Ingen utvärdering av en harmoniserad
DMSP–VIIRS-serie med lägesmått (centroid, spridning) vid skarven. Sökningarna
gällde "harmonized nighttime light spatial consistency 2013", "centroid shift
sensor change", "gravity center nighttime light DMSP VIIRS". Detta är ett
*negativt* sökresultat och ska formuleras försiktigt ("vi har inte funnit").
Bör kompletteras med sökning i Web of Science/Scopus.

## 3. Ljus per invånare, inkomst och befolkningsmodeller (P2)

| Källa | Vad den visar | Status | Betydelse |
| --- | --- | --- | --- |
| Mirza m.fl. (2021), "Global inequality remotely sensed", *PNAS* 118(18), PMC8106331. Författare och doi att fastställa. | Rumslig variation i ljus per invånare återspeglar inkomstfördelningen; Gini-koefficienter av ljus per invånare. | Kontrollerad (abstract, PMC). | Mekanismen (ljus per invånare följer inkomst) är känd. Vårt bidrag är följden för läge. |
| Galimberti, Pichler & Pleninger, "Measuring inequality using geospatial data" (arbetspapper 2020; publicerad version att fastställa). | Ojämlikhetsmått ur ljus per invånare. | Kontrollerad (arbetspapper). | Samma som ovan. |
| Gaughan m.fl.? (2019), "Evaluating nighttime lights and population distribution as proxies for mapping anthropogenic CO2 emission in Vietnam, Cambodia and Laos", *Environ. Res. Commun.* 1:091006, doi:10.1088/2515-7620/ab3d91. Författare att fastställa. | Fördelning med nattljus koncentrerar utsläpp annorlunda än fördelning med befolkning, på subnationell nivå. | Kontrollerad (abstract). | Analog till P2 för en annan tillämpning. |
| Price (1972), "Extension of covariance selection mathematics", *Annals of Human Genetics* 35:485–490. | Identiteten ave_w(x) = x̄ + cov(w, x)/w̄ (ekv. A7). | Kontrollerad (faksimil av s. 486). Doi att fastställa. | **Identiteten i P2 är inte ny**; den är Price-identiteten med w = ljus per invånare. Vårt bidrag är tillämpningen och tolkningen. |
| Niedomysl m.fl. (2017) (se 4) | Lägesfelet i Sverige hade bara signifikant samband med enhetens yta (R² = 0,12). | Enligt HARMONIZATION.md 11.4; att kontrollera i fulltext. | P2 ger en förklaring till varför yta, form och täthet inte förklarar felet. |

## 4. Nattljus som mått på befolkningens läge och förändring (P3)

| Källa | Vad den visar | Status | Betydelse |
| --- | --- | --- | --- |
| Hall, Archila Bustos, Boke Olén & Niedomysl (2019), *Scientific Data* 6:235, doi:10.1038/s41597-019-0250-z. | v1. Validering mot SCB:s rutnät 2007–2008 (Sverige); jämförelse av Aboufadel–Austin och Barmore för 2008 (skillnaden växer med ytan); ingen prövning av projektioner. Anger att data kan användas för att analysera förändringar i befolkningens fördelning och som mått på flyttavstånd. Maskning norr om 65° N berör 7 länder och 0,036 % av befolkningen. | Kontrollerad (fulltext, PMC). | P3 begränsar v1:s påstående om förändring. |
| Niedomysl, Hall, Archila Bustos & Ernstson (2017), "Using satellite data on nighttime lights intensity to estimate contemporary human migration distances", *Annals of the AAG* 107(3):591–605, doi:10.1080/24694452.2016.1270191. | Ljusets centroider ligger nära befolkningens; kan ge flyttavstånd. | Kontrollerad (abstract). | Lägen i tvärsnitt (flyttavstånd) påverkas inte av P3; förändring gör det. |
| Gibson m.fl. (Waikato arbetspapper 2407) | Ljus förutsäger lokala förändringar i ekonomisk aktivitet dåligt, även där tvärsnittet förutsägs väl. | Att kontrollera (arbetspapper; titel att fastställa). | Parallell till P3 för ekonomisk aktivitet. |

## 5. Att göra

- Fastställa tidskrift, volym och doi för Tuttle m.fl. 2013, Stathakis 2016,
  Mirza m.fl. 2021, Price 1972, ESSD 18:3449 och Gaughan m.fl. 2019.
- Läsa Stathakis 2016 och Ash & Mazur 2020 i fulltext: storlek och riktning
  för förskjutningarna, jämfört med config/dmsp_offsets.csv.
- Systematisk sökning (Web of Science/Scopus) efter utvärdering av
  harmoniserade serier med lägesmått.
- Litteratur om befolkningsmodeller som använder nattljus som förklaring
  (WorldPop, LandScan, GHS-POP): finns studier av skevhet mot rikare områden?
