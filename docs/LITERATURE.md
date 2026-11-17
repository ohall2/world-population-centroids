# Litteratur till artikeln: genomgång och kontroll

Påbörjad 2026-11-17. Syftet är att avgöra vad som redan är känt, så att
artikeln bara påstår det som är nytt. Varje post anger vad källan visar,
var det har kontrollerats och hur den berör våra påståenden. **Status**:
*kontrollerad* = uppgiften lästs i källan (abstract, fulltext eller
förlagets sida); *att kontrollera i fulltext* = uppgiften kommer från
sammanfattning, referat eller andrahandskälla.

Artikelns tre kärnpåståenden (se diskussionen 2026-11-17):

- **P1.** *(omformulerad efter HARMONIZATION.md §12)* De harmoniserade
  serierna hänger ihop globalt men inte enhet för enhet: vid skarven hoppar
  alla mått per enhet (läge, spridning, koncentration, tänd yta, total
  ljusmängd) 2–3 gånger mer än ett vanligt år. Aggregerad eller pixelvis
  utvärdering avslöjar det inte. För lägen ger bakåtskarvning förankrad i
  VIIRS bäst resultat mot facit.
- **P2.** Lägesfelet hos ljus som mått på befolkning bestäms av hur ljus per
  invånare samvarierar med läget: x_ljus − x_bef = Cov_P(x, r) / E_P(r).
- **P3.** Ljuset har ett brusgolv: förändring över tio år går inte att mäta på
  nivå 1–2 i de facit vi har.

## 1. Lägesfel i DMSP-OLS

| Källa | Vad den visar | Status | Betydelse för oss |
| --- | --- | --- | --- |
| Tuttle, Anderson, Sutton, Elvidge & Baugh (2013), "It used to be dark here: geolocation calibration of the DMSP OLS". Tidskrift, volym och DOI att fastställa (citeras som *Photogrammetric Engineering & Remote Sensing* 79(3):287–297 i andrahandskällor). | GPS-placerade lampor i Colorado och New Mexico 2009–2011 (F16, F18): systematisk förskjutning i medel 2,9 km, nästan rakt norrut (bäring −0,05°), "about one smooth pixel". Rekommenderar att data flyttas en smooth pixel söderut. Gäller överflygningsdata, inte årskompositerna; 28 punkter. | Kontrollerad (inskickat manus, P. Suttons webbsida). Tidskrift att kontrollera. | Riktningen (norr) stämmer med vår samregistrering (alla 35 kompositer ligger norrut om VIIRS, dy −0,14 … −1,38 px). Lägesfel i DMSP är alltså känt; vårt bidrag är mätningen per komposit mot VIIRS och följden för en serie av lägen. |
| Stathakis (2016), "Intercalibration of DMSP/OLS by parallel regressions", *IEEE Geoscience and Remote Sensing Letters* 13:1421–1424 (enligt referenslistan i Li & Zhou 2017, *Remote Sensing* 9:637). | Citeras för att F162009, F182010 och F182011 är förskjutna en halv pixel. | Att kontrollera i fulltext. | Se nästa rad. |
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
| Cheng, Geng, Li, Li, Zhao, Lin, Wang, Gong & Zhou (2026), "A temporally consistent global 500 m-resolution monthly VIIRS-like nighttime light dataset (1992–2024)", *ESSD* 18:3449, doi:10.5194/essd-18-3449-2026. | R² mot radianskalibrerad DMSP på stadsnivå; nämner att hjälpdata har "spatial misalignments" mot DMSP. | Att kontrollera i fulltext. | Visar att lägesfel uppmärksammas som indata-problem, inte som utvärderingsmått. |
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
| Mirza, Xu, van Bavel, van Nes & Scheffer (2021), "Global inequality remotely sensed", *PNAS* 118(18):e1919913118, doi:10.1073/pnas.1919913118. | Rumslig variation i ljus per invånare återspeglar inkomstfördelningen; motiveras med att ljus per person ökar med inkomsten (större, bättre upplysta hus och områden). Gini-koefficienter av ljus per invånare. | Kontrollerad (PMC). | Mekanismen (ljus per invånare följer inkomst) är känd. Vårt bidrag är följden för läge. |
| Galimberti, Pichler & Pleninger, "Measuring inequality using geospatial data" (arbetspapper 2020; publicerad version att fastställa). | Ojämlikhetsmått ur ljus per invånare. | Kontrollerad (arbetspapper). | Samma som ovan. |
| Gaughan m.fl.? (2019), "Evaluating nighttime lights and population distribution as proxies for mapping anthropogenic CO2 emission in Vietnam, Cambodia and Laos", *Environ. Res. Commun.* 1:091006, doi:10.1088/2515-7620/ab3d91. Författare att fastställa. | Fördelning med nattljus koncentrerar utsläpp annorlunda än fördelning med befolkning, på subnationell nivå. | Kontrollerad (abstract). | Analog till P2 för en annan tillämpning. |
| Price (1972), "Extension of covariance selection mathematics", *Annals of Human Genetics* 35(4):485–490. Jfr Price (1970), "Selection and covariance", *Nature* 227:520–521, doi:10.1038/227520a0. | Identiteten ave_w(x) = x̄ + cov(w, x)/w̄ (ekv. A7). | Kontrollerad (faksimil av s. 485–486). Doi för 1972 att fastställa (Scilit anger en uppenbart felaktig doi). | **Identiteten i P2 är inte ny**; den är Price-identiteten med w = ljus per invånare. Vårt bidrag är tillämpningen och tolkningen. |
| Niedomysl m.fl. (2017) (se 4) | Lägesfelet i Sverige hade bara signifikant samband med enhetens yta (R² = 0,12). | Enligt HARMONIZATION.md 11.4; att kontrollera i fulltext. | P2 ger en förklaring till varför yta, form och täthet inte förklarar felet. |

## 4. Nattljus som mått på befolkningens läge och förändring (P3)

| Källa | Vad den visar | Status | Betydelse |
| --- | --- | --- | --- |
| Hall, Archila Bustos, Boke Olén & Niedomysl (2019), *Scientific Data* 6:235, doi:10.1038/s41597-019-0250-z. | v1. Validering mot SCB:s rutnät 2007–2008 (Sverige); jämförelse av Aboufadel–Austin och Barmore för 2008 (skillnaden växer med ytan); ingen prövning av projektioner. Anger att data kan användas för att analysera förändringar i befolkningens fördelning och som mått på flyttavstånd. Maskning norr om 65° N berör 7 länder och 0,036 % av befolkningen. | Kontrollerad (fulltext, PMC). | P3 begränsar v1:s påstående om förändring. |
| Niedomysl, Hall, Archila Bustos & Ernstson (2017), "Using satellite data on nighttime lights intensity to estimate contemporary human migration distances", *Annals of the AAG* 107(3):591–605, doi:10.1080/24694452.2016.1270191. | Ljusets centroider ligger nära befolkningens; kan ge flyttavstånd. | Kontrollerad (abstract). | Lägen i tvärsnitt (flyttavstånd) påverkas inte av P3; förändring gör det. |
| Gibson m.fl. (Waikato arbetspapper 2407) | Ljus förutsäger lokala förändringar i ekonomisk aktivitet dåligt, även där tvärsnittet förutsägs väl. | Att kontrollera (arbetspapper; titel att fastställa). | Parallell till P3 för ekonomisk aktivitet. |

## 5. Konsistens per enhet i harmoniserade serier (P1, tillägg 2026-11-17)

| Källa | Vad den visar | Status | Betydelse |
| --- | --- | --- | --- |
| Li, S., Cao, Zhao, Jie, Liu, Chen & Cui (2023), "Developing a pixel-scale corrected nighttime light dataset (PCNL, 1992–2021) combining DMSP-OLS and NPP-VIIRS", *Remote Sensing* 15(16):3925, doi:10.3390/rs15163925. Data: Zenodo record 7612389. | Jämför PCNL med Li m.fl. 2020 och Chen m.fl. 2021: regional summa av ljus (TSOL) stiger abrupt 2013–2014 i Li-serien och 2010–2013 i Chen-serien. PCNL utvärderas för konsistens på global, nationell, stads- och pixelnivå. | Kontrollerad (abstract och sammanfattning i sökresultat; fulltext kunde inte hämtas). | **Delvis föregångare till P1**: ett hopp i regionala summor i Li-serien vid skarven är redan rapporterat. Vårt bidrag: alla enheter på nivå 1 i världen, fem mått, två oberoende serier, och jämförelse med ett vanligt årshopp. PCNL är tillagd i HARMONIZATION.md §12: skarven per enhet är nästan borta (kvot 0,85–1,6 utom tänd yta 4,1). |
| Li, X., Lei, Shen, Xin & Qu (2026), "Spatiotemporal dynamics and relative growth of county-level CO2 emissions in China (1997–2021) derived from harmonized nighttime lights", *Environmental Monitoring and Assessment* 198(12):1287, doi:10.1007/s10661-026-16129-9. | Konstaterar att den harmoniserade serien saknar hopp 2013–2014 i aggregat, men att det "does not, by itself, establish temporal stability at the county level". | Kontrollerad (förlagets sida). | Visar att problemet uppmärksammas som förbehåll men inte mäts. |
| Zhao, N., Zhou & Samson (2015), "Correcting incompatible DN values and geometric errors in nighttime lights time-series images", *IEEE TGRS* 53(4):2039–2049, doi:10.1109/TGRS.2014.2352598. | Geometrisk korrektion som ett av fyra steg för DMSP-tidsserier i Kina och USA 1992–2008. | Kontrollerad (abstract); metod och storlek på felen att läsa i fulltext. | Föregångare för geometrisk korrektion av DMSP mellan år; inte mot VIIRS och inte globalt. |
| Shen, Zhu, Cao & Chen (2019), "Measurement of blooming effect of DMSP-OLS nighttime light data based on NPP-VIIRS data", *Annals of GIS* 25(2):153–165, doi:10.1080/19475683.2019.1570336. | Registrerar DMSP mot VIIRS för Kina innan blooming mäts; blooming underskattar rumslig heterogenitet med 90 % där VIIRS > 10. | Kontrollerad (abstract och sidhuvud; storleken på registreringsförskjutningen att läsa i fulltext). | Visar att DMSP och VIIRS inte ligger kongruent; registrering görs men redovisas inte som felkälla för tidsserier. |
| Chen, Z. m.fl. (2021), "An extended time series (2000–2018) of global NPP-VIIRS-like nighttime light data from a cross-sensor calibration", *ESSD* 13:889. | Inga plötsliga hopp i total ljusmängd och tända pixlar 2010–2014 (globalt). | Kontrollerad (förlagets sida). | Ännu ett exempel på global utvärdering. |

## 6. Ljus som förklaring i befolkningsmodeller (P2, tillägg)

| Källa | Vad den visar | Status | Betydelse |
| --- | --- | --- | --- |
| Nieves, Bondarenko, Kerr m.fl. & Tatem (2021), "Measuring the contribution of built-settlement data to global population mapping", *Social Sciences & Humanities Open* 3:100102, doi:10.1016/j.ssaho.2020.100102. | I WorldPops slumpskogsmodeller för 172 länder 2000–2012 är nattljus den viktigaste förklaringsvariabeln för befolkningstäthet, med bebyggelse på andra plats. | Kontrollerad (abstract, PMC). | Om ljus per invånare följer inkomst (Mirza m.fl.) och nattljus är den viktigaste variabeln, kan befolkningsmodeller ärva samma skevhet som P2 beskriver. **Hypotes, inte prövad av oss.** |
| Thomson, Leasure, Bird, Tzavidis & Tatem (2021, preprint), "How accurate are WorldPop-Global gridded population data at the cell-level?" | Stora fel på cellnivå i en simulering för Namibia, särskilt i städer. | Kontrollerad (preprint). | Bakgrund. |

## 7. Tvärsnitt mot förändring (P3, tillägg)

Flera studier i ekonomi visar samma mönster som P3: nattljus förklarar
nivåer i tvärsnitt väl men förändring dåligt.

| Källa | Vad den visar | Status |
| --- | --- | --- |
| Asher, Lunt, Matsuura & Novosad (2021), "Development research at high geographic resolution", *World Bank Economic Review* 35(4):845–871, doi:10.1093/wber/lhab003. | På by-nivå i Indien är nattljus signifikanta mått på befolkning, sysselsättning och konsumtion, men elasticiteterna är mycket lägre i tidsserie än i tvärsnitt. | Kontrollerad (abstract). |
| Goldblatt, Heilmann & Vaizman (2020), "Can medium-resolution satellite imagery measure economic activity at small geographies? Evidence from Landsat in Vietnam", *World Bank Economic Review* 34(3):635–653, doi:10.1093/wber/lhz001. | Varken nattljus eller Landsat förutsäger förändringar över tid väl. | Kontrollerad (abstract). |
| Martinez, MacManus, Stokes, Wang & de Sherbinin (2023), "Suitability of NASA's Black Marble daily nighttime lights for population studies at varying spatial and temporal scales", *Remote Sensing* 15(10):2611. | Minskande landsbygdsbefolkning gav ingen minskning i nattljus trots betydande utflyttning. | Kontrollerad (abstract). |
| Bluhm & McCord (2022), "What can we learn from nighttime lights for small geographies? Measurement errors and heterogeneous elasticities", *Remote Sensing* 14(5):1190, doi:10.3390/rs14051190. | Ljusets förmåga att fånga förändringar beror på lokala förhållanden. | Kontrollerad (abstract). |
| Bennett & Smith (2017), "Advances in using multitemporal night-time lights satellite imagery to detect, estimate, and monitor socioeconomic dynamics", *Remote Sensing of Environment* 192:176–197, doi:10.1016/j.rse.2017.01.005. | Översikt; hävdar att nattljus i tidsserier kan vara ett tillförlitligt mått på förändring. | Kontrollerad (abstract). **Motsatt ståndpunkt**, som P3 måste förhålla sig till. |

**Följd för P3.** Påståendet är inte nytt som allmän iakttagelse. Det nya
är att det prövas för befolkningens *läge* mot folkräkningar i tre system,
med ett mått på brusgolvet (fel i förändring jämfört med verklig förflyttning).

## 8. Tyngdpunkter som mått (bakgrund)

| Källa | Innehåll | Status |
| --- | --- | --- |
| Quah (2011), "The global economy's shifting centre of gravity", *Global Policy* 2(1):3–9. | Världsekonomins tyngdpunkt ur BNP på omkring 700 platser; flyttade från Atlanten 1980 till öster om Helsingfors och Bukarest 2008. | Kontrollerad (LSE). |
| Grether & Mathys (2011), "On the track of the world's economic center of gravity", i *Economic Growth and Development* (Frontiers of Economics and Globalization 11), s. 261–287, doi:10.1108/S1574-8715(2011)0000011015. | Ekonomisk tyngdpunkt 1950–2008 ur G-Econ; produktionen är mindre koncentrerad än befolkningen. | Kontrollerad (förlaget). |

Visar att tyngdpunkter är ett etablerat sätt att sammanfatta fördelningar,
och att de beräknas ur rutnätsdata som själva kan vara modellerade.

## 9. Vad som är nytt, påstående för påstående (läge 2026-11-17)

- **P1.** *Känt:* lägesfel i DMSP (Tuttle m.fl. 2013; Ash & Mazur 2020;
  Geng m.fl. 2020), att DMSP och VIIRS måste registreras mot varandra
  (Shen m.fl. 2019), hopp i regionala summor i Li-serien (Li, S. m.fl. 2023)
  och förbehåll om stabilitet per enhet (Li, X. m.fl. 2026). *Nytt:* en
  global mätning per enhet med fem mått och två serier, som visar att
  skarven gäller alla mått, även summor; samregistrering per komposit mot
  VIIRS; att bakåtskarvning förankrad i VIIRS ger bäst lägen mot tre facit.
- **P2.** *Känt:* identiteten (Price 1970, 1972); att ljus per invånare
  följer inkomst (Mirza m.fl. 2021). *Nytt:* att identiteten förklarar
  lägesfelet hos ljus som mått på befolkning, att felet därför är
  förutsägbart och att det går åt motsatta håll i Mexiko och Norden under
  samma omvandling. *Prövat 2026-11-17 (HARMONIZATION.md §14):* skillnader i
  ljus per invånare mellan regioner står för 110–146 av 129–149 km i Mexiko,
  Sverige och Norge (rutnät 10′).
- **P3.** *Känt:* nattljus fångar förändring sämre än nivåer (Asher m.fl.
  2021; Goldblatt m.fl. 2020; Martinez m.fl. 2023). *Nytt:* att det gäller
  befolkningens läge, mätt mot folkräkningar, och med ett brusgolv i km.

## 10. Att göra

- Fastställa tidskrift och doi för Tuttle m.fl. 2013 (Suttons CV anger
  "in press"; ett kinesiskt index anger s. 287–297), doi för Price 1972,
  författare för CO2-artikeln i *Environ. Res. Commun.* 2019, doi för
  Stathakis 2016.
- Läsa Li, S. m.fl. 2023 i fulltext: hur stort är hoppet per region, och
  vilka enheter?
- Läsa Stathakis 2016 och Ash & Mazur 2020 i fulltext: storlek och riktning
  för förskjutningarna, jämfört med config/dmsp_offsets.csv.
- Systematisk sökning (Web of Science/Scopus) efter utvärdering av
  harmoniserade serier med lägesmått.
- Litteratur om befolkningsmodeller som använder nattljus som förklaring
  (WorldPop, LandScan, GHS-POP): finns studier av skevhet mot rikare områden?
