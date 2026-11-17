# Metod för v2 och skillnader mot v1

Det här dokumentet beskriver hur v2 beräknas och varför varje val gjordes. Det är också underlag för metodavsnittet i en kommande artikel.

## 1. Ljusdata: DMSP-OLS v4 Stable Lights

- **Källa.** De 35 satellit-års-kompositerna hämtas från Google Earth Engine (`NOAA/DMSP-OLS/NIGHTTIME_LIGHTS`, bandet `stable_lights`) med [`tools/download_ee.py`](../tools/download_ee.py). Exporten görs i källans eget rutnät utan omsampling. Varje fil får en metadatafil med checksumma.
- **Rutnät.** DMSP v4 har 43 201 × 16 801 pixlar à 30″, med pixelmittpunkter på hela 30″-steg från 180°V till 180°O och från 75°N till 65°S. Pixelhörnet ligger därför en halv pixel utanför, i (−180,0041667°, 75,0041667°). Earth Engine anger i stället hörnet till (−180°, 75°) och skalan avrundad till 0,0083333333°, vilket flyttar hela rutnätet en halv pixel (cirka 460 m) mot sydost. v2 använder EOG:s rutnät ([`wpc/lights.py`](../wpc/lights.py), `DMSP_TRANSFORM`).
- **F152008.** Den kompositen finns i Earth Engine men saknas i v1:s kalibreringstabell och används inte.

## 2. Kalibrering

- Koefficienterna är desamma som i v1 (Hall m.fl. 2019, tabell 2): IRQR enligt Elvidge m.fl., med F18 2010 som referens och Sicilien som invariant region. Se [`config/calibration_v1.csv`](../config/calibration_v1.csv).
- Kalibrerat värde = c0 + c1·DN + c2·DN². **Pixlar med DN = 0 förblir 0.** Annars skulle varje mörk pixel få värdet c0 (0,01–0,23) och få vikt i centroiden. Det framgår inte av v1-artikeln hur detta hanterades.
- År med två satelliter: medelvärdet av de kalibrerade värdena.
- Observera att andragradsfunktionen inte är monoton för alla satelliter. För F14 1997 når den sitt maximum vid DN ≈ 54. Detta ärvs från v1:s koefficienter.

## 3. Gasfacklor

- Polygonerna från Elvidge m.fl. (2009), version 2009-06-18, en fil per land. De hämtas och slås ihop med [`tools/prepare_flares.py`](../tools/prepare_flares.py). Det är samma källa som v1 använde.
- Tre av 65 filer (Elfenbenskusten, Ghana, Mauretanien) har koordinater som inte är grader och utesluts. Alla tre gäller facklor till havs, utanför landpolygonerna, och påverkar därför inga centroider.
- Pixlar som berörs av en fackelpolygon sätts till 0.
- Polygonerna täcker facklor som syntes fram till 2008. Nya facklor efter det fångas inte.
- **Huvudserien** använder i stället masken `flares_viirs_r2km.gpkg` för alla år och båda sensorerna: 2009-polygonerna plus en cirkel med radien 2 km kring varje fackelplats i VIIRS Nightfire (avsnittet *Fackelmask för VIIRS* nedan). Masken är statisk, så samma pixlar tas bort alla år. 2009-polygonerna ensamma används bara i jämförelsen med v1 och för CCNL.

## 4. Gränser

- **v2:** GADM 4.1, lagren ADM_0–ADM_2 i `gadm_410-levels.gpkg` (263, 3 662 respektive 47 217 enheter). Samma gränser används för alla år.
- **Jämförelse med v1:** GADM 2.8. GADM:s server för 2.8 finns inte längre. Geometrierna kommer från NOAA GML:s spegel, adm0-tabellen från Zenodo (doi:10.5281/zenodo.5103739), och attributen för adm1–2 återskapas ur v1:s utdata. Radordningen är kontrollerad: varje v1-centroid som v1 placerade i sin egen enhet ligger i polygonen med samma OBJECTID. Se [`tools/prepare_gadm28.py`](../tools/prepare_gadm28.py).

## 5. Vikter

För varje enhet och pixel:

> vikt = ljus × andel av pixeln som ligger i enheten × pixelns yta på sfären

- **Andel av pixeln** beräknas med exactextract. v1 räknade en pixel helt eller inte alls, beroende på om pixelns mittpunkt låg i polygonen. Med andelar får små enheter alltid sina pixlar med, och en pixel som delas mellan två enheter delas proportionellt.
- **Pixelns yta** = R²·Δλ·(sin φ_nord − sin φ_syd), med R = 6 371 km. Ljusvärdet är en intensitet per yta, och en pixel på 30″ krymper mot polerna (vid 69°N är den cirka 0,6 gånger så stor som vid 55°N). v1 viktade inte efter yta, vilket ger nordliga pixlar för stor vikt. För Sverige 2010 flyttar ytviktningen centroiden 30 km söderut.

## 6. Centroider

- **Aboufadel–Austin:** pixelmittpunkterna omvandlas till 3D-enhetsvektorer, och det viktade medelvärdet projiceras tillbaka till sfären. Metoden är sluten och hanterar datumlinjen korrekt.
- **Barmore:** viktat medel av x och y i en azimutal ekvidistant projektion (sfär, R = 6 371 km), centrerad i föregående skattning. Steget upprepas tills förflyttningen är under 1 m.
  - v2 startar i Aboufadel–Austin-centroiden. v1 startade i medelvärdet av longitud och latitud, vilket blir fel för enheter som korsar datumlinjen.
  - v2 stoppar när *avståndet* är under 1 m. v1 stoppade när *någon* av komponenterna var under 1 m (`&&` i while-villkoret).
- Metoderna skiljer sig mest för stora enheter: Ryssland cirka 25 km, USA cirka 8 km och Kanada cirka 2 km 2010, medan medianen är några meter. Samma mönster finns i v1.
- **Stora enheter** räknas ruta för ruta (2 048 × 2 048 pixlar). Rutorna följer pixelgränserna, så resultatet blir identiskt med en beräkning i ett svep. Det är testat.

## 7. Flaggor

Samma kodning som `OUT_FLAG` i v1: 0 = centroiden ligger i den egna enheten, 1 = utanför alla enheter (till exempel i havet), 2 = i en annan enhet på samma nivå. v2 lägger till −1 = ingen centroid, eftersom enheten saknar ljus.

## 8. Förhållandet till v1

v1 går inte att återskapa exakt med Earth Engine-rastren. Ett test för 2010, där kalibreringen är identitet, gav avvikelser på 50–1 000 m även i länder utan facklor, och avvikelserna gick åt olika håll i olika länder. Det tyder på att pixelvärdena skiljer sig från dem v1 använde, till exempel genom en annan bearbetningsversion av DMSP v4. Koden följer v1:s steg rad för rad ([`wpc/v1method.py`](../wpc/v1method.py)). Skillnaden mellan v1 och v2 redovisas därför som en jämförelse, inte som en reproduktion.

## 9. Huvudserien

Den publicerade serien (`tools/build_release.py`) byggs i fyra steg. Resonemang, prövade alternativ och validering står i [`docs/HARMONIZATION.md`](HARMONIZATION.md).

1. **Samregistrering av DMSP** (avsnitt 7 där): varje komposit flyttas med subpixel-noggrannhet till VIIRS-läget 2013 ([`config/dmsp_offsets.csv`](../config/dmsp_offsets.csv)). Källan heter `dmsp-irqr-reg`.
2. **Rå VIIRS VNL v2 2013–2025**, utan omvandling till DMSP-skala (källan `viirs`).
3. **Skarvning bakåt per enhet** (avsnitt 8–10 där): DMSP-åren 1992–2012 flyttas med enhetens skillnad mellan rå VIIRS och DMSP 2013, i km österut och norrut. VIIRS-åren lämnas orörda.
4. **Maskning norr om 65° N på nivå 0**, som i v1 (avsnitt 11 där): nationella centroider beräknas utan ljus norr om 65° N (`wpc.run --max-lat 65`). Nivå 1 och 2 använder allt ljus.

Valideringen mot US Census, Eurostat GEOSTAT och INEGI ger medianfel på 0,7–2,9 km på nivå 1 utom i Mexiko (7–9 km), 0,7–1,5 km på nivå 2 och 11–44 km på nivå 0 utom i Mexiko (omkring 130 km). Tabellerna finns i docs/validation/ och i datamängdens metadata.json.

## Referenser

- Hall, O., Archila Bustos, M. F., Boke Olén, N. & Niedomysl, T. (2019). Population centroids of the world administrative units from nighttime lights 1992–2013. *Scientific Data* 6, 235. doi:10.1038/s41597-019-0250-z.
- Elvidge, C. D. m.fl. (2009). A fifteen year record of global natural gas flaring derived from satellite data. *Energies* 2, 595–622.
- Elvidge, C. D. m.fl. (2014). National trends in satellite-observed lighting 1992–2012. I *Global Urban Monitoring and Assessment through Earth Observation*. CRC Press.
- Chiovelli, G., Michalopoulos, S., Papaioannou, E. & Regan, T. (2026). Illuminating the Global South. *The Economic Journal*. doi:10.1093/ej/ueaf134.
- Li, X., Zhou, Y., Zhao, M. & Zhao, X. (2020). A harmonized global nighttime light dataset 1992–2018. *Scientific Data* 7, 168. doi:10.1038/s41597-020-0510-y.
- Niedomysl, T., Hall, O., Archila Bustos, M. F. & Ernstson, U. (2017). Using satellite data on nighttime lights intensity to estimate contemporary human migration distances. *Annals of the American Association of Geographers* 107(3), 591–605. doi:10.1080/24694452.2016.1270191.
- Pandey, B., Zhang, Q. & Seto, K. C. (2017). Comparative evaluation of relative calibration methods for DMSP/OLS nighttime lights. *Remote Sensing of Environment* 195, 67–78.
- Zhao, C., Cao, X., Chen, X. & Cui, X. (2022). A consistent and corrected nighttime light dataset (CCNL 1992–2013) from DMSP-OLS data. *Scientific Data* 9, 424. doi:10.1038/s41597-022-01540-x.


## CCNL (Zhao m.fl. 2022) som alternativ ljuskälla

- Källa: originalfilerna på Zenodo (doi:10.5281/zenodo.6644980), 30″, en fil
  per år 1992–2013. Versionen i Earth Engine (`BNU/FGS/CCNL/v1`) är omsamplad
  till 0,008983° (1 km vid ekvatorn) och används inte.
- Filerna har 43 201 × 16 802 pixlar (en rad mer än DMSP v4, längst i söder),
  float64, nodata −1,797·10³⁰⁸ (förekommer inte i 2010-filen). Skalan anges
  avrundad till 0,0083333333; koden använder exakt 1/120° och hörnet
  (−180,0041667; 75,0041667), dvs. samma rutnät som DMSP.
- Värdena är DN efter korrigering för mellanårsskillnader, mättnad och
  blooming; de är inte heltal och kan överstiga 63 (max 238 år 2010). De
  används direkt som ljusvikt, utan ytterligare kalibrering.
- Korrigeringen lämnar många mycket små positiva värden (10⁻¹⁵–10⁻⁶). De
  påverkar inte centroiderna märkbart men räknas som lysande pixlar i
  kolumnen `n_lit`.
- Gasfacklor maskas med samma polygoner som för DMSP (Elvidge m.fl. 2009).

### Jämförelse CCNL – DMSP-IRQR (GADM 4.1, Barmore, 1992–2013)

| Nivå | Medianhopp per år, DMSP | Medianhopp, CCNL | 99 %-percentil DMSP / CCNL | Avstånd DMSP–CCNL samma år, median (90 %) |
| --- | --- | --- | --- | --- |
| adm0 | 4,2 km | 4,3 km | 164 / 157 km | 4,4 km (24,6 km) |
| adm1 | 1,3 km | 1,3 km | 46 / 47 km | 0,9 km (5,5 km) |
| adm2 | 0,53 km | 0,58 km | 14,3 / 14,8 km | 0,3 km (1,3 km) |

CCNL minskar inte hoppen mellan åren. De stora hoppen kommer samma år i båda
källorna (t.ex. Norge 1993, 1997, 2000; Sverige 2003; Kanada 1995–1998), och
bara en tredjedel av enheterna får lägre medianhopp med CCNL. Bruset ligger
alltså i de underliggande DMSP-kompositerna och inte i IRQR-kalibreringen.
Avståndet mellan källorna samma år är ett mått på hur känsliga centroiderna
är för valet av korrigering; det är störst där mättnad och blooming är stora
(Brasilien 55 km, Kanada 48 km, Ryssland 27 km på nationell nivå).

## VIIRS VNL v2 (2013–2025)

- Källa: Earth Engine `NOAA/VIIRS/DNB/ANNUAL_V21` (2013–2021) och
  `ANNUAL_V22` (2022–2025), bandet `average_masked` (nW/cm²/sr), i källans
  eget rutnät: 86 401 × 33 601 pixlar à 15″, hörn (−180,0020833; 75,0020833).
  Skalan anges avrundad till 0,0041666667; koden använder exakt 1/240°.
  2012 saknas i Earth Engine (och bygger bara på april–december).
- Negativa värden (47 905 pixlar 2013) sätts till 0.
- Ljuset vägs med cellarea som för DMSP; ingen kalibrering mellan åren.

### Fackelmask för VIIRS

Gasfacklor dominerar VIIRS-ljuset i oljeländer (2013: Nigerdeltat tre gånger
Stockholmsområdets ljussumma; högsta pixel 104 816 nW/cm²/sr). Polygonerna
från 2009 täcker 95 % av fackelljuset i Nigerdeltat men 0 % i Bakken och 10 %
i Permian. Masken byggs därför av fackelplatser från VIIRS Nightfire
(MYC25Q4 2012–2025 och EOG:s årskataloger 2012–2024, inklusive LNG,
gasbehandling och raffinaderier) plus 2009-polygonerna; varje plats får en
sfärisk cirkel med radien r (tools/prepare_flares_viirs.py). Masken är
statisk, som för DMSP.

Val av radie (VIIRS 2013, andel av ljuset i rutan som maskas):

| Område | 1 km | 2 km | 3 km | 5 km |
| --- | --- | --- | --- | --- |
| Bakken (fält) | 91 % | 99 % | 99 % | 100 % |
| Permian (fält) | 72 % | 83 % | 89 % | 96 % |
| Edmonton | 3 % | 7 % | 12 % | 24 % |
| Houston | 7 % | 10 % | 15 % | 26 % |
| Oklahoma City | 2 % | 6 % | 13 % | 31 % |
| Midland-Odessa (i fältet) | 37 % | 65 % | 84 % | 100 % |

Över 2 km växer andelen maskat stadsljus snabbare än andelen fackelljus.
Nationella centroider (adm0, Barmore) flyttar sig mest mellan 2009-polygonerna
och 1 km (2024: Ryssland 310 km, Kanada 91 km, Mexiko 96 km); mellan 1 och
2 km högst 58 km (2013) och 31 km (2024); därefter jämnt med radien i länder
där facklor ligger nära städer (Kanada, USA, Kina). **Standard är 2 km**;
1 och 3 km redovisas som känslighetsanalys (data/sensitivity/).

## Harmonisering VIIRS → DMSP (källan `viirs-dmsp`)

DMSP och VIIRS ger olika centroider samma år. För 2013 är medianavståndet
mellan dem 12,5 km för land, 3,6 km för adm1 och 1,3 km för adm2; för stora
länder upp till 150–250 km (Indien 146, Brasilien 159, Kanada 253 km).
Orsaken är sensorerna: DMSP är mättat i stadskärnor, har en
detektionströskel och smetar ut ljus över flera kilometer; VIIRS har högre
upplösning, ingen mättnad och ser svagare ljus. Att foga ihop serierna
direkt skulle ge ett falskt hopp mellan 2013 och 2014.

VIIRS förs därför över till DMSP-IRQR:s skala (tools/harmonize_viirs.py,
wpc/harmonize.py), i linje med Li m.fl. (2020, Scientific Data 7:168) men med
en icke-parametrisk överföringsfunktion:

1. VIIRS (15″) aggregeras ytviktat till DMSP:s 30″-rutnät. Rutnäten är
   förskjutna en halv VIIRS-cell; varje DMSP-cell består av en hel VIIRS-cell
   och halva grannceller (faltning med [¼, ½, ¼] i båda led).
2. Bilden smetas ut med en gaussisk kärna (σ i DMSP-pixlar).
3. För 2013 och för klasser av log10(utsmetad radians) beräknas andelen
   pixlar som är tända i DMSP-IRQR och medelvärdet bland de tända. Båda görs
   monotona (viktad pool-adjacent-violators). En pixel får medelvärdet om
   andelen är minst ½, annars 0, vilket efterliknar DMSP:s tröskel.
4. σ väljs med tvåfaldig rumslig korsvalidering (10°-band i longitud,
   varannat band i varje vikning). Fackelmasken (2 km + 2009-polygoner)
   tillämpas på VIIRS före aggregeringen och på DMSP i anpassningen.
5. Samma σ och funktion tillämpas på alla VIIRS-år 2013–2025.

Provkörning på ett band 55,7–60° N: korsvaliderat R² 0,80 (σ = 0,5), 0,84
(1), 0,87 (1,5), 0,89 (2), 0,892 (2,5) och 0,88 (3); vald σ = 2,5 (ca 2,3 km).
I bandet blev ljussumman 97 % av DMSP:s, korrelationen per pixel 0,94 och den
ljusviktade medelbreddgraden över Sverige och västra Ryssland lika inom
0,02°.

### Resultat för hela jorden (körning 2026-11-17)

Korsvaliderat R² per σ (DMSP-pixlar): 0 → 0,790; 0,5 → 0,820; 1 → 0,868;
1,5 → 0,895; **2 → 0,900**; 2,5 → 0,891; 3 → 0,875. Vald σ = 2,0 (ca 1,9 km).

Avstånd mellan centroider (Barmore, km; median / 90 %):

| | adm0 | adm1 | adm2 |
| --- | --- | --- | --- |
| DMSP 2013 – VIIRS rå 2013 | 12,5 / 65,5 | 3,6 / 19,9 | 1,34 / 5,2 |
| DMSP 2013 – VIIRS harmoniserad 2013 | 7,7 / 44,2 | 2,2 / 12,6 | 0,80 / 3,7 |
| Hopp DMSP 2013 → harmoniserad 2014 | 8,2 / 44,2 | 2,2 / 13,1 | 0,85 / 3,9 |
| DMSP, vanliga årshopp 2004–2013 | 3,9 / 23,8 | 1,2 / 8,0 | 0,53 / 2,5 |
| Harmoniserad, årshopp 2013–2025 | 2,3 / 14,8 | 0,7 / 5,1 | 0,25 / 1,7 |

Harmoniseringen minskar brottet mellan sensorerna med omkring 40 % på alla
nivåer, men hoppet 2013 → 2014 är fortfarande ungefär dubbelt så stort som
ett vanligt DMSP-årshopp. Stora länder: Brasilien 159 → 48 km, Ryssland
148 → 62 km, Argentina 85 → 19 km, Indien 146 → 127 km, Kina 106 → 77 km;
Kanada blir sämre (253 → 325 km) och Sverige något sämre (33 → 48 km). En
enda global överföringsfunktion räcker alltså inte. Resonemang, litteratur
och nästa steg: se docs/HARMONIZATION.md.
