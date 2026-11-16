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
- Polygonerna täcker facklor som syntes fram till 2008. Nya facklor efter det fångas inte. För VIIRS-åren behövs en annan källa.

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

## Referenser

- Hall, O., Archila Bustos, M. F., Olén, N. B. & Niedomysl, T. (2019). Population centroids of the world administrative units from nighttime lights 1992–2013. *Scientific Data* 6, 235.
- Elvidge, C. D. m.fl. (2009). A fifteen year record of global natural gas flaring derived from satellite data. *Energies* 2, 595–622.
- Elvidge, C. D. m.fl. (2014). National trends in satellite-observed lighting 1992–2012. I *Global Urban Monitoring and Assessment through Earth Observation*. CRC Press.


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
