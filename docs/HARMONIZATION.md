# Skarven DMSP–VIIRS: resonemang, prövade lösningar och jämförelser

Det här dokumentet följer arbetet med att foga ihop DMSP-OLS (1992–2013) och
VIIRS (2013–2025) till en serie av befolkningscentroider. Det beskriver vad
vi har prövat, varför, och vad det gav, i den ordning det gjordes. Metodernas
tekniska detaljer står i METHODS.md.

## 1. Problemet

Samma år, 2013, ger DMSP-IRQR och VIIRS VNL v2 olika centroider: median
12,5 km för länder, 3,6 km för adm1 och 1,3 km för adm2, och 150–250 km för
stora länder (Indien, Brasilien, Kanada). Ett vanligt årshopp i DMSP är
3,9 / 1,2 / 0,53 km. Att foga ihop serierna rakt av ger alltså ett falskt
hopp mellan 2013 och 2014 som är tre gånger större än ett vanligt år.

Orsaken är sensorerna. DMSP är mättat i stadskärnor (DN 63), har en
detektionströskel och smetar ut ljus över flera kilometer (blooming). VIIRS
har dubbelt så hög upplösning, ingen mättnad, ombordkalibrering och ser
svagare ljus, inte minst på landsbygden.

## 2. Första försöket: egen harmonisering VIIRS → DMSP (2026-11-16/17)

Val: att föra VIIRS till DMSP:s skala, inte tvärtom, eftersom DMSP-serien
är längst och eftersom en omvandling från grövre till finare data kräver
antaganden om detaljer som inte finns i DMSP.

Metod (wpc/harmonize.py, tools/harmonize_viirs.py): ytviktad aggregering
till 30″, gaussisk utsmetning, monoton överföringsfunktion i två steg (andel
tända pixlar och värde bland tända) anpassad på 2013, σ vald med tvåfaldig
rumslig korsvalidering. Tröskelsteget lades till för att en ren
medelvärdesfunktion skulle ge alla nästan mörka pixlar ett litet positivt
värde och dra centroiderna mot ytornas geometriska mitt.

Resultat för hela jorden: σ = 2,0 DMSP-pixlar, korsvaliderat R² = 0,90.
Skarven minskade med omkring 40 % (12,5 → 7,7 / 3,6 → 2,2 / 1,34 → 0,80 km)
men är fortfarande ungefär dubbelt så stor som ett vanligt DMSP-årshopp.
Brasilien, Ryssland och Argentina blev mycket bättre; Indien och Kina något
bättre; Kanada och Sverige sämre. Den harmoniserade serien har däremot de
minsta årshoppen av alla källor (2,3 / 0,7 / 0,25 km 2013–2025).

## 3. Hur andra har gjort

Genomgång 2026-11-17. Två huvudvägar finns.

**VIIRS → DMSP-likt.** Li m.fl. (2020, Scientific Data 7:168) är den mest
använda globala serien: DMSP kalibreras stegvis mellan åren och VIIRS görs
DMSP-likt med kärntäthet och punktspridningsfunktion, följt av en
sigmoidfunktion. De noterar själva att de simulerade DMSP-data får större
utbredning än DMSP, särskilt där ljuset är svagt, och rekommenderar att bara
använda pixlar med DN > 7. Regionala varianter finns (Zhao m.fl. 2020,
Sydostasien, sigmoid; Zheng m.fl. 2019, Kina, geografiskt viktad regression).
Chiovelli, Michalopoulos, Papaioannou och Regan (2026, Economic Journal,
doi:10.1093/ej/ueaf134) har
en global serie 1992–2025 i 30″ med korskalibrering enligt Li och Zhou
(2017), korrigering av mättnad med trunkerad Pareto (Bluhm och Krause 2022)
och modellering av blooming; den är validerad mot ekonomisk utveckling.

**DMSP → VIIRS-likt.** Nyare arbeten använder djupinlärning: Chen m.fl.
(2021, ESSD; autoencoder, 2000–2018), Nechaev m.fl. (2021, residualt U-Net),
Chen m.fl. (2024, Scientific Data 11:1380; NTLSRU-Net med NDVI, global
15″-serie 1992–2023, SVNL) och senare CUT- och diffusionsmodeller.

**Globala eller regionala modeller?** En systematisk jämförelse av
kalibreringsmetoder för DMSP (Pandey, Zhang och Seto 2017, Remote Sensing of
Environment 195:67–78) fann att globala
metoder fungerar bättre än regionala. Den gällde kalibrering mellan
DMSP-år, inte mellan sensorer, men är ett skäl att vara försiktig med
regionala överföringsfunktioner, som kan ge nya skarvar vid regiongränserna.

**Vad ingen har gjort.** Vi har inte hittat någon som utvärderar
harmoniseringen med centroider eller andra lägesmått. Det vanliga är summor
per pixel, stad eller region, eller samband med ekonomiska mått. Vårt mått
– skarven 2013 → 2014 jämfört med ett vanligt årshopp – är strängare,
eftersom ett systematiskt fel i svagt landsbygdsljus flyttar en centroid även
när summan stämmer.

## 4. Nästa steg: pröva de publicerade serierna (beslut 2026-11-17)

Innan vi bygger vidare på vår egen metod körs tre publicerade serier genom
samma pipeline, som objektivt jämförelseunderlag:

| Serie | Riktning | År | Upplösning | Källa |
| --- | --- | --- | --- | --- |
| Li m.fl. 2020, v7 | VIIRS → DMSP-likt | 1992–2021 | 30″ | Earth Engine `projects/sat-io/open-datasets/Harmonized_NTL/{dmsp,viirs}` (figshare doi:10.6084/m9.figshare.9828827), CC BY 4.0 |
| Chiovelli m.fl. 2026 | VIIRS → DMSP-likt, med korrigering av mättnad och blooming | 1992–2025 | 30″ | Zenodo doi:10.5281/zenodo.20558766, CC BY 4.0 |
| Chen m.fl. 2024, SVNL | DMSP → VIIRS-likt | 1992–2023 | 15″ | Earth Engine `projects/sat-io/open-datasets/srunet-npp-viirs-ntl` (figshare doi:10.6084/m9.figshare.22262545), licens att kontrollera |

Måtten är desamma för alla: skarvens storlek (hoppet vid sensorbytet jämfört
med vanliga årshopp före och efter) och årshoppens storlek över hela serien.
För SVNL ligger skarven mellan 2011 (simulerad) och 2012 (verklig VIIRS).

Möjliga följdsteg beroende på resultatet: (a) använda en publicerad serie
som huvud- eller jämförelseserie; (b) förbättra vår metod för svagt ljus,
t.ex. en tröskel som beror på lokal ljusmiljö; (c) regionala funktioner, med
korsvalidering och kontroll av skarvar vid regiongränserna.

## 5. Iakttagelser vid nedladdningen (2026-11-17)

**Li m.fl. 2020.** DMSP-delen (1992–2013) och den DMSP-lika VIIRS-delen
(2013–2021) är uint8 i 30″. VIIRS-delen innehåller även 2013, så serien har
ett överlappsår. För 2012 och 2013 är Li-filernas pixlar identiska med de
råa F18-kompositerna (kontrollerat för 25–60° N); andra år, t.ex. 2009, är
kalibrerade (10 % av pixlarna ändrade, högst ±6 DN). 2009–2011 har hörnet
(−180; 75) i Earth Engine, övriga år (−180,0042; 75,0042), dvs. en halv
pixel förskjutet.

**Lägesfel mellan DMSP-år.** Samma hörn (−180; 75) har våra egna
DMSP-filer för F16 2009, F18 2010 och F18 2011. Pipelinen lägger alla
DMSP-år på EOG:s rutnät oavsett metadata. För att se om pixlarna verkligen
är förskjutna jämfördes korrelationen mellan på varandra följande år
(25–60° N, 30° W–120° E) när det ena året flyttas en pixel i varje
riktning. Toppen ligger ofta snett, och olika för olika år:

| Par | −1 kol | 0 | +1 kol | −1 rad | +1 rad |
| --- | --- | --- | --- | --- | --- |
| 2009 F16 mot 2008 F16 | 0,973 | 0,970 | 0,955 | 0,944 | 0,972 |
| 2010 F18 mot 2009 F16 | 0,949 | 0,949 | 0,938 | 0,945 | 0,931 |
| 2011 F18 mot 2010 F18 | 0,923 | 0,949 | 0,965 | 0,935 | 0,940 |
| 2012 F18 mot 2011 F18 | 0,972 | 0,970 | 0,955 | 0,978 | 0,937 |

Kompositernas geografiska läge varierar alltså med omkring en halv pixel
(ca 0,5 km) mellan åren, även mellan år med samma metadata (2010 och 2011).
Det är lika stort som det vanliga årshoppet på adm2-nivå (0,53 km) och är
troligen en del av DMSP-bruset. En möjlig förbättring är att passa ihop
varje DMSP-år med ett referensår (subpixel-samregistrering) före
centroidberäkningen. Ej gjort ännu.

**Chiovelli m.fl. 2026.** Fyra varianter, 1992–2025: clean_dmsp (uint8),
topcode_fix, bloom_fix och bloomtopcode_fix (int32, upp till 2 000 efter
korrigering av mättnad). Rutnätet är 43 201 × 15 709 pixlar från (−180; 75),
alltså söderut bara till 55,9° S. Andelen tända pixlar (25–60° N) är 10,6 %
2013 och 19,9 % 2014 i clean_dmsp, men 7,5 % och 9,8 % i bloomtopcode_fix:
skarven ligger mellan 2013 och 2014, och korrigeringen för blooming
minskar ökningen av svagt ljus vid sensorbytet.

**Fel i den första jämförelsen.** dmsp-irqr har beräknats med
2009-polygonerna och viirs-dmsp med 2 km-masken (som även tar bort facklor
som har tillkommit efter 2009, t.ex. i Bakken). En del av skarven i avsnitt 2
kan alltså bero på olika masker. Alla serier jämförs därför nu med samma
mask, `flares_viirs_r2km.gpkg`, i `data/results_mask2km/`.

## 6. Jämförelse med samma fackelmask (2026-11-17)

Alla serier beräknade med `flares_viirs_r2km.gpkg` (GADM 4.1, Barmore).
Siffrorna räknas fram med

    python tools/compare_series.py --levels 0 1 --results data/results_mask2km data/results \
      --series ours=dmsp-irqr:1992-2013+viirs-dmsp:2014-2025 \
               li=li2020-dmsp:1992-2013+li2020-viirs:2014-2021 \
               chiov_cln=chiovelli-cln:1992-2025 chiov_bltc=chiovelli-bltc:1992-2025 \
      --same-year dmsp-irqr,viirs-dmsp,2013 li2020-dmsp,li2020-viirs,2013 --out docs/comparison

och ligger i docs/comparison/. Tabellen visar medianen av årshoppen (km):
"DMSP-år" är medianen av årsmedianerna 1993–2013, "skarv" är hoppet
2013 → 2014 och "VIIRS-år" medianen av årsmedianerna från 2015.

| Serie | Nivå | DMSP-år | Skarv | VIIRS-år | Skarv / DMSP-år | Skarv, 90 % |
| --- | --- | --- | --- | --- | --- | --- |
| Vår (dmsp-irqr + viirs-dmsp) | adm0 | 4,16 | 7,94 | 2,24 | 1,9 | 42 |
| | adm1 | 1,24 | 2,19 | 0,67 | 1,8 | 12,6 |
| Li m.fl. 2020 | adm0 | 4,61 | 18,90 | 4,67 | 4,1 | 187 |
| | adm1 | 1,40 | 4,75 | 1,17 | 3,4 | 37,0 |
| Chiovelli clean_dmsp | adm0 | 4,68 | 7,45 | 1,92 | 1,6 | 49 |
| | adm1 | 1,43 | 2,40 | 0,62 | 1,7 | 13,8 |
| Chiovelli bloomtopcode_fix | adm0 | 4,75 | 6,71 | 2,59 | 1,4 | 42 |
| | adm1 | 1,43 | 2,20 | 0,70 | 1,5 | 11,5 |

Slutsatser:

1. **Maskfelet i avsnitt 2 betydde lite.** Med samma mask är avståndet
   DMSP–harmoniserad VIIRS 2013 fortfarande 7,9 km (adm0) och 2,1 km (adm1).
2. **Ingen serie tar bort skarven.** Hoppet vid sensorbytet är 1,4–4 gånger
   ett vanligt DMSP-årshopp i alla serier.
3. **Li m.fl. 2020 har den största skarven** (19 km för länder, 90 % upp till
   187 km) och dessutom stort brus under VIIRS-åren. Serien passar inte för
   lägesmått.
4. **Vår serie och Chiovelli bloomtopcode_fix är jämbördiga.** Skarven är
   nästan lika stor i km (adm1 2,19 mot 2,20; adm0 7,9 mot 6,7). Kvoten blir
   lägre för Chiovelli delvis för att deras DMSP-år är brusigare (4,75 mot
   4,16 km). Korrigeringen för blooming och mättnad minskar skarven något
   jämfört med clean_dmsp.
5. **VIIRS-åren är lugnare än DMSP-åren i alla serier utom Li**, ungefär
   hälften så stora hopp.

Avståndet mellan Li-seriens DMSP och dess DMSP-lika VIIRS samma år (2013):
17,5 km (adm0) och 4,7 km (adm1), mot 7,9 och 2,1 km för vår serie.

Följder för det fortsatta arbetet: skarven finns kvar med dagens metoder och
måste hanteras uttryckligen. Möjligheter:
(a) skarva per enhet med överlappsåret – lägg till skillnaden
DMSP 2013 − harmoniserad VIIRS 2013 för varje enhet på alla VIIRS-år, så att
serien blir sammanhängande; osäkerheten blir då DMSP-brusets storlek 2013;
(b) redovisa serierna var för sig och skarven som en känd osäkerhet;
(c) minska DMSP-bruset (samregistrering mellan år, avsnitt 5), vilket både
minskar årshoppen och gör skarvningen i (a) säkrare.

## 7. Samregistrering av DMSP mot VIIRS (beslut 2026-11-17)

Val efter avsnitt 6: först minska DMSP:s eget brus genom att passa ihop
kompositerna geografiskt, sedan skarva per enhet vid överlappsåret.

**Metod** (wpc/coregister.py, tools/coregister_dmsp.py). Referensen är den
harmoniserade VIIRS-bilden 2013: läget kommer enbart från VIIRS, som har god
geolokalisering, och bilden ligger på DMSP:s rutnät och är utsmetad som
DMSP. För varje komposit och ruta om 5° × 5° med minst 1 500 tända pixlar i
båda bilderna beräknas korskorrelationen av log(1 + värde) med Hann-fönster
via FFT; toppen inom ±3 pixlar förfinas med en parabel. Kompositens
förskjutning är den viktade medianen över rutorna. Kompositen flyttas efter
kalibreringen med bilinjär interpolation, som bevarar ljussumman och flyttar
tyngdpunkten exakt lika mycket som förskjutningen (kontrollerat: Australien,
väntat +0,28/+0,78 px, uppmätt +0,29/+0,72 px; skillnaden kommer från ljus
vid fönstrets kant). Källan heter `dmsp-irqr-reg`.

**Provresultat (9 kompositer).**

| Komposit | dx | dy | Spridning (MAD) | Rutor |
| --- | --- | --- | --- | --- |
| F10 1992 | −0,96 | −0,85 | 0,15 / 0,14 | 461 |
| F14 2000 | −1,02 | −0,88 | 0,08 / 0,17 | 516 |
| F15 2000 | −0,19 | −0,21 | 0,11 / 0,17 | 538 |
| F16 2008 | −0,30 | −0,85 | 0,20 / 0,27 | 547 |
| F16 2009 | −0,82 | −0,30 | 0,27 / 0,23 | 519 |
| F18 2010 | −1,44 | −0,53 | 0,18 / 0,13 | 536 |
| F18 2011 | +0,11 | −0,37 | 0,09 / 0,08 | 559 |
| F18 2012 | −0,46 | −1,12 | 0,08 / 0,09 | 563 |
| F18 2013 | −0,28 | −0,78 | 0,09 / 0,08 | 565 |

(pixlar à 30″; negativt dx = kompositen ligger för långt västerut, negativt
dy = för långt norrut.)

Tolkning:
1. Förskjutningarna är upp till 1,4 pixlar (ca 1,3 km) och väl bestämda.
2. Nästan alla kompositer ligger nordväst om VIIRS, i genomsnitt omkring en
   halv pixel i båda led. ~~Det stämmer med att DMSP-filerna egentligen hör
   till rutnätet med hörnet (−180; 75), som Earth Engine anger.~~ *Rättat
   2026-11-17, se 7.1: Earth Engine anger samma rutnät som EOG för 32 av 35
   kompositer, så rutnätet förklarar inte förskjutningen.* Ovanpå det
   varierar varje komposit.
3. Förskjutningen är i stort sett densamma i Amerika och resten av världen
   (skillnad ≤ 0,1 px), utom F16 2008 (dy −1,23 mot −0,76). En förskjutning
   per komposit räcker som första steg; rutvärdena sparas för att kunna pröva
   regionala förskjutningar.
4. F18 2013, som harmoniseringen anpassades mot, ligger 0,78 px för långt
   norrut. Det förklarar troligen en del av skarven på adm2-nivå (0,8 km).
   Harmoniseringen görs därför om med samregistrerad DMSP
   (`tools/harmonize_viirs.py --dmsp-offsets config/dmsp_offsets.csv`,
   källan `viirs-dmsp-reg`).

### 7.1 Kontroll av rutnätet (2026-11-17, inför artikeln)

Påståendet i punkt 2 ovan kontrollerades mot metadata från nedladdningen
(`data/raw/dmsp/*.json`, fältet `crs_transform` som Earth Engine returnerar):

| Kompositer | Hörn enligt Earth Engine |
| --- | --- |
| 32 av 35 | (−180,0041667; 75,0041667), samma som EOG |
| F162009, F182010, F182011 | (−180; 75) |

Earth Engine och EOG anger alltså samma rutnät för nästan alla kompositer,
och pipelinen (`DMSP_TRANSFORM` i wpc/lights.py) använder det. Den
genomsnittliga förskjutningen åt nordväst mot VIIRS kan därför inte
förklaras av rutnätet. Den är en verklig skillnad mellan DMSP-kompositernas
och VIIRS-produktens lägen. Att DMSP ligger för långt norrut är känt:
Tuttle m.fl. (2013) fann med GPS-placerade lampor en förskjutning på i medel
2,9 km nästan rakt norrut i överflygningsdata (se docs/LITERATURE.md).

De tre undantagen är ett känt problem: rastren saknar TFW-filer och är
förskjutna av okänd anledning (Geng m.fl. 2020; Stathakis 2016).
Pipelinen läser dem med EOG:s rutnät, som de övriga. Det spelar ingen roll
för huvudserien, eftersom samregistreringen mäter och korrigerar varje
komposits förskjutning empiriskt. Relativt grannkompositerna (medlet av
F162008 och F182012, samma referens som Geng m.fl.) ger vår samregistrering
x −0,43 / −1,06 / +0,49 px och y +0,69 / +0,45 / +0,61 px för F162009 /
F182010 / F182011. Geng m.fl. fick x −0,64 / −1,27 / +0,27 och y +0,55 /
±0,09 / +0,36 (positivt y = söderut). Tecknen stämmer i x för alla tre och
i y för två; för F182010 är tecknet i y motstridigt hos Geng m.fl. Det är ett
oberoende stöd för att samregistreringen fungerar.

Kommentaren i wpc/lights.py har rättats.

### Resultat för alla 35 kompositer (2026-11-17)

Förskjutningarna står i config/dmsp_offsets.csv (dx −1,44 … +0,11 px, dy
−1,38 … −0,14 px; spridning mellan rutorna 0,07–0,39 px). Alla kompositer
utom F18 2011 (dx +0,11) ligger västerut och alla norrut om VIIRS.

Harmoniseringen mot samregistrerad DMSP 2013 gav högre korsvaliderat R²
för alla σ (σ = 2: 0,900 → 0,923); vald σ är fortfarande 2,0.

Årshopp (median, km; samma mask, Barmore):

| | adm0 | adm1 | adm2 |
| --- | --- | --- | --- |
| DMSP-år, utan samregistrering | 4,16 | 1,24 | 0,48 |
| DMSP-år, samregistrerad | 4,07 | 1,17 | 0,44 |
| Skarv 2013 → 2014, utan | 7,94 | 2,19 | 0,83 |
| Skarv 2013 → 2014, samregistrerad | 7,56 | 2,05 | 0,70 |
| Avstånd DMSP–VIIRS 2013, utan | 7,89 | 2,10 | 0,78 |
| Avstånd DMSP–VIIRS 2013, samregistrerad | 7,57 | 2,00 | 0,66 |

Störst effekt på adm2 och för vissa år: hoppet till 2011 minskar från 0,75
till 0,45 km och till 2012 från 0,63 till 0,41 km. Hoppen 2009 → 2010 (från
F16 till F18) är kvar (0,59 km). Siffrorna: docs/comparison_reg/ och
docs/comparison_reg_adm2/.

**Regionala skillnader.** Rutvärdena jämfördes med kompositens globala
förskjutning i sex regioner (Amerika, Europa/Afrika, Asien/Oceanien; norr
eller söder om 35° N). För de flesta kompositer är avvikelsen liten (median
0,08 px, 90 % under 0,27 px). F15 2002–2008 och F16 2005–2009 avviker däremot
upp till 0,7 px, ofta med motsatt tecken norr och söder om 35° N – de
satelliter vars banor drev mot slutet av sin livstid. För dem skulle ett
förskjutningsfält som varierar över jorden ge bättre resultat. F18 2013, som
skarvningen utgår från, har ingen sådan struktur (spridning 0,08 px).

## 8. Skarvning per enhet (beslut 2026-11-17)

Samregistreringen minskar skarven men tar inte bort den. Kvar är den
systematiska skillnaden mellan sensorerna, som är olika för varje enhet
(t.ex. hur mycket svagt landsbygdsljus VIIRS ser). Den skarvade serien
behåller DMSP 1992–2013 och flyttar varje enhets VIIRS-centroider
2014–2025 med enhetens skillnad mellan DMSP och VIIRS 2013 (i km österut
och norrut). Serien blir då sammanhängande vid skarven; osäkerheten i
skarvningen är lika stor som bruset i DMSP 2013 (ett vanligt årshopp).

### Resultat (2026-11-17)

Skarvat: dmsp-irqr-reg 1992–2013 + viirs-dmsp-reg 2014–2025
(`python tools/splice_series.py`, källan `spliced-reg` i
data/results_mask2km/). Skarvade enheter: 245 av 263 (adm0), 3 472 av 3 662
(adm1), 43 441 av 47 217 (adm2); övriga saknar centroid i någon källa 2013.

Flyttens storlek (Barmore, km): median 7,6 / 2,0 / 0,66, 90 % 42 / 12 / 3,5
(adm0 / adm1 / adm2). Störst för önationer med få, spridda ljus (Kiribati
455 km, Mikronesien 223 km, Sankt Helena 167 km) och för stora länder
(Kanada 224, Sydsudan 137, Indien 130, Australien 96, Kina 80 km).

Årshopp i den skarvade serien (median, km):

| | 2011 | 2012 | 2013 | **2014** | 2015 | 2016 |
| --- | --- | --- | --- | --- | --- | --- |
| adm0 | 3,94 | 3,54 | 3,61 | **2,32** | 2,37 | 1,97 |
| adm1 | 1,17 | 1,12 | 1,04 | **0,77** | 0,68 | 0,66 |

Hoppet 2013 → 2014 är nu lika stort som ett vanligt VIIRS-årshopp, som det
ska vara. Siffrorna: docs/comparison_spliced/.

**Kvarstående problem.**

1. *Skarvningen bygger på ett enda DMSP-år.* Ett slumpmässigt hopp i DMSP
   2013 för en enhet följer med till alla VIIRS-år. Osäkerheten är lika stor
   som ett vanligt DMSP-årshopp (median 4 / 1,2 / 0,44 km). Den kan minskas
   genom att uppskatta Δ från utjämnade banor, t.ex. en linjär trend i DMSP
   2009–2013 och i VIIRS 2013–2017, utvärderad 2013.
2. *Fler punkter hamnar utanför sin enhet.* På adm2 ligger 1 070 skarvade
   centroider 2020 i en annan enhet, mot 792 utan skarvning (adm1: 74 mot
   57). Flytten gäller små eller oregelbundna enheter, där DMSP:s
   utsmetning drar centroiden mot en granne; flaggan anger det.
3. *Regionala lägesfel för F15 2002–2008 och F16 2005–2009* (avsnitt 7) är
   inte korrigerade.

## 9. Skarvning: hur andra har gjort, och vårt val (2026-11-17)

**Litteratur.**

*Skarvning av ekonomiska serier.* Att skala om den äldre serien med
kvoten mellan serierna i överlappsåret kallas retropolering och är det
vanliga sättet att länka nationalräkenskaper. Alternativet, interpolering,
behandlar båda serierna som riktiga vid sina referensår och fördelar
skillnaden jämnt över åren däremellan; retropolering kan ge snedvridna
tillväxttal (Prados de la Escosura 2016). De la Fuente (2014) föreslår en
blandning där korrektionen delas mellan nivå och tillväxt efter hur
bestående orsaken till skillnaden antas vara. Vår skarvning per enhet är en
retropolering. Interpolering kräver två referensår; vi har ett överlapp,
men externa referenser (befolkningens mittpunkter från US Census 2000,
2010, 2020; SCB:s rutnät) kan spela den rollen och visa vilken sensor som
ligger närmast den verkliga befolkningen (gjort i avsnitt 10).

*Homogenisering av klimatserier.* Brott när mätinstrument byts är den
vanligaste formen av inhomogenitet, och de korrigeras med relativa metoder
som delar skillnaden mot en referens i segment med olika medelvärden. Brottet
skattas alltså från flera år på var sida, inte från ett år.

*Lägesfel i DMSP.* Ash och Mazur (2020, Remote Sensing 12:2219) visar
systematiska lägesförskjutningar under enskilda satelliters livstid (F16
2004–2009) och efterlyser automatisering; Tuttle m.fl. (2013) fann ett
lägesfel på 2,9 km, och för utsläppskartor har ett fel på omkring 1,7 km
korrigerats med vektorer per stad. Vår samregistrering (avsnitt 7) är en
automatiserad, global variant med subpixel-noggrannhet.

**Val av skattning för Δ.** Med ett enda år följer DMSP-brusets och
VIIRS-brusets slumpfel i skarvåret med till alla VIIRS-år. Tre alternativ
prövades analytiskt (standardavvikelsen för Δ i enheter av VIIRS-bruset,
med DMSP-bruset dubbelt så stort, vilket stämmer med årshoppen):

| Skattning | Osäkerhet i Δ |
| --- | --- |
| Ett år (skarvåret) | 2,24 |
| Två separata linjer, 5 + 5 år, utvärderade skarvåret | (skarvårets fel väger 0,6) |
| Gemensam trend + nivåskillnad, 5 + 5 år | 1,51 |
| **Gemensam trend + nivåskillnad, DMSP 2004–2013 + VIIRS 2013–2025, viktad** | **1,01** |
| Bara nivåskillnad (inga trender), samma år | 0,69 (men snedvriden om enheten flyttar sig) |

Separata linjer hjälper lite, eftersom skarvåret ligger i kanten av båda
fönstren. En gemensam modell – en rät linje för båda källorna och en
konstant skillnad mellan dem – utnyttjar alla år, och längre fönster gör
trenden säkrare. Vi väljer den gemensamma modellen med 2004–2013 och
2013–2025 och vikter 1/σ² (σ-kvoten skattas ur medianerna av årshoppen).
Den halverar osäkerheten och förutsätter bara att enhetens centroid rör
sig ungefär jämnt över 22 år. Ett test med en syntetisk enhet med känd
nivåskillnad och ett slumpfel i DMSP 2013 bekräftar att felet i Δ blir
mindre än hälften av felet med ett år (tests/test_splice.py).

Kommando: `python tools/splice_series.py --method trend --name spliced-reg-trend`.

### Resultat och omprövning (2026-11-17)

Den gemensamma trendmodellen gav större flyttar än skarvning på ett år
(median 11,2 mot 7,6 km för adm0; 2,7 mot 2,0 km för adm1; 0,93 mot 0,66 km
för adm2) och ett hopp 2013 → 2014 på 6,2 km (adm0) och 1,7 km (adm1). En
opartisk skarv borde ge omkring 3 km respektive 0,9 km, en kombination av ett
vanligt DMSP- och VIIRS-hopp. Modellen är alltså snedvriden. Ryssland flyttas
till exempel 136 km i stället för 50.

Orsaken syns när enheternas rörelseriktning under DMSP-åren (2004–2013)
jämförs med VIIRS-åren (2013–2025): medianen av cosinus för vinkeln mellan
lutningsvektorerna är 0,13 (adm0) och −0,02 (adm1), alltså ingen samvariation.
Lutningarna är 0,8–0,9 km/år för adm0 och 0,23–0,26 km/år för adm1, men det
som ser ut som en trend över tio år är till största delen brus och
sensordrift, inte bestående rörelse. **Antagandet om en gemensam trend
förkastas.** Metoden finns kvar i koden (--method trend) men används inte.

Skillnaden i nivå med medelvärden av flera år på var sida, utan trend:

| Hopp 2013 → 2014, median (km) | adm0 | adm1 |
| --- | --- | --- |
| Ett år (skarvåret) | 2,32 | 0,77 |
| Medel av 2 år per sida (DMSP 2012–13, VIIRS 2013–14) | 1,77 | 0,68 |
| Medel av 3 år per sida | 2,72 | 0,82 |

Skillnaderna är små, och måttet gynnar skarvning på ett år, eftersom VIIRS
då flyttas så att det träffar DMSP 2013 exakt: hoppet visar bara
VIIRS-bruset. Den verkliga osäkerheten i skarven kan inte avläsas ur serien
själv. Den är av samma storlek som DMSP:s årliga brus 2013 (vanligt årshopp
3,6 / 1,0 / 0,40 km för adm0 / adm1 / adm2).

**Beslut.** Skarvning på ett år behålls som standard: den är enkel, öppen
och motsvarar vanlig retropolering. Vilken skattning av skillnaden som ligger
närmast verkligheten avgörs med extern validering mot befolkningsdata
(US Census befolkningsmittpunkter 2000, 2010 och 2020; SCB:s rutnät;
Eurostat GISCO), som kan jämföra både DMSP- och VIIRS-delen med kända lägen.

## 10. Validering mot folkräkningar (påbörjad 2026-11-17)

Syfte: avgöra hur nära befolkningens verkliga centroider nattljusets
centroider ligger, före och efter samregistrering, harmonisering och
skarvning, och om serierna fångar verklig förändring.

Facit (tools/download_validation.py, tools/validation_centroids.py):

- **USA:** U.S. Census Bureaus befolkning och befolkningsmittpunkt per
  kvartersgrupp (omkring 220 000 punkter) 2000, 2010 och 2020.
- **Europa:** Eurostat GEOSTAT, folkräkningsbefolkning per km²-ruta, 2011
  (bara rutor som bygger på aggregerade folkräkningsdata, METHD_CL = A) och
  2021 (version 3, januari 2025). JRC-rutnätet 2018 används inte, eftersom
  det är modellerat.

Befolkningscentroiderna beräknas med samma metoder som nattljusets
(Aboufadel–Austin och Barmore), med befolkningen som vikt, för alla
GADM-enheter som täcks helt.

Mått (tools/validate_series.py): (1) lägesfel samma år; (2) förändringsfel
– skillnaden mellan facits och seriens förflyttning mellan två folkräkningar,
jämfört med facits egen förflyttning. Förändringen 2010 → 2020 (USA) och
2011 → 2021 (Europa) går över sensorbytet och prövar skarvningen direkt.

### Resultat (2026-11-17)

Facit: USA 128 av 131 delstater/motsvarande och 3 146 county fullständigt
täckta; Europa 2011: 27 länder (Italien, Grekland, Kroatien, Bulgarien och
Storbritannien saknar rutor som bygger direkt på folkräkningen), 2021: 33
länder. Enheter med färre än 10 facitpunkter utesluts.

Lägesfel (Barmore, median km; full tabell i docs/validation/):

| Serie | USA adm1 2000 / 2010 / 2020 | USA adm2 2000 / 2010 / 2020 | Europa adm1 2011 / 2021 | Europa adm2 2011 / 2021 |
| --- | --- | --- | --- | --- |
| Skarvad framåt (spliced-reg) | 1,76 / 1,59 / 1,57 | 2,23 / 2,09 / 2,18 | 3,58 / 4,10 | 0,99 / 1,17 |
| **Skarvad bakåt (spliced-back)** | **0,73 / 0,76 / 0,82** | **1,51 / 1,52 / 1,46** | **2,89 / 2,51** | **0,95 / 0,73** |
| Chiovelli bloomtopcode_fix | 1,18 / 1,08 / 0,95 | 2,01 / 1,89 / 1,84 | 3,15 / 3,13 | 1,10 / 0,94 |

Systematiskt fel (medelvärde öster/norr, USA adm1): framåt −3,9 / +4,1 km
(2000), bakåt −0,4 / +0,1 km. DMSP-centroiderna ligger systematiskt för långt
nordväst om befolkningen i USA; den omvända skarvningen tar bort det.

Förändringsfel (median km) jämfört med facits egen förflyttning:

| | Sann flytt | Framåt | Bakåt | Chiovelli |
| --- | --- | --- | --- | --- |
| USA adm1 2000 → 2010 | 0,20 | 0,35 | 0,35 | 0,70 |
| USA adm1 2010 → 2020 | 0,15 | 0,24 | 0,31 | 0,82 |
| USA adm2 2010 → 2020 | 0,26 | 0,81 | 0,80 | 1,34 |
| Europa adm1 2011 → 2021 | 0,41 | 1,66 | 1,92 | 2,87 |
| Europa adm2 2011 → 2021 | 0,13 | 0,57 | 0,64 | 0,84 |

Slutsatser:

1. **Rå VIIRS ger bäst lägen.** Att göra VIIRS DMSP-likt (harmoniseringen i
   avsnitt 2) försämrar lägena: utsmetning och mättnad flyttar centroiden.
2. **DMSP:s fel per enhet är till stor del systematiskt och bestående.**
   Skillnaden mot VIIRS 2013 gäller även 2000. Den omvända skarvningen
   (`tools/splice_series.py --direction backward --viirs viirs`) för därför
   över VIIRS noggrannhet till DMSP-åren och halverar ungefär lägesfelet på
   adm1 i USA. Den är bäst i alla jämförelser av läge, också mot Chiovelli.
3. **För förändring är framåt och bakåt ungefär lika**, och båda är klart
   bättre än Chiovelli och än att foga ihop rå VIIRS och DMSP utan skarvning
   (adm2 USA 2010 → 2020: 1,99 km).
4. **Begränsning.** I USA och Europa är alla seriers fel i förändringen
   över ett decennium större än befolkningens faktiska förflyttning. Att anta
   att ingen rörelse alls skett hade gett ett mindre fel. Subnationella
   centroider ur nattljus kan alltså visa var befolkningen bor, men inte små
   förändringar över tio år i länder där befolkningen förändras långsamt.
   Om de fångar större förändringar, till exempel i snabbt växande regioner i
   Afrika och Asien, kan inte prövas med dessa facit.

**Beslut.** Huvudserien blir den bakåt skarvade: samregistrerad DMSP-IRQR
1992–2012, flyttad per enhet med skillnaden mot rå VIIRS 2013, och rå VIIRS
VNL v2 2013–2025 (`spliced-back`). Harmoniseringen VIIRS → DMSP behövs inte
för huvudserien men redovisas som jämförelse.

**Förbehåll.** Valet av skarvriktning bygger på samma facit som används
för att utvärdera det. Två oberoende facit (USA och Europa, olika
folkräkningar och olika enhetsindelning) ger samma rangordning, vilket
talar för att resultatet håller, men en oberoende kontroll i ett tredje
område (t.ex. ett land i Afrika eller Asien med folkräkning på fin nivå) vore
värdefull.

### Tredje facit: Mexiko (INEGI, 2026-11-17)

Facit: INEGI:s ortsregister (ITER) från folkräkningarna 2010 (192 247 orter,
112,3 miljoner invånare) och 2020 (189 432 orter, 126,0 miljoner). Varje ort
är en punkt med befolkning; summaraderna tas bort (orter med en eller två
bostäder finns också var för sig). Stora städer är ofta en ort per kommun,
så facit är grovt inom städer. 32 delstater och 1 950 kommuner är
fullständigt täckta.

Lägesfel (Barmore, median km):

| Serie | Mexiko adm0 2010 / 2020 | adm1 2010 / 2020 | adm2 2010 / 2020 |
| --- | --- | --- | --- |
| Skarvad framåt (spliced-reg) | 35 / 42 | 9,0 / 11,6 | 1,42 / 1,62 |
| Skarvad bakåt (spliced-back) | 134 / 133 | 8,9 / 7,4 | **1,12 / 1,03** |
| Chiovelli bloomtopcode_fix | 95 / 119 | 11,3 / 10,6 | 1,81 / 1,71 |
| Li m.fl. 2020 | 41 / 60 | 9,1 / 17,0 | 1,44 / 2,06 |

Landsnivå i alla tre facit (median km):

| Serie | USA 2000 / 2010 / 2020 | Europa 2011 / 2021 | Mexiko 2010 / 2020 |
| --- | --- | --- | --- |
| Skarvad framåt | 59 / 73 / 54 | 14,5 / 16,9 | **35 / 42** |
| Skarvad bakåt | **34 / 40 / 43** | **13,0 / 11,4** | 134 / 133 |

Den bakåt skarvade serien är bäst på adm2 i alla tre facit och på adm1 i USA,
Europa och Mexiko 2020 (lika 2010). På landsnivå är den bäst i USA och Europa
men klart sämst i Mexiko: VIIRS-centroiden ligger omkring 80 km för långt
västerut och 110 km för långt norrut.

**Tolkning.** VIIRS mäter strålning som ökar med inkomst och ekonomisk
aktivitet. I ett land med stora regionala inkomstskillnader – det rikare,
industrialiserade norra Mexiko mot det folkrika, fattigare södra – drar det
centroiden mot de rikare områdena. DMSP:s mättnad sätter ett tak för de
ljusaste områdena och dämpar effekten. I USA och Europa är skillnaderna i
ljus per invånare mellan regionerna mindre, och där vinner VIIRS genom sin
skärpa. Det är en hypotes som behöver prövas, eftersom Mexiko bara är ett land.

**Följd.** För världens medel- och låginkomstländer, där inkomstskillnaderna
mellan regioner ofta är stora, kan rå VIIRS ge skevare nationella centroider än
DMSP. Nästa steg är (a) ett facit till med stora regionala inkomstskillnader –
Brasilien (IBGE:s rutnät 2010 och 2022) – och (b) att pröva om en konkav
omvandling av VIIRS-strålningen (t.ex. log eller ett tak som efterliknar
DMSP:s mättnad) före viktningen ger bra lägen på alla nivåer.

Utfall: (b) prövades inte. En konkav omvandling dämpar de ljusaste
områdena, vilket skulle hjälpa i Mexiko, men ökar samtidigt vikten av svagt
ljus i glesbygd. Det förvärrar det nordiska felet, där ljus per invånare är
högt just i glest befolkade områden (avsnitt 11). Felen i Mexiko och Norden
går alltså åt motsatta håll under samma omvandling. Det nordiska felet
hanteras med maskningen i avsnitt 11, och Mexikos redovisas som en
begränsning. (a) Brasilien är inte gjort, se 11.4.

## 11. Maskning norr om 65° N (som i v1)

v1 maskade ljuset norr om 65° N, eftersom nattljuset där är mindre
tillförlitligt (sena solnedgångar, norrsken); det berörde 7 länder och
0,036 % av världens befolkning (Hall m.fl. 2019). Det fördes inte över till
v2 från början. Valideringen visar att det spelar stor roll: ljus per
invånare är i Norrbotten omkring tre gånger riksgenomsnittet (7,5 % av
ljuset mot 2,4 % av befolkningen 2021), och de nordiska ländernas
nationella centroider hamnar 50–150 km för långt norrut.

Uppskattning från kommunernas centroider (Aboufadel–Austin, ljusviktat,
rå VIIRS 2021 / bakåt skarvad 2011), fel mot befolkningens centroid i km:

| Land | Alla ljus 2011 / 2021 | Utan ljus norr om 65° N 2011 / 2021 |
| --- | --- | --- |
| Sverige | 105 / 98 | 24 / 33 |
| Norge | 161 / 149 | 64 / 61 |
| Finland | 58 / 71 | 26 / 48 |

Hela kedjan kördes om med `--max-lat 65` (`wpc.lights.LatitudeMask`;
samregistrerad DMSP, rå VIIRS, bakåtskarvning; data/results_m65/).

### 11.1 Resultat på nationell nivå

Bakåt skarvad serie mot GEOSTAT, fel i km (Barmore):

| Land | Utan maskning 2011 / 2021 | Med maskning 2011 / 2021 |
| --- | --- | --- |
| Sverige | 79,4 / 98,1 | 14,4 / 32,7 |
| Norge | 95,6 / 149,3 | 67,2 / 61,1 |
| Finland | 47,7 / 71,2 | 27,7 / 42,0 |
| Danmark | 10,1 / 12,0 | 10,1 / 12,0 |
| Estland | 15,7 / 15,2 | 15,7 / 15,2 |

Uppskattningen ovan höll. Länder söder om 65° N påverkas inte, och
medianerna för USA, Europa och Mexiko ändras inte på nivå 0. Norge har kvar
ett fel på omkring 60 km, eftersom ljus per invånare stiger norrut även söder
om 65° N (Trøndelag, Nordland söder om polcirkeln).

*Rättat 2026-11-17 (se §14):* med maskning ligger Norges ljuscentroid
37 km för långt västerut och 49 km för långt **söderut**, inte norrut.
Förklaringen ovan stämmer alltså inte. Maskningen tar bort ljuset från
de nordligaste fylkena men inte befolkningen där, och det drar
ljuscentroiden söderut.

### 11.2 Resultat på nivå 1 och 2

Här gör maskningen mer skada än nytta:

- Enheter som ligger helt norr om 65° N förlorar sin centroid. Det gäller till
  exempel Norrbotten, Lappland, Finnmark, Troms och Svalbard (11 fler enheter
  på nivå 1), och på nivå 2 16 svenska och 89 norska kommuner (220 fler
  enheter i världen).
- Enheter som skärs av 65° N får sin centroid flyttad till den södra delen.
  Nordland hamnar 274 km fel (12–15 km utan maskning), och Storuman, den
  sämsta svenska kommunen, 64 km fel (12–14 km utan).
- Medianfelet för de nordiska länderna sjunker bara skenbart, eftersom de
  nordligaste enheterna faller bort.

Utan maskning var de subnationella felen redan små: median 0,8–1,1 km för
svenska kommuner och 3,8–4,1 km för län. Felet uppstår bara när en enhet
spänner över den nord–sydliga gradienten i ljus per invånare.

### 11.3 Beslut

Den publicerade datamängden tar nivå 0 från körningen med maskning norr om
65° N, som i v1, och nivå 1–2 från körningen utan maskning
(`tools/build_release.py --results-adm0 data/results_m65 --results
data/results_mask2km`). Kolumnen `max_lat` anger vilken som gäller. En
följd är att ett lands centroid inte beräknas på samma ljus som dess
regioner. Centroider summeras ändå inte mellan nivåer, och varje nivå
valideras för sig.

Valideringen görs nu direkt på de publicerade filerna (`validate_series.py
--release`), serien `release` i docs/validation/. Lägesfel, median i km
(Barmore):

| Facit | Nivå 0 | Nivå 1 | Nivå 2 |
| --- | --- | --- | --- |
| USA 2000 / 2010 / 2020 | 33,3 / 40,7 / 44,4 | 0,73 / 0,76 / 0,82 | 1,51 / 1,52 / 1,46 |
| Europa 2011 / 2021 | 13,0 / 11,4 | 2,89 / 2,51 | 0,95 / 0,73 |
| Mexiko 2010 / 2020 | 134,4 / 133,0 | 8,88 / 7,38 | 1,12 / 1,03 |

Felet i förändringen över tio år är fortfarande större än befolkningens
faktiska förflyttning på nivå 1 och 2 i alla tre områdena. I Mexiko är det
till exempel 0,85 mot 0,25 km på nivå 2.

### 11.4 Kvarstående begränsningar i v2.0.0

1. **Nationella centroider i länder med stora regionala skillnader i ljus per
   invånare** kan ligga långt fel: Mexiko omkring 130 km, Norge omkring 60 km
   trots maskningen. Inget facit finns utanför Nordamerika och Europa (ett
   fjärde facit, t.ex. IBGE:s rutnät för Brasilien, har inte prövats).
2. **Förändring över tio år** kan inte mätas på nivå 1 och 2 i något av de
   tre facit: felet är omkring 2–5 gånger befolkningens faktiska förflyttning.
3. **Regionala lägesfel för F15 2002–2008 och F16 2005–2009** (upp till 0,7
   px, avsnitt 7) är inte korrigerade. Hela samregistreringen minskade
   årsbruset på adm2 med omkring 8 % (0,48 → 0,44 km); den regionala resten är mindre. Den
   påverkar DMSP-årens brus men inte skarvningen, som utgår från F18 2013.
4. **Skarvningen bygger på ett år** (2013). Ett slumpfel i DMSP- eller
   VIIRS-centroiden det året följer med till alla DMSP-år; storleken ges per
   enhet av `noise_dmsp_km` och `noise_viirs_km`.
5. **Enheter utan skarvning** (`splice_ok = 0`; 13 / 166 / 3 198 på nivå
   0 / 1 / 2) saknar centroid i en av källorna 2013, och deras DMSP-år är
   inte jämförbara med VIIRS-åren.

Mexikos nationella fel (avsnitt 10) har en annan orsak – regionala
skillnader i ljus per invånare som följer inkomst – och redovisas som en
begränsning. Det tidigare arbetet (Niedomysl m.fl. 2017) prövade felets
samband med enheternas yta, form och befolkningstäthet i Sverige (lokala
arbetsmarknader, n = 78): inget samband med form eller täthet och bara ett
mycket svagt med ytan, som författarna förklarade med sex mycket stora och
glest befolkade regioner i norra Sverige. *(Rättat 2026-11-18 efter läsning
av fulltexten, s. 598–599; uppgiften R² = 0,12 finns inte i artikeln.)*

## 12. Gäller skarven bara centroider? (kontroll inför artikeln, 2026-11-17)

**Fråga.** Är skarven vid sensorbytet ett problem för just centroider, eller
för alla mått på hur ljuset fördelar sig inom en enhet? Och klarar de
harmoniserade serierna åtminstone summorna?

**Metod.** Ett fristående skript, med egen kod och utan pipelinen
(`tools/check_seam_measures.py`, sammanfattning med
`tools/check_seam_summary.py`). Enheterna är GADM 4.1 nivå 1, förenklad med
0,002°, och pixlar räknas med mittpunkten i polygonen. Fackelmasken har
2 km buffert, och vikten är ljus × cos(latitud). För varje enhet och år
beräknas:
- centroid (Aboufadel–Austin);
- standardavstånd på sfären, R·√(2(1 − |m|)), där m är medelvektorn;
- koncentration som effektivt antal tända pixlar, Neff = exp(entropin av
  vikterna);
- tänd yta (antal pixlar med ljus);
- total ljusmängd.

Årshoppen 2010–2017 jämförs: skarven 2013 → 2014 mot medianen av de övriga
sex hoppen. Bara enheter med minst 10 tända pixlar alla år tas med. Serier:
Li m.fl. 2020 (alla pixlar, och DN > 7 som författarna rekommenderar) och
Chiovelli m.fl. bloomtopcode_fix.

**Kontroll av skriptet.** Centroidskarven blir 4,87 km (Li) och 2,26 km
(Chiovelli). Pipelinen gav 4,75 och 2,20 km (avsnitt 6), trots annan kod,
rastrering och andra år.

**Resultat** (median över enheterna; |Δlog| för alla mått utom centroiden):

| Serie | Enheter | Mått | Vanligt hopp | Skarv | Kvot |
| --- | --- | --- | --- | --- | --- |
| Li, alla pixlar | 3 336 | Centroid (km) | 1,28 | 4,87 | 3,8 |
| | | Standardavstånd | 0,019 | 0,085 | 4,4 |
| | | Koncentration (Neff) | 0,070 | 0,686 | 9,7 |
| | | Tänd yta | 0,074 | 0,640 | 8,7 |
| | | Total ljusmängd | 0,088 | 0,534 | 6,1 |
| Li, DN > 7 | 3 004 | Centroid (km) | 1,25 | 2,27 | 1,8 |
| | | Standardavstånd | 0,022 | 0,045 | 2,0 |
| | | Koncentration (Neff) | 0,081 | 0,235 | 2,9 |
| | | Tänd yta | 0,083 | 0,251 | 3,0 |
| | | Total ljusmängd | 0,094 | 0,238 | 2,5 |
| Chiovelli bloomtopcode_fix | 3 171 | Centroid (km) | 1,04 | 2,26 | 2,2 |
| | | Standardavstånd | 0,026 | 0,051 | 2,0 |
| | | Koncentration (Neff) | 0,096 | 0,243 | 2,5 |
| | | Tänd yta | 0,080 | 0,234 | 2,9 |
| | | Total ljusmängd | 0,122 | 0,220 | 1,8 |

Kvoterna ändras lite om bara enheter med minst 100 eller 1 000 tända
pixlar tas med (centroid 1,8–2,1, total ljusmängd 1,7–3,0).

Den **globala** ljussumman hoppar inte ovanligt vid skarven. Förändringen i
log 2013 → 2014 är +0,13 (Li, DN > 7) och −0,06 (Chiovelli), mot −0,16 …
+0,29 respektive −0,06 … +0,12 för de övriga åren.

**Slutsatser.**
1. **Skarven gäller inte bara centroider.** Alla mått på enhetsnivå hoppar
   2–3 gånger mer än ett vanligt år vid sensorbytet (Li DN > 7 och
   Chiovelli). Det gäller även enhetens totala ljusmängd. Centroiden är inte
   känsligare än de andra måtten.
2. **Problemet ligger på enhetsnivå, inte globalt.** De harmoniserade
   serierna är sammanhängande i den globala summan, som de är anpassade och
   utvärderade för, men inte enhet för enhet. Formuleringen av P1 i
   artikeln ändras därför: skarven är en inkonsistens per enhet som
   aggregerad eller pixelvis utvärdering inte avslöjar. Den är inte specifik
   för lägen.
3. Li m.fl. utan tröskel har en mycket större skarv i alla mått, särskilt i
   tänd yta och koncentration. Det stämmer med författarnas egen varning för
   svagt ljus.
4. **Öppen fråga.** Hoppet 2016 → 2017 är ovanligt stort i båda serierna
   (Li DN > 7: centroid 2,41 km, global summa +0,29; Chiovelli: global summa
   +0,12). Det kan bero på en förändring i VIIRS-underlaget det året. Det
   ska kontrolleras i vår egen serie.

   *Kontrollerat 2026-11-17:* den publicerade v2.0.0-serien har också ett
   större hopp 2016 → 2017. Median för Barmore i km:

   | Nivå | 2014–15 | 2015–16 | **2016–17** | 2017–18 | 2018–19 |
   | --- | --- | --- | --- | --- | --- |
   | 0 | 2,28 | 2,31 | **3,06** | 2,16 | 2,31 |
   | 1 | 0,68 | 0,70 | **0,97** | 0,67 | 0,65 |

   Hoppet är omkring 1,4 gånger ett vanligt VIIRS-år och finns i alla tre
   serierna, som alla bygger på VIIRS. Orsaken är troligen en förändring i
   VIIRS-underlaget 2016 eller 2017. Det är inte utrett. Skarven 2013 → 2014
   i vår serie (0,80 km på nivå 1) är, som avsett, ungefär ett vanligt
   VIIRS-hopp.

**Tillägg: PCNL (Li, S. m.fl. 2023), 2026-11-17.** PCNL bygger på CCNL för
DMSP och VNL för VIIRS, och anpassningen görs pixel för pixel. Data:
Zenodo 7612389, åren 2010–2017 (`tools/download_zenodo.py --record 7612389
--out data/raw/pcnl --pattern 'PCNL201[0-7].tif'`). Samma kontroll, 3 259
enheter:

| Mått | Vanligt hopp | Skarv | Kvot |
| --- | --- | --- | --- |
| Centroid (km) | 0,72 | 0,94 | 1,3 |
| Standardavstånd | 0,013 | 0,016 | 1,25 |
| Koncentration (Neff) | 0,047 | 0,074 | 1,6 |
| Tänd yta | 0,020 | 0,084 | 4,1 |
| Total ljusmängd | 0,086 | 0,073 | 0,85 |

Årshoppen för centroiden är 1,05 / 1,24 / 0,84 km under DMSP-åren
2010–2013, 0,94 km vid skarven och 0,45 / 0,42 / 0,61 km under VIIRS-åren.

**Följd.** PCNL tar nästan bort skarven per enhet. Skarven är inte större än
ett vanligt DMSP-hopp, för alla mått utom tänd yta (kvot 4,1, troligen
tröskeln för vad som räknas som tänt). Påstående P1 måste därför
formuleras om. *(Karakteristiken i detta stycke är rättad i
docs/LITERATURE.md §11: Chiovelli m.fl. är anpassad per pixel, och PCNL är en
bakåtskarvning per pixel, så dess kontinuitet följer av konstruktionen.)*
De två serierna som är anpassade för summor (Li m.fl.;
Chiovelli m.fl.) hänger inte ihop enhet för enhet. En serie som anpassats
pixel för pixel gör det i stort sett, utom i tänd yta. Utvärderingen per
enhet skiljer alltså mellan serierna, och det är ett argument för att
använda den som standard. Hoppet 2016 → 2017 finns också i PCNL (global
summa +0,10, median |Δlog| 0,12 mot omkring 0,04).

Data: docs/check_seam/*.csv.gz (ett värde per enhet och år).

## 13. Förändring över 20 år (USA, kontroll inför artikeln, 2026-11-17)

Avsnitt 10–11 visade att felet i förändring över tio år är större än
befolkningens faktiska förflyttning på nivå 1 och 2. Gäller det också över
längre tid? I USA finns facit för 2000 och 2020, alltså 20 år och över
sensorbytet. Kommando (körs på några sekunder, kräver bara pandas):

    python tools/validate_series.py --dataset census --years 2000 2010 2020 \
        --change 2000 2020 --change 2000 2010 --change 2010 2020 \
        --release data/release/v2.0.0 --levels 1 2

Förändringsfel, Barmore, median i km (publicerade filer):

| Nivå | Period | Enheter | Verklig förflyttning | Fel i förändring | Kvot | Korrelation |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2000–2010 | 121 | 0,196 | 0,346 | 1,8 | 0,18 |
| 1 | 2010–2020 | 120 | 0,149 | 0,309 | 2,1 | 0,01 |
| 1 | **2000–2020** | 121 | **0,314** | **0,456** | **1,5** | 0,09 |
| 2 | 2000–2010 | 2 381 | 0,370 | 0,741 | 2,0 | 0,31 |
| 2 | 2010–2020 | 2 383 | 0,257 | 0,795 | 3,1 | 0,18 |
| 2 | **2000–2020** | 2 373 | **0,540** | **1,048** | **1,9** | 0,34 |

**Slutsats.** Över 20 år är felet fortfarande större än den verkliga
förflyttningen, 1,5 gånger på nivå 1 och 1,9 gånger på nivå 2. Kvoten
minskar något jämfört med tio år, eftersom den verkliga förflyttningen växer
medan felet växer långsammare. Felet ligger alltså nära ett golv. Det
bekräftar P3 för USA: förflyttningar av storleksordningen några hundra meter
till en kilometer går inte att mäta, ens över 20 år. Hur stor förflyttning
som krävs för att den ska synas kan uppskattas ur felets storlek, omkring
1 km på nivå 2. Om nattljuset fångar större förflyttningar, i snabbt
växande regioner, går inte att pröva med dessa facit.

## 14. Varför ligger landscentroiderna fel? Prövning av identiteten (2026-11-17)

**Fråga.** Förklaras de stora nationella felen (Mexiko omkring 130 km,
Norden 70–150 km) av att ljus per invånare samvarierar med läget, som
identiteten x_ljus − x_bef = Cov_P(x, r) / E_P(r) säger (Price 1972)?

**Metod** (`tools/check_identity.py`). Ljuset (rå VIIRS 2020 för Mexiko och
2021 för Norden, fackelmask 2 km, vikt ljus × cos(lat)) och facit (INEGI:s
orter 2020; GEOSTAT:s rutor 2021) läggs på ett gemensamt rutnät. Avvikelsen
D mellan ljusets och befolkningens centroid (Aboufadel–Austin, som 3D-
medelvektorer) delas då exakt i fyra delar:

- **A**: ljus i celler utan befolkning;
- **Bb**: skillnader i ljus per invånare *mellan* regioner (GADM nivå 1);
- **Bw**: skillnader i ljus per invånare *inom* regioner;
- **C**: att ljusets och befolkningens tyngdpunkt skiljer sig inom en cell.

Skriptet kontrollerar att A + Bb + Bw + C = D (exakt, i vektorform).

**Kontroll av skriptet.** Det totala felet (Aboufadel–Austin) blir
131,2 km (Mexiko), 97,9 / 32,6 km (Sverige utan / med maskning), 148,8 /
60,8 km (Norge) och 71,3 / 42,0 km (Finland). Valideringen i §10–11, med
pipelinen och den publicerade serien (Barmore), gav 133,0 / 98,1 / 32,7 /
149,3 / 61,1 / 71,2 / 42,0 km. *(Siffrorna uppdaterade 2026-11-18, efter att
facitpunkterna tilldelas land genom punkt-i-polygon i GADM nivå 1 i stället
för via VIIRS-rutnätet; skillnaden var högst 2 km.)*

**Resultat, rutnät 10′ (omkring 18 km i nord–syd).** Längder i km och, för
D och Bb, riktning (öst, nord):

| Land | Maskning | D (öst, nord) | A | **Bb** (öst, nord) | Bw | C |
| --- | --- | --- | --- | --- | --- | --- |
| Mexiko | – | 131 (−75, +107) | 0,5 | **127** (−77, +102) | 5,3 | 0,3 |
| Sverige | – | 98 (+25, +95) | 0,7 | **95** (+26, +91) | 3,0 | 0,1 |
| Sverige | 65° N | 33 (+4, +32) | 0,2 | **32** (+4, +31) | 1,0 | 0,1 |
| Norge | – | 149 (+43, +142) | 1,6 | **146** (+42, +139) | 1,5 | 0,2 |
| Norge | 65° N | 61 (−37, −49) | 0,1 | **61** (−36, −49) | 0,6 | 0,2 |
| Finland | – | 71 (−30, +65) | 0,5 | 38 (−12, +36) | 34 | 0,4 |
| Finland | 65° N | 42 (−36, +22) | 0,0 | 18 (−17, −7) | 35 | 0,4 |

Ljus per invånare i förhållande till landets medel (rutnät 10′): i Mexiko
från 0,54 (Guerrero), 0,62 (Distrito Federal) och 0,65 (Oaxaca) till
1,72 (Chihuahua) och 1,76 (Sonora). I Sverige är Norrbotten 2,33 och Skåne
0,77, i Norge är Finnmark 3,00, Troms 1,95 och Oslo 0,34, och i Finland är
Lappland 1,56 och Södra Finland 0,64 (rutnät 1′). På rutnätet 10′, som artikeln använder: Norrbotten 3,07, Finnmark 3,98, Oslo 0,29.

**Rutnätets betydelse.** Identiteten gäller exakt för varje rutnät, men
fördelningen mellan A och B beror på cellstorleken:

| Land | Rutnät | Andel ljus i befolkade celler | A | Bb | Bw |
| --- | --- | --- | --- | --- | --- |
| Mexiko | 1′ | 0,48 | 115 | 21 | 2,5 |
| Mexiko | 2′ | 0,75 | 87 | 48 | 2,9 |
| Mexiko | 5′ | 0,97 | 20 | 112 | 3,1 |
| Mexiko | 10′ | 1,00 | 0,5 | 127 | 5,3 |
| Mexiko | 30′ | 1,00 | 0,0 | 129 | 3,9 |
| Sverige | 1′ | 0,90 | 34 | 63 | 1,0 |
| Norge | 1′ | 0,81 | 78 | 71 | 0,3 |
| Finland | 1′ | 0,84 | 23 | 25 | 24 |

I Mexiko är A stor i fina rutnät, eftersom INEGI anger stora städer som en
enda punkt per kommun. Stadens ljus hamnar då i celler utan facitpunkt. Från
5′ och uppåt ligger nästan allt ljus i befolkade celler. Då förklarar
skillnaderna mellan delstaterna 112–129 av 131 km. I Norden, där facit är
rutor om 1 km², ligger 10–20 % av ljuset i rutor utan invånare (vägar,
industri, gruvor). Det ljuset ligger mest norrut, och i 1′-rutnätet står det
för 22–78 km. Men även där står skillnaderna mellan regioner för den största
delen (25–71 km).

**Slutsatser.**
1. **Identiteten förklarar felen.** De nationella felen beror nästan helt på
   att ljus per invånare samvarierar med läget. Det mesta ligger mellan
   regioner (delstater, län, fylken). Ljusets läge inom cellerna (C) betyder
   ingenting. Spridningen inom regioner (Bw) betyder lite, utom i Finland,
   där GADM:s nivå 1 bara har fem stora regioner som spänner över stora
   delar av landet.
2. **Riktningen förutsägs av vilka regioner som är ljusa.** I Mexiko lyser de
   rikare nordliga delstaterna 1,5–1,8 gånger medlet, och felet pekar mot
   nordväst. I Norden lyser de glesa nordliga länen 1,5–3 gånger medlet, och
   felet pekar norrut. Samma mekanism ger olika riktningar beroende på var
   ljus per invånare är högt.
3. **Maskningen norr om 65° N tar inte bort orsaken, den förskjuter den.** I
   Norge vänder felet: utan maskning ligger centroiden 142 km för långt
   norrut, med maskning 49 km för långt söderut, eftersom befolkningen norr
   om 65° N då saknar ljus. Påståendet i §11.1 om orsaken har rättats.
4. **Uppdelningen i "ljus utan invånare" och "olika ljus per invånare" beror
   på skalan.** Det är två sidor av samma sak: ljus som inte följer
   befolkningen. Artikeln bör redovisa identiteten på en skala där i stort
   sett allt ljus ligger i befolkade celler (10′), och visa skalberoendet
   som känslighetsanalys.

Resultat: docs/check_identity/*.json.

## 15. Brasilien som fjärde facit, och rättad täckningskontroll (2026-11-18)

### 15.1 Facit

IBGE:s Grade Estatística för Censo Demográfico 2022, rutor om 1 km² för hela
Brasilien (`BR1KM_20251002.zip`, se docs/DATA.md). 1 501 555 befolkade rutor
med 203 002 064 personer; det officiella resultatet för Censo 2022 är 203,08
miljoner. Rutans tyngdpunkt används. Kommando:

    python tools/validation_centroids.py --dataset ibge --years 2022

### 15.2 Rättad täckningskontroll

Under arbetet framkom att enheter i grannländer kunde räknas som
fullständiga, om de bara fått några facitrutor eller -punkter vid gränsen.
Kontrollen av täckning (andel av ljuset i underenheter utan befolkning)
räknade då enheten som täckt. Det gällde Paraguay och Guyana för IBGE,
Ukraina, Ryssland, Vitryssland, Turkiet, Marocko, Serbien och Bosnien för
GEOSTAT, och Bulgarien och Luxemburg 2011, där GEOSTAT bara har aggregerade
folkräkningsdata för en del av landet. **Ny regel**
(`MIN_COUNTRY_COVERAGE` i tools/validation_centroids.py): ett land räknas
som täckt först om minst 90 % av dess finaste GADM-enheter har befolkning i
facit. Antal felaktigt fullständiga nivå 2-enheter före rättelsen: 1–2
(USA), 5–6 (Mexiko), 231 (Europa 2011), 80 (Europa 2021), 55 (Brasilien);
på nivå 1 tre av 30 för Brasilien.

Valideringen av den publicerade serien kördes om med den nya regeln
(docs/validation/). USA och Mexiko ändras inte (de få enheterna föll redan
bort genom kravet på minst tio facitpunkter). För Europa:

| Nivå | År | Enheter före / efter | Median före / efter (km) |
| --- | --- | --- | --- |
| 1 | 2011 | 365 / 362 | 2,89 / 2,90 |
| 1 | 2021 | 435 / 432 | 2,51 / 2,47 |
| 2 | 2011 | 6 393 / 6 360 | 0,954 / 0,949 |
| 2 | 2021 | 7 246 / 7 234 | 0,733 / 0,732 |

Slutsatserna påverkas inte. `metadata.json` i den publicerade datamängden
innehåller valideringen och ska byggas om (`tools/build_release.py`) innan
taggen sätts.

### 15.3 Resultat för Brasilien 2022

Lägesfel, Barmore (docs/validation/ibge_position.csv):

| Serie | Nivå 0 | Nivå 1 (27 delstater) | Nivå 2 (5 567 kommuner) |
| --- | --- | --- | --- |
| Publicerad (v2.0.0) | 94,4 km (−53 öst, −78 nord) | median 16,5, p90 39,6 | median 1,11, p90 4,30 |
| Chiovelli bloomtopcode_fix | 110,9 | 13,8 | 1,53 |

Nivå 2 ligger i samma storleksordning som i USA, Europa och Mexiko
(0,7–1,5 km). Nivå 1 är sämre (16,5 km mot 0,8–7,4), men delstaterna är
mycket stora: Amazonas och Pará är 1,6 och 1,3 miljoner km², och de har de
största felen (76 och 82 km). I nordöstra Brasiliens kuststater ligger
ljuscentroiden 15–32 km för långt västerut, alltså inåt land (Bahia,
Pernambuco, Paraíba, Alagoas, Rio Grande do Norte); i São Paulo ligger den
45 km mot nordväst, bort från storstadsområdet.

### 15.4 Identiteten för Brasilien

`tools/check_identity.py --facit points` med VIIRS 2022 (docs/check_identity/id_BRA_*.json):

| Rutnät | Andel ljus i befolkade celler | D | A | Bb (mellan delstater) | Bw (inom) | C |
| --- | --- | --- | --- | --- | --- | --- |
| 10′ | 1,00 | 94,6 | 0,6 | 91,8 | 13,4 | 0,1 |
| 2′ | 0,99 | 94,6 | 3,3 | 93,5 | 12,8 | 0,1 |
| 1′ | 0,95 | 94,6 | 6,2 | 95,2 | 11,5 | 0,0 |

(Delarna är vektorer; längderna summeras inte.) Skillnaderna i ljus per
invånare mellan delstaterna förklarar felet, oberoende av rutnätet.
Ljus per invånare i förhållande till landets medel: Mato Grosso 1,69,
Distrito Federal 1,54, Sergipe 1,45, Santa Catarina 1,32, Goiás 1,29;
Amazonas 0,60, Pará 0,72, Maranhão 0,73, Acre 0,77, Pernambuco 0,77,
Rio de Janeiro 0,83, São Paulo 0,95.

**Slutsatser.**
1. Brasilien prövades efter att metoden och identiteten utvecklats på de
   andra facit, och resultatet är detsamma: det nationella felet är ett fel
   i ljus per invånare mellan regioner.
2. Ljus per invånare följer inte bara inkomst. De största storstads-
   områdena (Rio de Janeiro, São Paulo, och Distrito Federal i Mexiko med
   0,62) har lågt ljus per invånare trots hög inkomst, sannolikt för att tät
   bebyggelse ger mindre ljus per person. Glest befolkade jordbruks- och
   gruvregioner (Mato Grosso; Norrbotten, Finnmark) har högt. Mekanismen i
   P2 bör därför beskrivas som rumslig variation i ljus per invånare, med
   inkomst och bebyggelsetäthet som två kända orsaker, inte som enbart
   inkomst.

## 16. Kan nattljusets centroider mäta förflyttning? Statistisk prövning (2026-11-18)

**Fråga.** §10 och §13 visade att felet i förändring är större än den verkliga
förflyttningen (medianer). Här prövas det ordentligt, per enhet och med
osäkerhetsintervall, i den publicerade serien och i Chiovelli m.fl.
(bloomtopcode_fix) som jämförelse. Inga nya data: facit och serierna är de
som redan finns. GEOSTAT 2011 är omräknad med Bulgarien och Tjeckien
uteslutna (blandat underlag).

**Metod** (`tools/change_detection.py`). För varje enhet: verklig
förflyttning t (facit) och uppmätt e (ljus), i km öster och norr. Mått:
andel enheter där e pekar åt rätt håll (cos > 0; slump = 0,5), Spearmans
rangkorrelation och OLS-lutning mellan e och t (komponenterna sammanslagna),
och andel enheter där den verkliga förflyttningen är större än felet.
95 %-intervall med blockbootstrap (1 000 dragningar; block = nivå 1-enhet för
nivå 2, land för Europa nivå 1, annars enheten). Kommando:

    python tools/change_detection.py --dataset census --pairs 2000-2010 2010-2020 2000-2020 \
        --series chiov_bltc=data/results_mask2km/gadm41/chiovelli-bltc --out docs/change_detection
    # likadant för geostat (2011-2021) och inegi (2010-2020)

**Resultat, publicerad serie** (docs/change_detection/):

| Facit | Nivå | Period | Enheter | Verklig flytt (median km) | Rätt riktning | Rangkorrelation | Lutning | Flytt > fel |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| USA | 1 | 2000–2010 | 121 | 0,20 | 0,54 (0,45–0,63) | 0,05 (−0,12–0,22) | 0,73 (−0,42–2,04) | 26 % |
| USA | 1 | 2010–2020 | 120 | 0,15 | 0,51 (0,43–0,59) | 0,06 (−0,13–0,25) | 0,02 (−1,0–0,90) | 18 % |
| USA | 1 | 2000–2020 | 121 | 0,31 | 0,64 (0,55–0,72) | 0,24 (0,08–0,39) | 0,32 (−0,25–0,96) | 33 % |
| USA | 2 | 2000–2010 | 2 381 | 0,37 | 0,61 (0,58–0,64) | 0,24 (0,18–0,29) | 0,45 (0,30–0,60) | 26 % |
| USA | 2 | 2010–2020 | 2 383 | 0,26 | 0,58 (0,56–0,60) | 0,18 (0,14–0,21) | 0,49 (0,19–0,77) | 16 % |
| USA | 2 | 2000–2020 | 2 373 | 0,54 | 0,66 (0,64–0,69) | 0,34 (0,29–0,38) | 0,48 (0,36–0,59) | 28 % |
| Europa | 1 | 2011–2021 | 329 | 0,42 | 0,49 (0,43–0,53) | −0,08 (−0,16–−0,01) | −0,14 (−0,70–0,53) | 5 % |
| Europa | 2 | 2011–2021 | 6 161 | 0,13 | 0,51 (0,49–0,53) | 0,02 (−0,01–0,05) | 0,22 (−0,03–0,47) | 8 % |
| Mexiko | 1 | 2010–2020 | 32 | 2,07 | 0,72 (0,56–0,88) | 0,42 (0,09–0,73) | 0,38 (−0,26–2,44) | 19 % |
| Mexiko | 2 | 2010–2020 | 1 915 | 0,25 | 0,56 (0,53–0,58) | 0,13 (0,08–0,18) | 0,26 (0,02–0,86) | 14 % |

Chiovelli m.fl. ger samma bild eller sämre (rätt riktning 0,40–0,63,
rangkorrelation −0,21–0,27). Resultatet beror alltså inte på vår skarvning.

**Hur stor måste förflyttningen vara?** Andel enheter med rätt riktning, efter
den verkliga förflyttningens storlek (nivå 2, docs/change_detection/*_bins.csv):

| Facit, period | < 0,1 km | 0,1–0,25 | 0,25–0,5 | 0,5–1 | 1–2 | 2–5 | 5–10 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| USA 2000–2010 | 0,53 | 0,56 | 0,60 | 0,64 | 0,75 | 0,77 | 0,77 (13) |
| USA 2000–2020 | 0,54 | 0,58 | 0,60 | 0,69 | 0,73 | 0,84 | 0,91 (22) |
| Europa 2011–2021 | 0,51 | 0,52 | 0,50 | 0,50 | 0,51 | 0,69 (13) | – |
| Mexiko 2010–2020 | 0,50 | 0,56 | 0,57 | 0,54 | 0,67 | 0,58 | 0,43 (7) |

(antal enheter inom parentes där de är färre än 25)

**Slutsatser.**
1. **Förflyttningar av den storlek som förekommer går inte att mäta.** På
   nivå 2 pekar den uppmätta förflyttningen åt rätt håll i 51–66 % av
   enheterna (slumpen ger 50 %). Rangkorrelationen är 0,02–0,34, och i 72–92 %
   av enheterna är felet större än den verkliga förflyttningen. På nivå 1 är
   intervallen breda, och rätt riktning skiljer sig inte från slumpen i
   tre av fem fall (undantagen är USA 2000–2020 och Mexiko).
2. **Det som syns är svagt och dämpat.** Där ett samband finns (USA nivå 2,
   Mexiko nivå 1) är lutningen 0,3–0,5: ljuset fångar mindre än hälften av
   förflyttningen.
3. **En förflyttning måste vara flera km för att riktningen ska bli
   tillförlitlig.** I USA krävs över 2–5 km för att 80–90 % av enheterna ska
   få rätt riktning; sådana förflyttningar finns i under 10 % av enheterna.
   I Europa och Mexiko nås inte 80 % i någon klass.
4. **Det gäller båda sensorerna.** USA 2000–2010 ligger helt inom DMSP-tiden
   och ger samma resultat (rätt riktning 0,61 på nivå 2) som perioderna över
   skarven. Förändring helt inom VIIRS-tiden har inte prövats, eftersom facit
   saknas; slutsatsen gäller därför DMSP-perioden och ihopfogade serier.
5. **Förbehåll för facit.** Europas verkliga förflyttningar är små (median
   0,13 km på nivå 2), och GEOSTAT 2011 och 2021 bygger på olika underlag, så
   en del av "den verkliga förflyttningen" kan vara skillnader i facit.
   Europas resultat är därför den svagaste delen av belägget; USA och Mexiko
   bygger på samma typ av underlag båda åren.

*Avrundning rättad 2026-11-18: tabellen ovan avrundades tidigare uppåt på fel sätt i fyra celler; värdena följer nu CSV-filerna med avrundning halvt uppåt (tools/write_release_docs.py).*

## 17. Rättelse: USA och Puerto Rico redovisas var för sig (2026-11-18)

Valideringen mot US Census slog ihop USA och Puerto Rico. I GADM är Puerto
Ricos 78 municipios nivå 1-enheter. De är fler än de 50 delstaterna och styr
därför medianen på nivå 1. På nivå 0 blev medianen av två länder ett
medelvärde av USA och Puerto Rico. Per enhet (tools/figures/unit_errors.py,
docs/figures/data/unit_errors.csv.gz; samma urval som tidigare, medianerna
återskapas exakt):

| | 2000 | 2010 | 2020 |
| --- | --- | --- | --- |
| USA nivå 0 (km) | 65,3 | 80,4 | 87,9 |
| Puerto Rico nivå 0 (km) | 1,32 | 0,96 | 1,00 |
| USA:s delstater, median (km) | 13,2 | 10,5 | 12,1 |
| Puerto Ricos municipios, median (km) | 0,36 | 0,33 | 0,33 |
| USA:s counties, median (km) | 1,51 | 1,52 | 1,46 |

Tidigare angivna 33–44 km (nivå 0) och 0,73–0,82 km (nivå 1) för "USA" ska
inte användas. Nivå 2 påverkas inte. Testet av förändring på nivå 1 för
census (§16) domineras av Puerto Rico och redovisas i artikeln inte som ett
resultat för USA. Slutsatserna om nivå 2 står fast; landsfelet för USA
(65–88 km) ligger i samma storleksordning som för Mexiko och Brasilien.

## 18. Ändringar efter granskningen av manuset (2026-11-18)

Granskningen (docs/review/review-manuscript-2026-11-18.md) ledde till följande
ändringar utan nya körningar: figurerna 1 och 2 bytte nummer (lägesfelet först);
Fig. 6 mäter skarven i varje series eget skarvår (2012→13 för v2.0.0, kvot 0,89;
2013→14 för övriga); antalet skarvade enheter i utgåvan är 250 av 263, 3 496 av
3 661 och 44 019 av 47 217 (splice_ok = 1 år 2012); 34 DMSP-kompositer används
(F15 2008 saknar kalibrering). Manusets påståenden om förändring, landsfel och
korrektion formulerades om enligt granskningen.

## 19. Fyra kontroller efter granskningen av manuset (2026-11-18)

`tools/robustness_checks.py`, resultat i docs/robustness/. Bara befintliga data.

1. **USA:s 50 delstater (utan Puerto Rico).** Verklig förflyttning 2,6–3,9 km
   (median), men rätt riktning bara 0,50–0,58 med intervall som omfattar 0,5;
   felet i förändringen är 5,5–6,6 km. Mexikos delstater (2,1 km, fel 4,7 km)
   får 0,72. Det avgörande är förflyttningen i förhållande till lägesfelet, inte
   dess storlek i km.
2. **Treårsmedel** kring varje folkräkningsår minskar medianfelet i uppmätt
   förflyttning på nivå 2 med 6–15 %, men ändrar andelen rätt riktning med högst
   två procentenheter. USA:s delstater 2000–2020: 0,58 → 0,64 (0,52–0,76).
3. **Permutationstest** (1 000 omblandningar inom block): nollfördelningen för
   rätt riktning på nivå 2 är 0,50 (0,48–0,52) och för rangkorrelationen ±0,04.
   USA:s counties och Mexikos kommuner p ≤ 0,001; Mexikos delstater p = 0,01;
   Europas kommuner p = 0,03; USA:s delstater p = 0,19–0,48; Europas regioner,
   negativ rangkorrelation p = 0,04.
4. **Befolkningsviktat medel av regionernas ljuscentroider** som landscentroid:
   USA 65–88 → 6–11 km, Mexiko 133–134 → 1,6–2,7 km, Brasilien 94 → 12 km,
   Europa median 11–13 → 1,1–1,3 km (mindre fel i 55 av 56 land-år; Finland kvar
   på 14–27 km eftersom felet där ligger inom regionerna). Bekräftar att
   landsfelet till största delen ligger mellan regionerna, även för USA som inte
   delades upp.
