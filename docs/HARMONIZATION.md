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
Chiovelli, Michalopoulos, Papaioannou och Regan (Economic Journal 2025) har
en global serie 1992–2025 i 30″ med korskalibrering enligt Li och Zhou
(2017), korrigering av mättnad med trunkerad Pareto (Bluhm och Krause 2022)
och modellering av blooming; den är validerad mot ekonomisk utveckling.

**DMSP → VIIRS-likt.** Nyare arbeten använder djupinlärning: Chen m.fl.
(2021, ESSD; autoencoder, 2000–2018), Nechaev m.fl. (2021, residualt U-Net),
Chen m.fl. (2024, Scientific Data 11:1380; NTLSRU-Net med NDVI, global
15″-serie 1992–2023, SVNL) och senare CUT- och diffusionsmodeller.

**Globala eller regionala modeller?** En systematisk jämförelse av
kalibreringsmetoder för DMSP (Pandey, Zhang och Seto 2017, Remote Sensing of
Environment; referensen ska kontrolleras) fann att globala
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
   halv pixel i båda led. Det stämmer med att DMSP-filerna egentligen hör
   till rutnätet med hörnet (−180; 75), som Earth Engine anger, och inte till
   EOG:s rutnät (−180,0042; 75,0042), som vi har använt. Ovanpå det varierar
   varje komposit.
3. Förskjutningen är i stort sett densamma i Amerika och resten av världen
   (skillnad ≤ 0,1 px), utom F16 2008 (dy −1,23 mot −0,76). En förskjutning
   per komposit räcker som första steg; rutvärdena sparas för att kunna pröva
   regionala förskjutningar.
4. F18 2013, som harmoniseringen anpassades mot, ligger 0,78 px för långt
   norrut. Det förklarar troligen en del av skarven på adm2-nivå (0,8 km).
   Harmoniseringen görs därför om med samregistrerad DMSP
   (`tools/harmonize_viirs.py --dmsp-offsets config/dmsp_offsets.csv`,
   källan `viirs-dmsp-reg`).

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
