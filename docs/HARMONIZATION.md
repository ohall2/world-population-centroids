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
