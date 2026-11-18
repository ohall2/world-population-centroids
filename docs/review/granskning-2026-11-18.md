# Granskning: befolkningscentroider ur nattljus 1992–2025, påståendena P1–P3

Granskare: oberoende granskning, 2026-11-18. Underlag: repot `/home/claude/wpc` (gren v2, senaste commit `dc4573e`), lokala data under `/mnt/user-data/uploads/wpc/data/`, egna sökningar på webben. Inget i repot har ändrats. Egna skript och utdata ligger i `/home/claude/review/scripts/` och `/home/claude/review/out/`.

**Vad jag körde om:**
- `tools/check_seam_summary.py` på alla fyra `docs/check_seam/*.csv.gz`. Siffrorna i HARMONIZATION.md §12 återskapas exakt.
- `tools/check_identity.py` för SWE (VIIRS 2021, rutnät 10′). Resultatet blev D = 97,81, A = 0,70, Bb = 94,71, Bw = 3,03 och C = 0,05 km, samma som §14.
- Två egna kontroller: (i) hur skarvkvoterna beror på vilket "vanligt hopp" som används som referens och på viktningen; (ii) hur uppdelningen Bb/Bw beror på regionindelningen, med `scripts/maup_identity.py` för SWE, MEX, FIN och BRA.

---

## (a) Sammanfattande bedömning

Arbetet är noggrant dokumenterat, och koden är i sak korrekt. Uppdelningen D = A + s(Bb + Bw + C) stämmer algebraiskt, och tangentplansprojektionen är rimlig. Rastreringen och vikten ljus × cos(lat) är konsekventa. Jag hittade inga kodfel som ensamma vänder slutsatserna.

Problemet ligger i hur resultaten tolkas och vad som påstås vara nytt. Alla tre påståendena är i sin nuvarande form starkare än resultaten bär.

- **P1** bygger på en karakteristik av serierna som inte stämmer. Chiovelli m.fl. är inte "anpassad för globala summor". Den är en maskininlärningsmodell på pixelnivå som utvärderas med pixel-RMSE. PCNL i sin tur är i grunden just en bakåtskarvning förankrad i VIIRS 2013, gjord pixel för pixel. Att PCNL saknar skarv följer därför av konstruktionen. Att bakåtskarvning per enhet är nytt faller också, eftersom PCNL redan gör det per pixel.
  - Kvoten "2–3 gånger ett vanligt hopp" beror på referensen. Mot enbart DMSP-år är den 1,5–1,6 för centroid, spridning och total ljusmängd i Chiovelli. Viktat med ljusmängd är skarven i total ljusmängd i Chiovelli inte större än hoppet 2010→2011.
  - I Li DN > 7 är hoppet 2016→2017 lika stort som skarven eller större i alla fem mått.
  - "Bakåtskarvning ger bäst lägen" gäller medianer på nivå 1–2 i USA och Europa (där metoden valdes) och nivå 2 i Mexiko och Brasilien. Det gäller inte på nivå 0 i Mexiko (133 mot 35–42 km) och Brasilien (94 mot 73 km), och inte på nivå 1 i Brasilien mot Chiovelli (16,5 mot 13,8 km).
- **P2** är en identitet. Den kan inte "förklara" felet i empirisk mening: när s ≈ 1 och C ≈ 0, vilket rutnätet 10′ i praktiken garanterar, blir D = Cov/E per definition. Det som har empiriskt innehåll är att variationen ligger mellan regioner, och det beror starkt på indelningen. Med samma data går Bb från 0 till 95 km i Sverige och från 0 till 127 km i Mexiko när antalet regioner går från 1 till alla. Att felet är "förutsägbart" har inte prövats. Påståendet att felet "går åt motsatta håll i Mexiko och Norden under samma omvandling" har aldrig testats; HARMONIZATION.md rad 582 säger själv "(b) prövades inte".
- **P3** blandar ihop sensorns brus med skarvens fel. Alla förändringsperioder med facit går över skarven, eller ligger helt i DMSP-åren (USA 2000→2010). Någon förändring helt inom VIIRS har inte prövats. Kvoten "1,5–3" stämmer inte med de egna tabellerna, som ger 1,5–5,0 (Europa 4,8–5,0).

Bedömning: **större revision**. Artikeln kan bli ett bra bidrag om P1 formuleras om till en jämförelse av metoder (global överföringsfunktion mot förankring per pixel eller enhet), P2 som en uppdelning efter skala med kontrollerat skalberoende och en verklig förutsägelse, och P3 begränsas till DMSP-baserade och skarvade serier eller prövas inom VIIRS.

---

## (b) Allvarliga problem

### B1. P1 karakteriserar serierna fel, och PCNL:s kontinuitet är mekanisk

- LITERATURE.md rad 12–17 och HARMONIZATION.md rad 810–813 kallar Li m.fl. 2020 och Chiovelli m.fl. "anpassade för summor" och PCNL "anpassad pixel för pixel".
  - **Chiovelli m.fl.** omvandlar VIIRS till DMSP-likt med en ensemble av "extremely randomized trees" på 30″-pixlar, med grannskapsstatistik som förklaringsvariabler. Modellen tränas på 2013 och utvärderas med RMSE per pixel: 1,50 mot 3,27 för Li m.fl. ([Chiovelli m.fl., arbetspapper GWU 2025-009](https://cer.columbian.gwu.edu/sites/g/files/zaxdzs4246/files/2025-11/2025-009.pdf)). Den är alltså anpassad och utvärderad på pixelnivå, inte för summor.
  - HARMONIZATION.md rad 52–56 säger att Chiovelli använder "korskalibrering enligt Li och Zhou (2017)". I arbetspappret nämns Li och Zhou i beskrivningen av Li m.fl. 2020. Uppgiften måste kontrolleras.
- **PCNL** (Li, S. m.fl. 2023):
  - Metoden är regression ln-linjär på 2013, följd av en *skillnadsraster för 2013 (regresserad VNL − CCNL) som läggs på alla DMSP-år 1992–2013*. VNL används från 2014 ([MDPI](https://mdpi.com/2072-4292/15/16/3925)).
  - Det är samma princip som projektets `spliced-back`, men per pixel. PCNL 2013 är därför i stort sett regresserad VIIRS 2013, och hoppet 2013→2014 i PCNL är ett hopp inom VIIRS. Att skarven "nästan är borta" (HARMONIZATION.md rad 807) följer av konstruktionen och visar inte att harmoniseringen är bättre.
  - Följd 1: Den verkliga skillnaden mellan serierna är *global överföringsfunktion* (Li m.fl.: sigmoid; Chiovelli m.fl.: ML-modell) mot *förankring lokalt i överlappsåret* (PCNL; projektets bakåtskarvning). Det är en annan och mer begriplig tes än "summor mot pixlar".
  - Följd 2: LITERATURE.md §9 (rad 125–132) räknar "att bakåtskarvning förankrad i VIIRS ger bäst lägen" som nytt, men PCNL har redan gjort bakåtskarvning förankrad i VIIRS (per pixel). Det som eventuellt är nytt är att metoden utvärderas mot folkräkningar.
  - Rad 83 i LITERATURE.md beskriver PCNL utan att nämna att metoden är en bakåtskarvning. Det måste rättas.
- Zenodo-posten som laddats ned (`data/raw/pcnl/PCNL2010.json`) heter "PCNL, 1992–2024". Den kan vara en senare version än den som beskrivs i artikeln (1992–2021). Versionen måste anges, och metoden kontrolleras för just den versionen.

### B2. Skarvens kvot beror på referensen, och skarven är inte unik

`tools/check_seam_summary.py` rad 6 och 37–38 jämför skarven med medianen av sex hopp: tre DMSP-hopp (2010–2013) och tre VIIRS-hopp (2014–2017). VIIRS-hoppen är ungefär hälften så stora. Kvoten blir därför ett blandat mått. Mina omräkningar på `docs/check_seam/*.csv.gz` (medianer över enheter, minst 10 tända pixlar):

| Serie | Mått | Kvot mot alla 6 (redovisad) | Mot DMSP-år | Mot VIIRS-år | Mot √((D²+V²)/2)* |
| --- | --- | --- | --- | --- | --- |
| Chiovelli | centroid | 2,17 | **1,51** | 2,93 | 1,90 |
| Chiovelli | spridning | 1,95 | **1,51** | 3,15 | 1,93 |
| Chiovelli | total ljusmängd | 1,81 | **1,59** | 2,78 | 1,95 |
| Li DN > 7 | centroid | 1,82 | **1,57** | 2,17 | 1,80 |
| Li DN > 7 | total ljusmängd | 2,54 | 2,24 | 2,92 | 2,51 |
| PCNL | centroid | 1,30 | 0,90 | 2,10 | 1,17 |

\* Det väntade hoppet när det ena året har DMSP-brus och det andra VIIRS-brus, om bruset är oberoende mellan åren.

Ytterligare iakttagelser:

1. **Hoppet 2016→2017 i Li DN > 7 är lika stort som skarven eller större i alla fem mått.** Kvot mot medianen av övriga år: centroid 1,92 mot 1,82, spridning 2,23 mot 2,04, Neff 4,24 mot 2,90, tänd yta 4,24 mot 3,02, ljusmängd 2,89 mot 2,54. Hoppet 2011→2012 i Li DN > 7 är nästan lika stort som skarven för Neff, tänd yta och ljusmängd (2,35 / 2,34 / 2,19). HARMONIZATION.md rad 771–789 noterar hoppet 2016–2017 men drar inte slutsatsen att skarven 2013–2014 i Li DN > 7 inte sticker ut.
2. **Viktning.** Medianen över enheter (rad 34) ger en liten, svagt upplyst enhet samma vikt som en storstadsregion. Viktat med enhetens ljusmängd 2013 blir medelvärdet av |Δlog T| i Chiovelli 0,169 vid skarven mot 0,179 för 2010→2011 och 0,142 för 2011→2012. Med den viktningen hoppar alltså inte total ljusmängd mer än ett vanligt DMSP-år. Påståendet "alla mått ... även total ljusmängd" (rad 758–760) gäller bara för medianen över enheter med lika vikt.
3. **"Tänd yta"** är i stort sett en tröskeleffekt. I Chiovelli bloomtopcode_fix (25–60° N, 10° V–40° Ö) stiger andelen pixlar > 0 från 0,25 (2013) till 0,39 (2014), medan andelen > 7 är oförändrad (0,091 / 0,088). Måttet säger något om tröskeln för svagt ljus, inte om att serien är inkonsistent per enhet.
4. Med sex referenshopp och ingen osäkerhetsskattning är "2–3 gånger" inte belagt statistiskt. Underlaget bör vara alla år, 1992–2021 eller 2025.

**Krav:** Redovisa kvoten mot DMSP-år, mot VIIRS-år och mot det väntade hoppet för korsade sensorer. Visa fördelningen av alla årshopp 1992–2021 (inte bara 2010–2017), och ange var 2013→2014 och 2016→2017 hamnar. Redovisa både viktat med ljus och med lika vikt. Formulera P1 utifrån det.

### B3. "Bakåtskarvning ger bäst lägen" gäller inte generellt, och prövningen är delvis cirkulär

Ur `docs/validation/*_position.csv` (Barmore, median i km):

| Facit, nivå | release (bakåt) | spliced (framåt) | Chiovelli |
| --- | --- | --- | --- |
| Mexiko 0, 2010 / 2020 | **134 / 133** | 35 / 42 | 95 / 119 |
| Brasilien 0, 2022 | 94 (−53 Ö, −78 N) | **73 (+33 Ö, +65 N)** | 111 |
| Brasilien 1, 2022 | 16,5 | 26,1 | **13,8** |
| Europa 0, 2011 | 13,0 | 14,5 | **7,6** |
| USA 1, p90 2000 | 22,1 | 40,8 | **16,3** |

- **Vad testas:** under VIIRS-åren (2020, 2021, 2022) är "release" rå VIIRS och "spliced" DMSP-lik VIIRS. Jämförelsen säger där bara att rå VIIRS slår utsmetad VIIRS, inte något om skarvningen. Skarvningens riktning prövas bara under DMSP-år: USA 2000/2010, Europa 2011 och Mexiko 2010. Brasilien har bara 2022 och prövar inte skarvningen alls.
- **Cirkularitet:** riktningen valdes på USA och Europa (HARMONIZATION.md rad 530–535). Mexiko och Brasilien är de enda oberoende facit, och där är bakåtskarvning sämst på nivå 0.
- **Förändring:** bakåtskarvning är sämre än framåtskarvning för förändring i Europa (nivå 1: 1,91 mot 1,66; nivå 2: 0,63 mot 0,56) och Mexiko (nivå 1: 4,74 mot 3,08) (`*_change.csv`).
- **Statistik:** inga osäkerhetsintervall och inga parvisa test. Skillnader som 0,949 mot 0,982 km (Europa nivå 2, 2011) är inte belagda.

**Krav:** Begränsa påståendet till nivå 1–2 under DMSP-år. Redovisa parvisa skillnader per enhet med bootstrap som tar hänsyn till rumsligt beroende (block per land eller region), och redovisa nivå 0 som ett uttryckligt undantag.

### B4. P2: identiteten är en tautologi, och mellan/inom-uppdelningen beror på skalan

- **Tautologin.** `check_identity.py` rad 159–178: med s = andelen ljus i befolkade celler och C = skillnaden i läge inom celler gäller D = A + s(Cov_q(u, r)/E(r) + C) exakt för varje datamängd. På 10′ är s ≈ 1 och C ≈ 0 nästan per konstruktion, eftersom stora celler gör att allt ljus hamnar i befolkade celler. Då *måste* D ≈ Cov/E. Att identiteten "förklarar" 95–146 av 98–149 km (LITERATURE.md rad 137–140) är alltså inget empiriskt resultat. Det är samma sak som att säga att ljuscentroiden ligger fel för att ljuset inte är proportionellt mot befolkningen.
- **Det enda med empiriskt innehåll är fördelningen Bb/Bw** (variansuppdelning, law of total covariance). Den beror på regionindelningen. Mina körningar (`scripts/maup_identity.py`, rutnät 10′, GADM nivå 1 sammanslagen med k-means till k grupper; Bb + Bw är konstant):

| Land | k = 1 | 2 | 3 | 5 | 8 | ½·n | n (GADM 1) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sverige (D = 97,8) | Bb 0 / Bw 97 | 30 / 70 | 80 / 18 | 79 / 18 | 77 / 20 | 87 / 10 | **95 / 3** (21) |
| Mexiko (D = 131) | 0 / 131 | 84 / 58 | 98 / 38 | 118 / 18 | 117 / 20 | 127 / 8 | **127 / 5** (32) |
| Brasilien (D = 94,6) | 0 / 95 | 67 / 30 | 56 / 44 | 81 / 24 | 77 / 23 | 87 / 14 | **92 / 13** (27) |
| Finland (D = 71,3) | 0 / 71 | 13 / 60 | 29 / 43 | **38 / 34** (5) | – | – | – |

  - Att "det mesta ligger mellan regioner" följer delvis av att ett rumsligt autokorrelerat r med 20–30 regioner alltid ger en stor andel mellan regionerna.
  - Finland avviker inte i grunden, eftersom alla länder ger stor Bw med få regioner. Men vid samma antal (k = 5) har Finland lägre andel mellan regioner (53 %) än Sverige (81 %) och Mexiko (90 %). Finland visar alltså också något verkligt: en större del av gradienten i ljus per invånare ligger inom de fem stora regionerna.
  - Förklaringen i HARMONIZATION.md rad 924–926 är riktig men bör generaliseras: P2:s påstående om "mellan regioner" är ett påstående om skalan för rumslig variation i r, och det kräver en analys av flera skalor.
- **Skalberoendet i A/Bb** (HARMONIZATION.md rad 897–918) redovisas ärligt, men det behandlas som en egenhet hos INEGI-facit. På 1′ står A för 115 av 131 km i Mexiko och 78 km i Norge. Valet av 10′ bestämmer alltså vad som syns som "ljus per invånare". Mexikos A på fint rutnät beror på att INEGI anger stora städer som en punkt, vilket är ett fel i facit, inte i ljuset. Det talar för att använda ett annat facit för Mexiko på fin skala, inte för att välja en grov skala.
- **"Förutsägbart" har inte prövats.** r_g beräknas ur samma ljus och samma facit som D. För att P2 ska säga något empiriskt behövs en prövning utanför urvalet: skatta log r_g ur yttre data (BNP eller inkomst per invånare, täthet, urbaniseringsgrad, latitud) på några länder, förutsäg D i andra, och jämför. Utan det är "felet är förutsägbart" (LITERATURE.md rad 135) ogrundat.
- **"Motsatta håll under samma omvandling" har inte prövats** (HARMONIZATION.md rad 582: "(b) prövades inte"). LITERATURE.md §9 rad 135–136 räknar det ändå som ett nytt resultat. Det måste strykas eller prövas.
- **P2 gäller VIIRS-radians, inte "ljus" i allmänhet.** I Brasilien pekar felet åt motsatta håll med DMSP-lik VIIRS (+33 Ö, +65 N) och rå VIIRS (−53 Ö, −78 N), och i Mexiko är felet 35–42 km mot 133 km. Identiteten "förklarar" båda, vilket visar att den inte gör några anspråk. Mekanismen (inkomst och täthet) måste knytas till den radiometriska omvandlingen.
- **Kvantifiering:** "Bb står för X av D km" jämför längder av vektorer som inte summeras (rad 1012). Projicera Bb på D̂ (Bb·D̂/|D|) och redovisa resten. I Sverige ger det 94,7 av 97,8, så det spelar liten roll där, men i Brasilien (Bw 13 km) bör det redovisas.

### B5. P3: brusgolvet är inte skilt från skarvens och DMSP:s fel, och siffrorna stämmer inte

- Alla förändringsperioder i `docs/validation/*_change.csv` har minst ett år under DMSP-perioden. I den publicerade serien betyder det DMSP + förskjutning från 2013. Felet i förändringen innehåller därför (i) DMSP-årets brus, (ii) skarvens fel (bruset i DMSP 2013 och VIIRS 2013, som enligt HARMONIZATION.md rad 443–444 är 1,0 / 0,40 km på nivå 1 / 2) och (iii) VIIRS-årets brus.
  - Skarvens fel ensamt är av samma storlek som de uppmätta felen i förändring: USA nivå 2 0,74–0,80 km; Europa nivå 2 0,63 km.
  - USA 2000→2010 ligger helt inom DMSP och ger kvoten 1,8–2,0.
  - Någon förändring helt inom VIIRS har inte prövats. Påståendet "ljuset har ett brusgolv" bygger alltså på DMSP och skarven, inte på ljus i allmänhet.
- **Siffrorna:** P3 i uppdraget säger "1,5–3 gånger". Ur de egna filerna: USA nivå 1: 1,5–2,1; USA nivå 2: 1,9–3,1; Mexiko nivå 1: 2,3; Mexiko nivå 2: 3,4; **Europa nivå 1: 4,8; Europa nivå 2: 5,0.** HARMONIZATION.md rad 683–684 säger "2–5", LITERATURE.md rad 20 "tio år". Det måste göras enhetligt, och intervallet blir 1,5–5.
- **Måttet:** median |e| delat med median |sann flytt| (`validate_series.py` rad 141–147) är inget detektionstest. Korrelationen (0,01–0,34 i USA och Mexiko, ≈0 eller negativ i Europa) är mer informativ och bör vara huvudmåttet tillsammans med regressionslutning och andel enheter med rätt riktning.
- **"Går inte att mäta"** gäller enstaka år. Medel över 3–5 år i var ände minskar bruset, och det har inte prövats.
- **Facits egen förändring:** GEOSTAT 2011 (bara METHD_CL = "A", 27 länder) och 2021 (v3, 33 länder, "unallocated" borttaget) bygger på olika metoder. En "sann flytt" på 0,13 km (Europa nivå 2) kan vara lika stor som facits egna metodskillnader. Detsamma gäller USA:s kvartersgrupper, som ritas om varje folkräkning, och INEGI:s orter som tillkommer eller försvinner. Facit måste beräknas på samma underlag båda åren, eller så måste osäkerheten skattas.

### B6. Litteraturen: delar av nyhetspåståendena är redan gjorda eller fel underbyggda

Se (d). Kort sagt:
- PCNL innehåller redan bakåtskarvning per pixel och rapporterar TSOL-hoppet i Li-serien 2013–2014 nationellt och för städer. Det motsäger LITERATURE.md rad 15–16 ("aggregerad eller pixelvis utvärdering avslöjar det inte").
- Chiovelli m.fl. visar att Li-serien passar sämre per pixel (RMSE 3,27 mot 1,50). Även pixelvis utvärdering pekar alltså ut Li.
- Att nattljus fångar förändring dåligt är väl belagt, också för befolkning: Dong m.fl. 2021 och Zhang & Gibson 2025.

---

## (c) Mindre problem

1. **`check_seam_measures.py` rad 47–48:** fackelmasken (2 km buffert) bygger på VIIRS och används också för DMSP, som smetar ut ljus över flera km. Rester av fackelljus i DMSP-åren kan ge skarvhopp i fackelregioner. Det påverkar inte medianen nämnvärt, men bör prövas med större buffert för DMSP-åren.
2. **`check_seam_measures.py` rad 21–23:** förenklingen 0,002° med `preserve_topology=False` kan ge överlapp. `rasterize` låter då den sist ritade polygonen vinna. Felet är detsamma alla år och påverkar inte hoppen, men bör nämnas.
3. **`check_seam_summary.py` rad 19:** kravet på minst 10 tända pixlar *alla år* ger olika enhetsurval för olika serier (3 004–3 336 enheter). Jämförelsen mellan serier bör göras på de enheter som finns i alla.
4. **2013–2014 som sensorskarv:** det är bekräftat för Li (filnamnen `calDMSP` 2013 och `simVIIRS` 2014), Chiovelli ("predict DMSP-like values for 2014 onward") och PCNL (VNL från 2014). Men LITERATURE.md rad 40 beskriver Chiovellis kontroll "ingen diskontinuitet 2012–2013" som om den gällde skarven. I själva verket är det en kontroll utanför urvalet i överlappsåret. Skarven 2013→2014 har Chiovelli m.fl. bara prövat på landsnivå (elektrifiering i sex länder).
5. **Verklig förändring 2013–2014:** den alternativa förklaringen avfärdas inte uttryckligen i §12, men kan avfärdas med data som redan finns. Rå VNL v2 ger ett vanligt hopp 2013→2014 (0,77 mot 0,68 km på nivå 1, HARMONIZATION.md rad 341–344). Det bör visas för alla fem måtten med rå VIIRS 2013/2014 och rå DMSP 2012/2013.
6. **Hoppet 2016→2017** finns i den publicerade serien (1,4 gånger ett VIIRS-år; HARMONIZATION.md rad 779–786). Det är en inhomogenitet i VNL v2 som påverkar både P1 (referensen) och P3 (förändring 2010→2020 över 2017) och måste utredas eller redovisas.
7. **`check_identity.py` rad 145–149:** en cell får regionen där flest invånare bor. På 10′–30′ hamnar gränsceller i fel region, vilket flyttar en del från Bb till Bw. Det är litet, men bör kontrolleras genom att dela cellerna längs regiongränserna.
8. **`check_identity.py` rad 91–92:** Norges ljusfönster begränsas till 32° Ö och 72° N, men befolkningen begränsas inte på samma sätt. Det betyder inget i dag, eftersom GEOSTAT saknar Svalbard, men regeln bör vara symmetrisk.
9. **`check_identity.py` rad 182–187:** projektionen R·T/|m_P| tar bort den radiella delen av D, som speglar skillnad i spridning, och skalar med |m_P| i stället för |m_L|. Det relativa felet är 1 − |m_P|/|m_L|, ungefär 10⁻³ på landsskala, alltså försumbart. Jämförelsen med storcirkelavståndet (97,81 mot 97,87 km) bekräftar det. Inget fel, men det bör nämnas i metodbeskrivningen.
10. **`validation_centroids.py` rad 240–252 (MIN_COUNTRY_COVERAGE):**
    - Täckningen mäts som andel *enheter* med minst en facitpunkt, inte som andel befolkning eller ljus. Ett land med många små, glesa adm2-enheter kan klara 90 %, men det fångas oftast av `missing_light_share` (rad 255–272).
    - **Svaghet:** för en adm2-enhet är `missing_light_share` 0 så fort enheten har en enda punkt (rad 259). Delvis täckta adm2-enheter i ett "täckt" land räknas då som fullständiga. Det berör GEOSTAT 2011 i länder med blandade METHD_CL ("A" och "D"/"M" i samma land). Kravet på minst 10 punkter (`validate_series.py` rad 103) räcker inte för stora enheter.
    - Förslag: andelen av enhetens ljus (eller yta) i 1 km-rutor utan facit, eller att GEOSTAT 2011 begränsas till länder med 100 % "A".
    - Själva rättelsen är korrekt kodad. `_gid0` används rätt på nivå 0, och regeln `complete &= covered` tillämpas på alla nivåer.
11. **`missing_light_share`** använder DMSP 2010 (rad 284–285) som referens för facit 2020–2022. Ny bebyggelse efter 2010 bedöms inte. Det spelar sannolikt liten roll.
12. **Nivå 0 är enstaka observationer:** Mexiko n = 1, Brasilien n = 1, USA n = 2. Jämförelser på nivå 0 bör inte redovisas som medianer.
13. **Brasiliens facit** är 1 km-rutor med rutans tyngdpunkt. Det är bra nationellt, men på nivå 2 i små kommuner (och för A och C på fint rutnät) ger det ett eget fel på upp till ±0,5 km. Det bör anges som en undre gräns för felet.
14. **INEGI:** stora städer anges som en punkt per kommun (rad 7–12). Det ger ett fel i facit på nivå 2 i storstadskommuner som kan vara lika stort som felet i serien (1,0–1,1 km). Felet bör redovisas separat för kommuner med och utan stora orter.
15. **GEOSTAT 2021 "unallocated":** andelen per land redovisas i `summary.json` men tas inte upp i analysen. Om andelen skiljer sig mellan regioner (t.ex. Norrland mot Stockholm) skjuts facitcentroiden. Det bör redovisas för de nordiska länderna i §14.
16. **Namnet "Price-identiteten":** uppdelningen är elementär (viktat medel = medel + cov/medelvikt, och law of total covariance för mellan/inom). Det går att hänvisa till Price, men artikeln bör inte ge intryck av att tillämpningen av en känd identitet är ett teoretiskt bidrag.

---

## (d) Litteratur som saknas eller motsäger

**Föregångare som måste diskuteras (P1)**
- Li, S., Cao, X., Zhao, C., Jie, N., Liu, L., Chen, X. & Cui, X. (2023). Developing a pixel-scale corrected nighttime light dataset (PCNL, 1992–2021) combining DMSP-OLS and NPP-VIIRS. *Remote Sensing* 15(16):3925. doi:10.3390/rs15163925. [MDPI](https://mdpi.com/2072-4292/15/16/3925).
  - Finns redan i LITERATURE.md, men felbeskriven. Metoden: skillnadsrastern för 2013 läggs på alla DMSP-år, vilket är en bakåtskarvning per pixel.
  - Rapporterar "an abrupt rise in TSOL for LiNTL in 2013–2014" globalt, nationellt och för städer. Det motsäger att aggregerad utvärdering inte avslöjar skarven.
- Chiovelli, G., Michalopoulos, S., Papaioannou, E. & Regan, T. (2025). *Illuminating the Global South*. GWU CER Working Paper 2025-009. [PDF](https://cer.columbian.gwu.edu/sites/g/files/zaxdzs4246/files/2025-11/2025-009.pdf).
  - ML-modell per pixel, tränad på 2013 och utvärderad med pixel-RMSE.
  - Li m.fl. 2020 får RMSE 3,27 mot 1,50, och författarna konstaterar att Li-serien "comparable only for pixels with high light values". Pixelvis utvärdering pekar alltså redan ut Li.
- Cheng, H., Geng, M., Li, X., Li, S., Zhao, M., Lin, C., Wang, J., Gong, P. & Zhou, Y. (2026). *ESSD* 18:3449. [ESSD](https://essd.copernicus.org/articles/18/3449/2026/).
  - Utvärderas på pixel-, stads-, lands- och kontinentnivå, men inte per län eller kommun. I granskningen (essd-2026-129, AR1) svarade författarna att de inte lagt till regionspecifika jämförelser.
  - Det stöder att utvärdering per enhet saknas, men det är ett negativt fynd.

**Nattljus och förändring (P3), saknas i LITERATURE.md**
- Zhang, X. & Gibson, J. (2025). How well do gridded population estimates proxy for actual population changes? Evidence from China. *Population and Development Review*. doi:10.1111/padr.70025. [Population Council](https://popcouncil.org/insight/the-limits-of-gridded-data-for-tracking-population-shifts-evidence-from-china/).
  - Rutnätsdata förklarar över 95 % av variationen i läge men under 20 % av förändringen per county.
  - Arbetspappret är Waikato WP 7/24, nr 2407 ([PDF](https://repec.its.waikato.ac.nz/wai/econwp/2407.pdf)). **LITERATURE.md rad 77 tillskriver det fel innehåll**: det handlar om befolkningsrutnät, inte om ljus och ekonomisk aktivitet.
- Dong, B., Ye, Y., You, S., Zheng, Q., Huang, L., Zhu, C., Tong, C., Li, S., Li, Y. & Wang, K. (2021). Identifying and classifying shrinking cities using long-term continuous night-time light time series. *Remote Sensing* 13(16):3142. doi:10.3390/rs13163142. [MDPI](https://mdpi.com/2072-4292/13/16/3142).
  - 85,4 % av städerna minskade i befolkning 2010–2020, men bara 15,7 % i ljus. Förändringen i ljus och befolkning följs inte åt. Det är direkt relevant för P3 och för brusgolvet i snabbt föränderliga regioner.
- Gibson, J. & Boe-Gibson, G. (2021). Nighttime lights and county-level economic activity in the United States: 2001 to 2019. *Remote Sensing* 13(14):2741 (preprint doi:10.20944/preprints202105.0356.v1). [Preprint](https://www.preprints.org/manuscript/202105.0356/v1/download).
  - R² inom county ≈ 0,12 för maskad VNL v2. Författarna avstår uttryckligen från att skarva DMSP och VIIRS ("no single adjustment factor ... can be most appropriate in all times and places"). Volym och artikelnummer bör kontrolleras.
- Gibson, J. (2021). Better night lights data, for longer. *Oxford Bulletin of Economics and Statistics* 83(3):770–791. doi:10.1111/obes.12417. [EconPapers](https://econpapers.repec.org/RePEc:bla:obuest:v:83:y:2021:i:3:p:770-791).
  - DMSP har brus som återgår mot medelvärdet och underskattar rumslig ojämlikhet. Relevant för att DMSP-delen dominerar felet i P3.
- Archila Bustos, M. F., Hall, O., Niedomysl, T. & Ernstson, U. (2020). A pixel level evaluation of five multitemporal global gridded population datasets: a case study in Sweden, 1990–2015. *Population and Environment* 42(2):255–277. doi:10.1007/s11111-020-00360-8. [Springer](https://link.springer.com/article/10.1007/s11111-020-00360-8).
  - Samma forskargrupp och samma facit (SCB). Bör citeras om tidsserier och algoritmbyten.

**Ljus per invånare och skala (P2)**
- Sutton, P., Roberts, D., Elvidge, C. & Baugh, K. (2001). Census from Heaven. *International Journal of Remote Sensing* 22(16):3061–3076. [PDF](https://urizenapw02-vlp.du.edu/~paul.sutton/AAA_Sutton_WebPage/Sutton/Publications/Sut_Pub_26_Fix.pdf).
  - Sambandet mellan ljus och befolkning varierar med BNP per invånare och med boendetäthet: "a small city of a given area in a low-income country will have a higher population". Det är en tidig föregångare till P2:s mekanism, både inkomst och täthet.
- Mellander, C., Lobo, J., Stolarick, K. & Matheson, Z. (2015). Night-time light data: a good proxy measure for economic activity? *PLoS ONE* 10(10):e0139779. [EconPapers-post](https://econpapers.repec.org/RePEc:plo:pone00:0139779).
  - Svensk finskalig jämförelse av ljus med befolkning och ekonomisk aktivitet. Fulltexten kunde inte hämtas (PMC spärrades), och innehållet bör kontrolleras innan den citeras.
- Fiaschi, D., Parenti, A. & Ricci, C. (2026). Aggregation bias in proxy measurement: nighttime lights and local economic activity. arXiv:2607.14825. [arXiv](https://arxiv.org/pdf/2607.14825).
  - Visar formellt att aggregering ändrar sambandet mellan ljus och aktivitet ("pulls elasticities toward one"). Det är relevant för skalberoendet i A/Bb/Bw.
- Elvidge m.fl. (1997), *IJRS* 18(6):1373–1379, om sambandet mellan ljusets utbredning, befolkning och BNP per land, bör citeras som grund för att ljus per invånare följer inkomst. *Ej kontrollerad i källan i denna granskning.*

**Lägesfel i DMSP (P1), kompletteringar**
- Tuttle, B. T., Anderson, S. J., Elvidge, C. D., Ghosh, T., Baugh, K. & Sutton, P. C. (2014). Aladdin's Magic Lamp: active target calibration of the DMSP OLS. *Remote Sensing* 6(12):12708–12722. doi:10.3390/rs61212708. [DOAJ](https://doaj.org/article/7ed482334d954d2d987a7b7869a56c96). Samma grupp som Tuttle m.fl. 2013.

---

## (e) Kontrollerade referenser i LITERATURE.md

| Referens (rad) | Utfall | Kommentar |
| --- | --- | --- |
| Tuttle m.fl. 2013 (27) | **Stämmer; kan nu fastställas** | *PE&RS* 79(3):287–297, doi:10.14358/PERS.79.3.287 ([register](https://data.zhenzhi.com/details/1/867062780)). 2,9 km norrut och 28 punkter stämmer med avhandlingen ([CORE](https://core.ac.uk/download/217242517.pdf)). |
| Geng m.fl. 2020 (29) | Stämmer bibliografiskt; **svag källa** | Värdena och motsägelsen i tecknet för F182010 stämmer ([Walsh Medical Media](https://www.walshmedicalmedia.com/open-access/shift-estimation-for-dmspols-nighttime-light-rasters-in-20092011-53391.html)). Förlaget brukar räknas som rovförlag; artikeln bör inte bära argumentet. |
| Stathakis 2016 (28) | Titeln bekräftad indirekt | Geng m.fl. citerar titeln och anger halv pixel. Volym och sidor (13:1421–1424) kunde inte bekräftas. |
| Ash & Mazur 2020 (30) | Delvis | *Remote Sensing* 12(14):2219, doi:10.3390/rs12142219 stämmer ([MDPI](https://mdpi.com/2072-4292/12/14/2219)). Förskjutningen är en lägesförskjutning, som rad 30 antyder, men den är bekräftad för F16 2007 mot 2009, inte "under 2004–2009". Artikeln anger själv motstridiga storlekar (17 km i inledningen, ≈ 0,016° ≈ 1,8 km i tabell 2). Citera försiktigt. |
| Li m.fl. 2020 (37) | Stämmer | *Scientific Data* 7:168. |
| Chen m.fl. 2024 (39) | Stämmer bibliografiskt | doi:10.1038/s41597-024-04228-6. |
| Chiovelli m.fl. 2026 (40) | **Delvis fel i sak** | Bibliografin är rimlig (doi ej kontrollerad hos förlaget). Men tolkningen i P1 ("anpassad för summor") är fel: metoden är ML per pixel och utvärderas med pixel-RMSE. "Ingen diskontinuitet 2012–2013" avser kontrollen i överlappsåret, inte skarven. HARMONIZATION.md rad 52–56 ("korskalibrering enligt Li och Zhou 2017") måste kontrolleras. |
| Cheng m.fl. 2026 (41) | Stämmer | Författare, *ESSD* 18:3449 och doi stämmer. "Spatial misalignments" avser DMSP mot hjälpdata enligt författarsvaret, inte DMSP mot VIIRS. |
| Li, S. m.fl. 2023 PCNL (83) | **Stämmer bibliografiskt, felbeskriven metod** | Se B1. Nedladdad version heter 1992–2024. |
| Li, X. m.fl. 2026 *Environ. Monit. Assess.* 198:1287, doi:10.1007/s10661-026-16129-9 (84) | **Kunde inte verifieras** | Varken titeln eller doi:n gav träff (två sökningar). Kontrollera hos förlaget innan den citeras. Citatet "does not, by itself, establish temporal stability at the county level" bär en del av P1:s motivering. |
| Price 1970 (68) | Stämmer | *Nature* 227:520–521, doi:10.1038/227520a0. |
| Price 1972 (68) | Volym och sidor rimliga, doi okänd | Den PDF som sökningen gav var en annan Price-artikel från 1972 (*Ann. Hum. Genet.* 36:129–140). Ekv. A7 i 35:485–490 kunde inte kontrolleras. Doi:n saknas fortfarande. |
| Gaughan m.fl. 2019 (67) | **Författarna kan nu fastställas** | Gaughan, A. E., Oda, T., Sorichetta, A., Stevens, F. R., Bondarenko, M., Bun, R., Krauser, L., Yetman, G. & Nghiem, S. V. (2019). *Environ. Res. Commun.* 1:091006, doi:10.1088/2515-7620/ab3d91 ([Southampton](https://eprints.soton.ac.uk/434120/)). |
| Hall m.fl. 2019 (75) | Stämmer | *Scientific Data* 6:235. |
| Niedomysl m.fl. 2017 (76) | Stämmer bibliografiskt | *Annals AAG* 107(3):591–605 ([IDEAS](https://ideas.repec.org/a/taf/raagxx/v107y2017i3p591-605.html)). R² = 0,12 (rad 69) har inte kontrollerats. |
| "Gibson m.fl., Waikato arbetspapper 2407" (77) | **Fel** | WP 2407 är Zhang & Gibson (2024) om befolkningsrutnät i Kina, publicerad i *PDR* 2025 (doi:10.1111/padr.70025). Den handlar inte om ljus och ekonomisk aktivitet. Rätt referens för det påståendet är troligen Gibson & Boe-Gibson 2021 (*Remote Sensing* 13:2741) eller Zhang & Gibson 2022. |
| Martinez m.fl. 2023 (105) | Stämmer bibliografiskt | *Remote Sensing* 15(10):2611. |
| Bluhm & McCord 2022 (106) | Stämmer bibliografiskt | *Remote Sensing* 14(5):1190 ([MDPI](https://mdpi.com/2072-4292/14/5/1190)). |
| Pandey, Zhang & Seto 2017 (HARMONIZATION.md rad 64) | Titeln bekräftad | "Comparative evaluation of relative calibration methods for DMSP/OLS nighttime lights" ([ORNL](https://impact.ornl.gov/en/publications/comparative-evaluation-of-relative-calibration-methods-for-dmspol/)). |
| Wang m.fl. 2025, doi:10.11834/jrs.20254567 (50) | Ej kontrollerad | Doi:n ser ut som en platshållare (löpande siffror) och måste kontrolleras. |

Övriga poster (Mirza 2021, Asher 2021, Goldblatt 2020, Bennett & Smith 2017, Nieves 2021, Quah 2011, Grether & Mathys 2011, Gibson m.fl. 2021 JDE, Zhao 2015, Shen 2019, Dong 2020) har inte kontrollerats i denna granskning.

---

## (f) Föreslagna ytterligare kontroller, i prioritetsordning

1. **Skriv om P1 med rätt karakteristik av serierna** (global överföringsfunktion mot förankring lokalt i överlappsåret). Ange att PCNL är en bakåtskarvning per pixel, och att projektets bidrag i så fall är skarvning per enhet och utvärdering mot folkräkningar.
2. **Hela årsserien för skarvkontrollen:** alla hopp 1992–2021 eller 2025 för Li, Chiovelli, PCNL och den egna serien. Redovisa fördelningen och rangen för 2013→2014 och 2016→2017, kvoter mot DMSP-år, VIIRS-år och det väntade hoppet för korsade sensorer, både med lika vikt och viktat med ljus, med bootstrap över enheter (block per land).
3. **Förändring helt inom VIIRS (P3):** SCB:s rutnät per år (2013–2024), Statistics Finlands och SSB:s årliga rutnät, eller GEOSTAT 2021 mot en senare GEOSTAT-omgång. Skilj sensorns brus från skarvens fel. Pröva medel över 3–5 år i var ände.
4. **Analys av flera skalor för P2:** Bb/Bw för hierarkin nivå 1 → nivå 2 → rutnät (som i tabellen i B4), plus en jämförelse där r permuteras slumpmässigt mellan celler med bevarad rumslig autokorrelation. Redovisa vid vilken skala variationen i r ligger.
5. **Gör P2 prövbart:** förutsäg r_g ur yttre data (BNP per invånare per region, täthet, urbaniseringsgrad) i en uppsättning länder och förutsäg D i de övriga (Mexiko, Brasilien, Norden, USA och Europa på nivå 0). Utan det bör "förutsägbart" strykas.
6. **Pröva eller stryk "motsatta håll under samma omvandling":** kör identiteten med log(1 + VIIRS), VIIRS med tak och DMSP-lik VIIRS för MEX, BRA, SWE, NOR och FIN. Det ger också svar på varför Brasiliens fel byter riktning mellan rå och DMSP-lik VIIRS.
7. **Statistisk prövning av "bakåt bäst":** parvisa skillnader per enhet med blockbootstrap. Redovisa skilt för år under DMSP (prövar skarvningen) och under VIIRS (prövar bara rå mot DMSP-lik VIIRS). Redovisa nivå 0 och Brasilien nivå 1 som undantag.
8. **Facits fel:**
   - (a) GEOSTAT 2011 bara i länder med 100 % "A", eller täckning per adm2 som andel ljus i 1 km-rutor utan facit;
   - (b) samma facitunderlag båda åren vid förändring;
   - (c) INEGI-kommuner med och utan stora orter för sig, och helst ett facit på rutnät för Mexiko (INEGI:s rutnät eller AGEB);
   - (d) andelen "unallocated" i GEOSTAT 2021 per region i Norden.
9. **Utred hoppet 2016→2017 i VNL v2** (finns i alla serier som bygger på VIIRS) och redovisa dess följder för förändringen 2010→2020.
10. **Mindre kontroller:** större fackelbuffert för DMSP-åren, gemensamt enhetsurval för alla serier, celler delade längs regiongränser i `check_identity.py`, och projektion av Bb på D̂.
11. **Referenserna:** rätta rad 77 (WP 2407), kontrollera Li, X. 2026 och Wang 2025 i källan, fastställ Tuttle 2013 och Gaughan 2019 enligt (e), och kontrollera Chiovellis kalibrering av DMSP.
