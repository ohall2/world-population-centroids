# Analysplan: utökade kontroller efter granskningen 2026-11-20

**Status: ** Planen låses genom att den committas innan någon analyskod
skrivs. Därefter ändras den bara genom en rad i avsnittet *Avvikelser* längst ned.

Bakgrund: `docs/review/granskning-2026-11-20.md`. Frågan flyttas från "kan
nattljus mäta förflyttning?" till "när, och med vilken skattning?".

## 1. Utforskande och bekräftande

- **Utforskande:** A1–A3 och B. De körs på facit som redan har använts för
  designval (USA, Europa, Mexiko, Brasilien). Resultaten redovisas som utforskande.
- **Senare moment:** D (flyttmönster och flöden) har en egen plan, se avsnitt 8.
- **Bekräftande:** C (Kanada). Facit hämtas först när skripten för A1–A3 och B är
  committade. Samma skript körs oförändrade.

## 2. Fryst

- Utgåvan `data/release/v2.0.0` och facit i `data/validation/`.
- Måtten i `tools/change_detection.py` (`metrics`): t = verklig flytt, e = uppmätt
  flytt (km öster, norr); rätt riktning (cos > 0), Spearmans ρ, median |e − t|,
  andel |t| > |e − t|.
- Urval: enheter som är `complete` med minst tio facitpunkter båda åren.
  USA på nivå 1 är delstaterna utan Puerto Rico.
- Intervall: blockbootstrap, 1 000 dragningar, block som i `block_ids`.
- Perioder: USA 2000–2010, 2010–2020, 2000–2020; Mexiko 2010–2020; Europa 2011–2021.

Nya skript läggs i `tools/`, tester i `tests/`, resultat i `docs/extended/`.

## 3. A1 – Hybridskattning

**Fråga.** Tillför ljuset något när befolkningstal per underenhet är kända?

**Skattningar** för en nivå 1-enhet och ett år:
1. *Ljus:* utgåvans ljusviktade centroid (nuvarande).
2. *Hybrid-ljus:* medel (3D) av nivå 2-enheternas ljuscentroider, viktat med
   nivå 2-enheternas befolkning i facit samma år.
3. *Hybrid-geometri:* samma vikter, men nivå 2-enheternas geometriska mittpunkter
   (ytcentroid i ytriktig projektion). Skiljer ut vad ljuset bidrar med.

Alla nivå 2-enheter med befolkning > 0 ingår. En nivå 1-enhet används om
nivå 2-enheter med ljuscentroid täcker minst 99 % av dess befolkning.

**Mått.** Läge: medianfel mot facit per facit och år, och parvis skillnad per enhet
med intervall. Förändring: samma mått som i avsnitt 2.

**Hypoteser.** H1a: hybrid-ljus har mindre lägesfel än ljus. H1b: hybrid-ljus ger
rätt riktning oftare än ljus. H1c: hybrid-ljus har mindre lägesfel än
hybrid-geometri.

**Förbehåll som ska stå i texten.** Hybriden använder facits egna befolkningstal
per nivå 2-enhet. Den prövar därför hur mycket fel som återstår från läget inom
enheterna, inte en oberoende skattning.

## 4. A2 – Trendskattning

**Fråga.** Minskar felet om hela tidsserien används i stället för två ändpunkter?

**Skattning.** Per enhet och koordinat (km öster och norr i lokalt plan kring
enhetens centroid 2013): Theil–Sen-lutning inom varje sensorperiod.
- Period helt inom DMSP (y1, y2 ≤ 2012): lutning över åren y1−2 … y2+2, högst till
  2012; flytt = lutning × (y2 − y1).
- Period över skarven: s_D över åren y1−2 … 2012 och s_V över 2013 … y2+2, högst
  till 2025; flytt = s_D × (2013 − y1) + s_V × (y2 − 2013). Skarvens förskjutning
  ingår inte.
- Minst fyra år krävs i varje del, annars utgår enheten.

**Jämförs med** ändpunkter (nuvarande) och treårsmedel, på samma enheter.

**Mått.** Median |e − t| och rätt riktning, nivå 1 och 2.

**Hypotes.** H2: trendskattningen minskar median |e − t| på nivå 2 med minst 20 %
mot ändpunkter i minst fyra av de fem perioderna (treårsmedel gav 6–15 %).

## 5. A3 – Detektionsmodell

**Fråga.** Går det att ange per enhet hur stor en förflyttning måste vara?

**Brusmått** σ per enhet ur utgåvan: `noise_dmsp_km` för perioder inom DMSP,
√(`noise_dmsp_km`² + `noise_viirs_km`²) över skarven.

**Steg.**
1. Spearmans ρ mellan σ och |e − t|, nivå 2, per facit och period.
2. Logistisk regression: rätt riktning ~ log(|t| / σ). Jämförelsemodell:
   rätt riktning ~ log |t|.
3. Anpassning på USA counties 2000–2020. Prövning på Mexiko 2010–2020 och USA
   2010–2020: förutsagd mot observerad andel i tio lika stora klasser, Brier-värde
   mot konstant andel, och AUC för båda modellerna.
4. Minsta upptäckbara förflyttning per enhet: det |t| där förutsagd sannolikhet
   är 0,8.

**Hypoteser.** H3a: ρ > 0 i alla facit. H3b: modellen med σ har högre AUC än
modellen med bara |t| i prövningen på Mexiko. H3c: förutsagd och observerad andel
skiljer sig med högst 0,05 i minst åtta av tio klasser i Mexiko.

## 6. B – Jämförelsemått

**Fråga.** Hur står sig ljuscentroiden mot det en användare kan få i dag?

**Källor** (tillgång, version och epoker kontrolleras innan något hämtas, och förs
in i `docs/DATA.md`): GHS-POP, WorldPop och bebyggd yta (GHS-BUILT-S), 30″ eller
närmaste upplösning. År: folkräkningsåret eller närmaste epok.

**Metod.** Varje källa läggs in som årsprodukt i `wpc/lights.py` och körs med
samma centroidkod, nivå 1–2, bara för facitländerna. Geometrisk mittpunkt läggs
till som fjärde jämförelse.

**Mått.** Samma läges- och förändringsmått som för ljuset.

**Förbehåll som ska stå i texten.** GHS-POP och WorldPop bygger på folkräkningar.
De är inte oberoende av facit, och förändring mellan epoker är delvis interpolerad.
Jämförelsen visar vad en användare kan få, inte ett oberoende test.

## 7. C – Kanada, förändring helt inom VIIRS-tiden (bekräftande)

**Villkor.** Inga kanadensiska facitdata hämtas eller öppnas innan skripten för
A1–A3 och B är committade. Commit-hashen antecknas under *Avvikelser*.

**Data.** Statistics Canada, Geographic Attribute File 2016 och 2021 (befolkning
per dissemination block). Samma geografiska nivå och samma typ av punkt används
båda åren; finns blockens punkter inte för 2016 används dissemination areas båda
åren. Läsaren läggs i `tools/validation_centroids.py` (`--dataset statcan`).

**Serie.** Utgåvan v2.0.0 som den är: VIIRS 2016 och 2021, GADM 4.1 CAN.

**Huvudmått.** Andel rätt riktning på nivå 2, 2016–2021, med blockbootstrap
(block = provins) och permutationstest.

**Tolkning, bestämd i förväg.**

| Utfall | Slutsats |
| --- | --- |
| Punktskattning ≥ 0,65 och intervallets nedre gräns > 0,50 | VIIRS mäter riktningen bättre än DMSP-tidens och de skarvade perioderna (0,51–0,66). |
| Punktskattning 0,50–0,65 och intervallbredd ≤ 0,15 | Resultatet från DMSP-tiden gäller också VIIRS. |
| Intervallbredd > 0,15 eller färre än 150 enheter | Inkonklusivt. Ett andra facit inom VIIRS-tiden läggs till enligt samma plan, valt innan dess data öppnas. |

**Bimått.** Lägesfel 2021 (median, nivå 1–2); median |e − t| mot median |t|;
detektionsmodellens förutsagda andel mot den observerade (A3); hybridskattningen
på provinsnivå, redovisad beskrivande (10–13 enheter); jämförelsemåtten (B).

## 8. D – Flyttmönster och flöden (senare moment, egen plan)

**Status.** Påbörjas först när A–C är klara. Kan bli en egen artikel. Momentet får
en egen plan (`docs/ANALYSIS_PLAN_D.md`) som låses på samma sätt innan något körs.
Här anges bara ramarna.

**Syfte.** Beskriva hur flyttmönster har förändrats 1992–2025, relatera det till
forskningslitteraturen och visa flöden och deras förändring i kartor och figurer.

**Beroende av A–C.** Resultaten hittills säger att en enhets centroidförflyttning
över 10–20 år för det mesta inte går att skilja från fel på nivå 1–2. Vad som går
att säga om förändring avgörs av detektionsmodellen (A3), hybridskattningen (A1)
och Kanada-testet (C). Momentet får inte bygga på förflyttningar som de testerna
visar är omätbara.

**Två möjliga upplägg, att välja mellan innan planen skrivs.**
1. *Centroidernas banor.* Förflyttningsvektorer per enhet och delperiod, visade som
   vektorfält och banor. Bara för enheter och perioder där förflyttningen
   överstiger enhetens detektionsgräns (A3), eller med hybridskattningen där
   befolkningstal finns. Måste ange om det som följs är befolkningen eller ljuset.
2. *Flöden mellan enheter.* Flyttningar mellan enheter ur yttre data
   (flyttmatriser ur register och folkräkningar), med centroiderna som start- och
   målpunkter för avstånd, som i Niedomysl m.fl. (2017). Använder centroiderna för
   läge ett givet år, vilket är den användning som valideringen stöder.

**Att bestämma i den egna planen.** Vad "flyttmönster" och "flöde" betyder
operativt; länder, nivåer och delperioder; vilka flyttdata som finns och med
vilken licens; vilka frågor ur litteraturen som ska besvaras (efter en egen
litteraturgenomgång); vilka figurer som är huvudresultat; och hur osäkerheten
visas i dem.

## 9. Ingår inte

Annan vikt än linjär radians, vad ljuset följer i stället för boendet
(arbetsplatser, bebyggelse), jämnt fält för samregistreringen, skalanalys av ljus
per invånare och två överlappsår. Se granskningen; de hör till en senare version
eller en egen artikel.

## 10. Avvikelser

| Datum | Avsnitt | Ändring | Skäl |
| --- | --- | --- | --- |
| 2026-11-21 | §2 | scipy lagt till i `requirements.txt` och installerat i miljön. | De frysta måtten i `tools/change_detection.py` använder pandas Spearman, som kräver scipy under pandas ≥ 3. Ingen ändring av koden eller måtten. |
