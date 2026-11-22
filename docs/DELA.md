# Dela materialet med läsare och fylla på

Läget den 22 november 2026. Läsarna får läsexemplaret (PDF). Rutan på första sidan länkar till
datamängden på Zenodo och till repot på GitHub. Läsarnas egen vägledning är
[`READERS.md`](../READERS.md) i repots rot.

## Repot är offentligt

`github.com/ohall2/world-population-centroids` är ett offentligt repo (skapat 2017 för v1). Allt som
pushas dit kan läsas av vem som helst. Läsarna behöver därför ingen inbjudan. Det betyder också:

- **Manuset ligger inte i repot.** `docs/manuscript/` (manus, bilaga, bilagetabeller, talfil,
  läsexemplar) och `docs/review/` (granskningar, talrevision) finns kvar på datorn men står i
  `.gitignore`. Den 22 november 2026 rensades de ur historiken på `v2` och `utokade-kontroller`;
  manuset hade då legat öppet sedan den 18 november. Säkerhetskopia från rensningen:
  `manus_sakerhetskopia_2026-11-22/`.
- **De har ingen versionshistorik längre.** Ta egna kopior, eller lägg dem i ett privat repo.
- **Lägg aldrig Zenodo-länken med nyckel i repot.** Läsexemplaret skickas bara till läsarna.
- Kod, resultattabeller, figurer, facit och analysplan ligger öppet. Ska också de vara stängda tills
  artikeln är inskickad: gör repot privat (*Settings → Danger Zone → Change visibility*; då slutar
  även v1-filerna att synas för utomstående) eller flytta arbetsgrenarna till ett privat repo.

## Vad som ligger var

| Del | Plats | Läge |
| --- | --- | --- |
| Kod, tester, plan, resultattabeller, figurer | GitHub, grenen `utokade-kontroller` | Uppladdat |
| Manus, bilaga, granskningar, läsexemplar | Bara på datorn: `docs/manuscript/`, `docs/review/` | Utanför git |
| Facit per enhet och jämförelsekällornas centroider (11 MB) | `replication/` | Uppladdat |
| Centroidserien v2.0.0 | Zenodo, utkast 24336701 | README och metadata rättade lokalt; laddas upp med steg 2 |
| Rådata, cirka 40 GB | Offentliga källor | Delas inte; källor och kommandon i `docs/DATA.md` |
| GADM-gränser | geodata.ucdavis.edu | Får inte spridas vidare; läsarna hämtar dem själva |

## Det här gör du

1. **Hämta ändringarna till datorn.**

       cd ~/wpc && git pull

2. **Ladda upp den rättade beskrivningen till Zenodo**, om det inte redan är gjort. Skriptet jämför
   kontrollsummor och laddar bara upp `README.md` och `metadata.json`; datafilerna rörs inte och
   utkastet publiceras inte.

       source ~/wpc-venv/bin/activate && python tools/zenodo_upload.py

3. **Kontrollera delningslänken på Zenodo.** Öppna utkastet, välj *Share* och fliken för länkar.
   Läsare ska ha *Can preview drafts*. En länk med *Can edit* låter den som har länken ändra
   utkastet. Länkar kan tas bort och ges ett slutdatum i samma dialog.
   Zenodos anvisning: https://help.zenodo.org/docs/share/link-sharing

4. **Skicka** läsexemplaret med länkrutan. Länken till repot:
   https://github.com/ohall2/world-population-centroids/tree/utokade-kontroller

Förslag till följebrev:

> Here is the reading copy of the manuscript. The box on the first page links to everything else:
> the data set on Zenodo (a private link to an unpublished draft) and the repository with scripts,
> result tables and benchmarks. Start with the guide for readers (READERS.md) in the repository; it
> says where things are and what can be rerun.
>
> The raw data are about 40 GB from public sources and are not redistributed; docs/DATA.md lists the
> sources and the download commands. Please do not circulate the draft or the Zenodo link.

## Lägga till senare

- **Nya skript, resultat eller text.** Committa och pusha till `utokade-kontroller`. Läsarna ser
  ändringen direkt; länkarna i PDF:en pekar på grenen, inte på en viss version.
- **Nytt läsexemplar** efter ändringar i manuset. Kräver pandoc och xelatex:

      python tools/build_review_pdf.py --data-link "<Zenodo-länken>"

  Filen hamnar i `docs/manuscript/manuscript_review.pdf`, som ligger utanför git.
- **Ändrad beskrivning av datamängden.** Texten skrivs av `tools/write_release_docs.py`, som räknar
  spannen ur tabellerna. Kör det och sedan `tools/zenodo_upload.py`. Redigera inte README för hand.
- **Nytt facit eller ny jämförelsekälla.** Kopiera mappen från `data/` till samma plats under
  `replication/data/`, lägg till en rad i `replication/README.md` och committa.

## Kvar att göra innan inskick

- Manuset är inte anonymiserat (GitHub-namn, DOI, versionsnummer), och koden ligger öppet under
  ditt namn. Tidskriften granskar anonymt.
- Villkoren för att sprida tabellerna i `replication/` är inte genomgångna källa för källa.
- En publicerad Zenodo-post kan inte tas bort. Publicera först när artikeln är antagen eller när
  datamängden ska ut oberoende av den.
- `docs/METHODS.md` avsnitt 9 är inte uppdaterat efter de utökade kontrollerna.
