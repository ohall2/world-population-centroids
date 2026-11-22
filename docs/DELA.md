# Dela materialet med läsare och fylla på

Läget den 22 november 2026. Läsarna får läsexemplaret (PDF). Rutan på första sidan länkar till
datamängden på Zenodo och till repot på GitHub. Läsarnas egen vägledning är
[`READERS.md`](../READERS.md) i repots rot.

## Repot är offentligt

`github.com/ohall2/world-population-centroids` är ett offentligt repo (skapat 2017 för v1). Allt som
pushas dit kan läsas av vem som helst, också grenarna `v2` och `utokade-kontroller` med manus,
granskningar och facittabeller. Läsarna behöver därför ingen inbjudan. Det betyder också:

- Manuset ligger öppet med namn, och tidskriften granskar anonymt.
- Lägg aldrig Zenodo-länken med nyckel i repot. Läsexemplaret i repot byggs utan länk; exemplaret
  med länk skickas bara till läsarna.
- Ska arbetet vara stängt tills artikeln är inskickad finns två vägar: gör repot privat
  (*Settings → Danger Zone → Change visibility*; då slutar också länkar till v1-filerna att fungera
  för utomstående), eller flytta arbetsgrenarna till ett nytt privat repo och ta bort dem här.

## Vad som ligger var

| Del | Plats | Läge |
| --- | --- | --- |
| Kod, tester, plan, granskningar, resultattabeller, manus, figurer | GitHub, grenen `utokade-kontroller` | Uppladdat |
| Läsexemplar utan Zenodo-länk | `docs/manuscript/manuscript_review.pdf` | Uppladdat |
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

      python tools/build_review_pdf.py --data-link "<Zenodo-länken>" --out ~/Desktop/WPC_reading_copy.pdf
      python tools/build_review_pdf.py

  Det första kommandot ger exemplaret till läsarna, med klickbar länk, utanför repot. Det andra
  uppdaterar exemplaret i repot, utan länk.
- **Ändrad beskrivning av datamängden.** Texten skrivs av `tools/write_release_docs.py`, som räknar
  spannen ur tabellerna. Kör det och sedan `tools/zenodo_upload.py`. Redigera inte README för hand.
- **Nytt facit eller ny jämförelsekälla.** Kopiera mappen från `data/` till samma plats under
  `replication/data/`, lägg till en rad i `replication/README.md` och committa.

## Kvar att göra innan inskick

- Manuset är inte anonymiserat (GitHub-namn, DOI, versionsnummer).
- Villkoren för att sprida tabellerna i `replication/` är inte genomgångna källa för källa.
- En publicerad Zenodo-post kan inte tas bort. Publicera först när artikeln är antagen eller när
  datamängden ska ut oberoende av den.
- `docs/METHODS.md` avsnitt 9 är inte uppdaterat efter de utökade kontrollerna.
