# Dela materialet med läsare och fylla på

Läget den 22 november 2026. Läsarna får tre saker: läsexemplaret (PDF), åtkomst till det privata
repot på GitHub och en förhandslänk till datamängden på Zenodo. Läsarnas egen vägledning är
[`READERS.md`](../READERS.md) i repots rot.

## Vad som ligger var

| Del | Plats | Läge |
| --- | --- | --- |
| Kod, tester, plan, granskningar, resultattabeller, manus, figurer | GitHub, grenen `utokade-kontroller` | Uppladdat |
| Läsexemplar med länkruta | `docs/manuscript/manuscript_review.pdf` | Uppladdat |
| Facit per enhet och jämförelsekällornas centroider (11 MB) | `replication/` | Uppladdat |
| Centroidserien v2.0.0 | Zenodo, utkast 24336701 | README och metadata är rättade lokalt men inte uppladdade (steg 2) |
| Rådata, cirka 40 GB | Offentliga källor | Delas inte; källor och kommandon i `docs/DATA.md` |
| GADM-gränser | geodata.ucdavis.edu | Får inte spridas vidare; läsarna hämtar dem själva |

## Det här gör du

1. **Hämta ändringarna till datorn.**

       cd ~/wpc && git pull origin utokade-kontroller && git branch -u origin/utokade-kontroller

2. **Ladda upp den rättade beskrivningen till Zenodo.** Skriptet jämför kontrollsummor och laddar bara
   upp `README.md` och `metadata.json`; datafilerna rörs inte och utkastet publiceras inte.

       source ~/wpc-venv/bin/activate && python tools/zenodo_upload.py

3. **Skapa förhandslänken på Zenodo.** Öppna utkastet (https://zenodo.org/deposit/24336701), välj
   *Share*, fliken för länkar, behörigheten *Can preview drafts*, och skapa länken. Etiketterna kan
   skilja sig något. En länk från förra omgången gäller fortfarande om den inte har gått ut.
   Zenodos anvisning: https://help.zenodo.org/docs/share/link-sharing

4. **Ge läsarna åtkomst till repot.** GitHub: *Settings → Collaborators → Add people*. Rollen *Read*
   räcker. Hoppa över steget om de redan är inbjudna.

5. **Skicka** PDF:en, förhandslänken och länken till repot:
   https://github.com/ohall2/world-population-centroids/tree/utokade-kontroller

Förslag till följebrev:

> Here is the reading copy of the manuscript. The box on the first page links to everything else.
> Start with the guide for readers (READERS.md) in the repository; it says where the scripts, result
> tables and benchmarks are and what can be rerun.
>
> Data set (draft, not public): [paste the Zenodo preview link]
> Repository (private, you have been invited): https://github.com/ohall2/world-population-centroids/tree/utokade-kontroller
>
> The raw data are about 40 GB from public sources and are not redistributed; docs/DATA.md lists the
> sources and the download commands. Please do not circulate.

## Lägga till senare

- **Nya skript, resultat eller text.** Committa och pusha till `utokade-kontroller`. Läsarna ser
  ändringen direkt; länkarna i PDF:en pekar på grenen, inte på en viss version.
- **Nytt läsexemplar** efter ändringar i manuset. Kräver pandoc och xelatex:

      python tools/build_review_pdf.py --data-link "<förhandslänken>"

  Med `--data-link` blir länken klickbar i rutan på första sidan; utan den hänvisar rutan till
  följebrevet. Committa sedan `docs/manuscript/manuscript_review.pdf`.
- **Ändrad beskrivning av datamängden.** Texten skrivs av `tools/write_release_docs.py`, som räknar
  spannen ur tabellerna. Kör det och sedan `tools/zenodo_upload.py`. Redigera inte README för hand.
- **Nytt facit eller ny jämförelsekälla.** Kopiera mappen från `data/` till samma plats under
  `replication/data/`, lägg till en rad i `replication/README.md` och committa.
- **Fler läsare.** Steg 4 för repot; samma förhandslänk fungerar för Zenodo.

## Innan något görs offentligt

- Manuset är inte anonymiserat (GitHub-namn, DOI, versionsnummer), och tidskriften granskar anonymt.
- Villkoren för att sprida tabellerna i `replication/` är inte genomgångna källa för källa.
- En publicerad Zenodo-post kan inte tas bort. Publicera först när artikeln är antagen eller när du
  har bestämt att datamängden ska ut oberoende av den.
- `docs/METHODS.md` avsnitt 9 är inte uppdaterat efter de utökade kontrollerna.
