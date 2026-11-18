# Review of the manuscript "Can Nighttime Lights Track Where People Move?"

Reviewer: independent review, 2026-11-18. Object: `docs/manuscript/manuscript.md` (branch v2). I checked the numbers against `docs/validation/`, `docs/figures/data/`, `docs/change_detection/`, `docs/check_identity/`, `docs/figures/data/fig6_seam_jumps.csv`, `config/`, `wpc/lights.py`, `tools/change_detection.py` and `tools/figures/*.py`. I checked the literature against `docs/LITERATURE.md` and ran web checks on five references. No existing file was changed. Rounding is half-up from the CSVs/JSONs. The release data (`data/release/v2.0.0`) is not in the checkout, so splice counts could not be recomputed.

## Summary verdict

**Major revision before submission.** Nearly all the numbers can be traced to the repository and are correctly rounded. Table 2 and the level-2 location results reproduce exactly. The main empirical findings are sound:

- light-weighted centroids lie about 1 km from census centroids at level 2;
- national centroids are displaced by tens to over a hundred kilometres;
- measured municipal shifts over 10–20 years agree only weakly with the true shifts.

Five problems remain:

1. **Overstated claims.**
   - "Cannot track how they move at regional scales" and "direction becomes reliable above 2–5 km" generalise from one benchmark and period (US 2000–2020).
   - At level 1, Mexico's 32 states do detect change (72 %, interval 56–88; ρ = 0.42, interval 0.09–0.73).
   - In Mexico, direction is never reliable at any size class.
2. **The national-bias explanation goes beyond the evidence.**
   - The United States is not decomposed, yet the abstract explains its 65–88 km error.
   - The Bb/Bw split depends on the zoning as well as the grid.
   - The size of the bias depends on the radiometric series (Mexico: 35–42 km forward-spliced, 95–119 km Chiovelli, 133–134 km release). The text does not say this.
   - The explanation of the Nordic mask is wrong for Norway.
   - The Discussion's statement that a correction "would require knowing light per resident by region, which is the information the method is meant to replace" is incorrect. Regional population totals are enough to reweight regional light centroids, and they are widely available.
3. **Five factual and numerical errors:**
   - "65–135 km" should be 65–134;
   - "the same" (61 vs 58 %);
   - 35 composites (34 are used);
   - splice counts taken from an abandoned forward splice;
   - country counts in Table 1.
4. **The seam metric in Fig. 6 does not measure the join for this series.** Its join is 2012→13, not 2013→14.
5. **Undisclosed design dependence.** The splice direction and the 65° N mask were chosen on the same US/European benchmarks that are then used for evaluation. Several robustness checks a referee will ask for are missing: naive baselines, multi-year averaging, a permutation null, US states without Puerto Rico, and zoning sensitivity.

The manuscript is not yet anonymised: there is a byline, internal notes, and "version 2.0.0" of "the first version of the dataset (Hall et al. 2019)".

## Numerical errors

| Location in text | Says | Should be | Source |
| --- | --- | --- | --- |
| Abstract; Recommendation 2 | national displacement "65–135 km" (USA, Mexico, Brazil) | **65–134 km**. USA 65.3 / 80.4 / 87.9; Mexico 134.4 / 133.0; Brazil 94.4. | `figures/data/unit_errors.csv.gz` (level 0) |
| Results, change, para. 3 | US 2000–2010 share 61 %, "the same as for the decade that spans the join" | **61 % (58–64) against 58 % (56–60)**. Not the same, and the intervals barely overlap. Write "61 against 58 percent". | `change_detection/census_change_detection.csv` |
| Data and Methods, Lights | "the 35 annual stable-lights composites" | **34 are used**. F15 2008 has no calibration coefficients (`config/calibration_v1.csv`: 34 rows; `wpc/lights.py` uses only satellites with coefficients). Offsets exist for 35, which may be the source of the confusion. | `config/calibration_v1.csv`, `config/dmsp_offsets.csv`, `wpc/lights.py` l.100–107 |
| Data and Methods, Splicing | spliced units "3,472 of 3,662 regions and 43,441 of 47,217 municipalities" | These are the counts for the abandoned **forward** splice (HARMONIZATION §8). For the release, §11.4 gives `splice_ok = 0` for 13 / 166 / 3,198 units, i.e. **3,496 regions and 44,019 municipalities**. Verify against the release file. | HARMONIZATION.md §8 l.329–330 vs §11.4 l.692 |
| Data and Methods, Splicing | DMSP–VIIRS distance in 2013: "12.5 km / 3.6 km / 1.3 km" | These come from raw, non-co-registered DMSP and possibly mismatched flare masks (HARMONIZATION §1, §5 "Fel i den första jämförelsen"; METHODS harmonisation table). They do not describe the offsets applied in the release (co-registered DMSP against raw VIIRS, 2-km mask). Recompute, or report the distribution of applied offsets. | METHODS.md l.154–156, 194; HARMONIZATION §5 l.136–140 |
| Table 1; Introduction ("26–33 European countries") | GEOSTAT coverage "26 (2011) and 33 (2021) countries" | GADM level-0 units with results: **2011: 25 (levels 0, 2), 26 (level 1, incl. Malta); 2021: 31 / 32 / 29 (levels 0 / 1 / 2)**. These include Åland, Kosovo, Liechtenstein and Akrotiri & Dhekelia, which are not countries. "33" is the pre-correction count from §10, before the 90 % rule in §15.2. | `unit_errors.csv.gz` (geostat); HARMONIZATION §10, §15.2 |
| Results, National Bias, para. 1 | Mexico "131 km north-northwest" | Bearing of D is **325° = northwest**. The 131 km is from the identity computation (raw VIIRS, Aboufadel–Austin, grid). The release error is 133.0 km (Barmore), cited two paragraphs earlier as 133–134. Brazil likewise: 95 here, 94 earlier. Say once that the decomposition uses Aboufadel–Austin on a grid and therefore gives 131 / 95 km. | `check_identity/id_MEX_0.1666667.json` (D = −75.41 E, 107.18 N); `validation/inegi_position.csv` |
| Results, National Bias, para. 2 | masking "reduces the Nordic errors … but does not remove them, because the contrast in light per resident persists south of the mask" | Holds for Sweden only (masked Bb 32 km, still northwards). In **Norway the masked error reverses**: D = (−37 E, −49 N), southwards, because population north of 65° N keeps its weight while its light is removed. HARMONIZATION §11.1 already corrected exactly this explanation, and §14 point 3 confirms it. In Finland the masked error is mainly within-region (Bw 35 of 42 km). | `check_identity/id_NOR_65.json`, `id_FIN_65.json`; HARMONIZATION §11.1 l.629–633, §14 l.938–941 |
| Results, Location, para. 1 | error "shrinks by two orders of magnitude from countries to municipalities" | Ratios: Europe 16× (11.4 → 0.73), USA 60× (87.9 → 1.46), Brazil 85×, Mexico 129×. Write "one to two orders of magnitude". | `unit_errors.csv.gz` |
| Results, Sensor Join; Fig. 6 caption | "In PCNL and in our series … the ratio is 0.8 and 0.7, so the join is no larger than a normal year by construction" | Numerically correct (0.85, 0.70). For **this** series, however, the join is **2012→13** (METHODS §9: DMSP 1992–2012, raw VIIRS from 2013), so 2013→14 is a VIIRS–VIIRS step. The actual join step is 1.08 km, i.e. 0.92× the median of the three DMSP-era steps. That conclusion still holds, but the metric and caption must use each series' own join year. For PCNL (regressed VNL in 2013, raw VNL from 2014) 2013→14 is appropriate. | `figures/data/fig6_seam_jumps.csv`; METHODS.md §9; LITERATURE §5 |
| Discussion, "A Good Snapshot", para. 1 | "unless those errors are strongly correlated between years, [the change error] is of the same order as the location error … they are not correlated enough" | The data imply the **opposite framing**. Median change error is 0.49–0.79 × the median location error, which (under isotropic normal errors) implies a year-to-year error correlation of about **0.7–0.9** (US 0.75–0.88, Mexico 0.68, Europe 0.71). The errors *are* strongly correlated. The residual still exceeds true shifts because those shifts are only 0.13–0.54 km. | `validation/*_position.csv`, `change_detection/*_change_detection.csv` (err_med) |
| Discussion, "Why National Centroids Fail" | "The Barmore and Aboufadel–Austin centroids differ by less than 0.01 km at the first and second levels" | Correct only **in median** (as the Methods says). Differences reach about 8 km for the USA and 25 km for Russia at level 0 (METHODS §6). The level-0 validation itself is the better evidence: Barmore and Aboufadel–Austin errors differ by < 0.2 km for Mexico and Brazil. | METHODS.md §6; `validation/*_position.csv` (aa_median) |
| Abstract; Results, Fig. 5 para. | "Direction becomes reliable only for shifts above 2–5 km" | Holds only for US 2000–2020 (2–5 km: 84 %, n = 170). US 2000–2010: 77 % (n = 78). US 2010–2020: 70 % (n = 33). **Mexico 2–5 km: 58 % (n = 62)**, below its 1–2 km class (67 %). Qualify. | `change_detection/*_bins.csv` |

**Checked and correct.** All other numbers check out:

- **Location** (`unit_errors.csv.gz`, `*_position.csv`):
  - level 2: US 1.46–1.52; Europe 0.73–0.95; Mexico 1.03–1.12; Brazil 1.11; p90 2.41–4.30;
  - level 1: Puerto Rico 0.33–0.36; Europe 2.47–2.97; Mexico 7.4–8.9; US states 10.5–13.1; Brazil 16.5;
  - national: European median 11.4–13.3.
- **Table 2**: all cells, including every bootstrap bound.
- **Change-detection text**:
  - shares 51–66 % right direction; rank correlations 0.02–0.34; error exceeds true shift in 72–92 %; slopes 0.26–0.49;
  - bins 54 / 69 / 84 %; 8.4 % of counties moved more than 2 km; Europe 50–52 %;
  - Chiovelli: right direction 0.45–0.62, rank correlation −0.12 to 0.25.
- **Decomposition**:
  - 127/131, 92/95, 95/98, 146/149; Finland 38 / 34 / 71; A and C < 2 km; A = 115 at 1′;
  - masked errors 33 / 61 / 42;
  - light per resident (r) values at 10′: 1.72 / 1.76, 0.54 / 0.62, 1.69 / 0.60, 3.07 / 3.98 / 0.29.
- **Seam**: ratios 3.3 / 1.6 / 1.5; 2016–17 step 1.42× in this series; n = 2,975.
- **Methods**: flare-mask percentages (METHODS, VIIRS 2013); co-registration up to 1.44 px; abstract word count 183.

## Literature and references

**Citations versus reference list.** Every in-text citation has an entry, and every entry is cited. Data citations are missing, as the draft itself notes:

- DMSP-OLS v4;
- VNL v2.1/v2.2;
- MYC25Q4 (Zenodo 10.5281/zenodo.23091822);
- Chiovelli data (Zenodo 10.5281/zenodo.20558766);
- Li et al. 2020 data (figshare);
- PCNL (Zenodo 7612389).

**Claims against LITERATURE.md.** Most claims match what is recorded: Asher, Goldblatt, Bluhm & McCord, Bennett & Smith, Mirza, Gaughan, Nieves, Archila Bustos, Hall 2019, Niedomysl 2017 (full text), PCNL, Chiovelli (corrected characterisation), Zhang & Gibson, Quah, Grether & Mathys. Issues:

1. **"That use rests on an assumption that has not been tested"** (Introduction). LITERATURE §2 and §12 record this as a *negative search result* and say it should be phrased "to our knowledge". Add the hedge.
2. **L. Li et al. 2026.** It is verified only through a third-party summary (synapsesocial). The publisher page was not read. The 53.25 % figure matches that summary. Before relying on it, which the Discussion does ("the similarity is close"), check the publisher page or full text. The comparison is also loose: direction of a *count* change has a chance baseline of 50 % only if increases and decreases are equally common.
3. **Plane and Rogerson 2015.** The abstract confirms "five … mean and median center points" and the decomposition into births, deaths, internal migration and immigration ([IDEAS](https://ideas.repec.org/a/taf/raagxx/v105y2015i5p968-986.html)). The claim is fine, but these details are not recorded in LITERATURE.md. "The U.S. Census Bureau has published the national mean centre … every decade since 1790" is imprecise: the early centres were computed retrospectively. Write "has computed … for every census since 1790".
4. **Barmore 1991, Aboufadel and Austin 2006, Elvidge et al. 2009 and 2021, Zhizhin et al. 2026** are not in LITERATURE.md. Verify them and record them there.
5. **"Our decomposition and Mirza et al. (2021) indicate [that light per resident varies with income]"** (Discussion). The decomposition uses no income data. HARMONIZATION §15.4 notes that rich metropolitan regions (Mexico City, São Paulo, Rio) have *low* r. Limit the decomposition's support to "varies by region".
6. **"An independent harmonized series"** (Introduction) for Chiovelli et al. Its DMSP part builds on the same composites as Li et al. 2020, and it is tested against the same benchmarks. "An alternative harmonized series" is accurate.
7. **"Studies that read 10–30-year movements of light-weighted centres … as shifts of population or activity"** (Discussion). The tests are against population only. Studies that track the centre of *light* as a measure of urban expansion are not invalidated by a census-population benchmark. Restrict the claim to population, or argue the link explicitly.
8. **Missing literature a referee will expect:**
   - Sutton et al. 2001, *IJRS* 22(16) (the light–population relation varies with income and density; a predecessor of the r mechanism);
   - Dong et al. 2021, *Remote Sensing* 13:3142 (shrinking cities rarely show declining light);
   - Gibson 2021, *OBES* 83(3); Martinez et al. 2023 (both in LITERATURE.md);
   - basic MAUP literature for the zoning dependence of Bb/Bw (e.g. Openshaw 1984).

**Reference entries:**

- **Wang et al. 2026.** The publisher lists **8 authors**: Wang J., Hou P., Liu Y., Zhai J., Chen Y., Liu Y., Xie L., **Ren X.** ([journal page](https://www.ygxb.ac.cn/en/article/doi/10.11834/jrs.20254567/)). Chicago author-date lists up to ten authors in full, so replace "et al." with Ren, X. Year 2026, vol. 30(1):61–78 and the DOI are correct. LITERATURE.md gives the year as 2025; make it consistent. The publisher abstract confirms centre of gravity and standard deviational ellipse but does not name the sensors.
- **Zhizhin et al. 2026.** Add **doi:10.3390/rs18020314** ([MDPI](https://www.mdpi.com/2072-4292/18/2/314)). The other details are correct.
- **X. Li et al. 2026.** The full names are Xiaoping Li, Yanhui Lei, **Junnan Shen, Rui Xin, Kunpeng Qu**. The quoted sentence is verified in the Fig. 4 caption ([Springer](https://link.springer.com/article/10.1007/s10661-026-16129-9)). The earlier review could not verify this reference; it is now verified.
- **L. Li et al. 2026.** With 13 authors, Chicago 17 author-date lists the first seven followed by "et al." Volume and article number are still missing.
- **Chiovelli et al. 2026.** Volume and pages are pending. Price 1972 has no DOI, which is acceptable.
- **Note above the reference list.** Remove it. Its statement that "every entry … has been checked against the publisher" is not true for L. Li et al. 2026 or Wang et al.

## Argument and evidence

**Reconstructed argument.** Light-weighted centroids are good cross-sectional locators within municipalities. They are biased at national scale because light per resident covaries with position. They cannot detect decadal shifts at levels 1–2 because the noise in a measured shift is larger than the true shift. Hence they should be used for location, not for change.

**Strengths.**

- Four census systems, with Brazil held out.
- Per-unit metrics with block-bootstrap intervals.
- A comparison series, and a DMSP-only period.
- Honest limitations: no change test entirely within VIIRS, a weak European benchmark, the 2016–17 step.
- The decomposition is explicitly called an identity, not a causal estimate, and its grid dependence is shown. This is presented honestly.

**Major concerns.**

1. **Generalisation of the change result.**
   - *What the evidence shows:* three benchmarks, all in upper-middle or high-income countries with slow redistribution, and only two census years each. The headline is stated as a general property ("cannot track how they move at regional scales"; Conclusion: "cannot detect how the population centre of a region or municipality moves").
   - *The data support a conditional claim:* detection depends on the size of the shift relative to the noise. Where shifts are about 2 km (Mexican states), direction is right in 72 %. Where shifts are 0.1–0.5 km, it is near chance.
   - *Detection threshold:* rests on US 2000–2020 alone. Mexico never exceeds 67 % in any class.
   - *Required:* restate the claim as "shifts of the size observed here (median 0.13–0.54 km) are not detectable". Report the threshold per benchmark. Bring the Mexico level-1 result into the abstract and Conclusion rather than leaving it as a qualifier.
2. **Missing baselines and robustness checks a referee will demand.**
   - (a) *Naive location baseline.* How much better is the light centroid than the polygon's geometric centroid, or a WorldPop or GHS-POP centroid? Without this, "locate populations well" has no reference point.
   - (b) *Naive change baseline.* "True shift > error" in only 5–28 % means that assuming *no movement* beats the lights in 72–95 % of units. Say so explicitly. It is the clearest statement of the result.
   - (c) *Multi-year averaging.* The Discussion attributes the residual to year-to-year noise. Averaging 3–5 annual composites at each end is the obvious remedy and has not been tested. For the US 2000–2020 it can be done with existing data. Without it, the conclusion applies to single-year composites only.
   - (d) *Null distribution.* The 50 % chance level assumes that the direction of e is independent of and uniform relative to t. Chiovelli in Europe gives 45 % (interval 43–46 %) and European regions give a *significantly negative* rank correlation (−0.16 to −0.01). That is systematic below-chance behaviour, pointing to benchmark or sensor artefacts. Use a within-block permutation null, and discuss the negative results instead of reading them as chance.
   - (e) *US states.* Report level 1 for the 50 states and DC, excluding Puerto Rico. The caption's reason for omitting it applies to the pooled sample, not to the states.
   - (f) *Interval method.* The Fig. 5 Wilson intervals ignore spatial clustering, whereas Table 2 uses a block bootstrap. Harmonise them, or note the difference.
3. **National bias: explanation, scope and the correction claim.**
   - (a) The abstract explains the US, Mexico and Brazil errors "because light per resident differs systematically between regions". The **USA is not decomposed**, and Fig. 1 does not show it. Decompose the USA (block-group population points exist) or restrict the claim.
   - (b) Bb/Bw depends on the zoning as well as the grid. With k-means merges of the same data, Bb goes from 0 to 95 km in Sweden as the number of regions goes from 1 to 21 (earlier review, B4). "Dominated by between-region contrasts" is a statement about 21–32 GADM regions. Show a zoning sensitivity, or rephrase it as "at the scale of first-level regions".
   - (c) The magnitude and even the direction depend on the light product. For Mexico: 133–134 km in the release, 35–42 km forward-spliced (DMSP-like), 95–119 km in Chiovelli. For Brazil: 94 km pointing SW in the release, 73 km pointing NE with DMSP-like VIIRS (`*_position.csv`; HARMONIZATION §10, §15; LITERATURE §11). The manuscript should say that the 65–134 km bias is a property of raw VIIRS radiance weights.
   - (d) **The Discussion's correction claim is wrong.** The bias can be largely removed by weighting regional light centroids by regional population totals, which are available almost everywhere even without fine-scale data. This would remove most of Bb by construction. Test it with existing data, or at least correct the sentence. It is also a useful practical recommendation.
   - (e) Regional r is computed over populated cells only, so it is grid-dependent. Chihuahua is 1.72 at 10′ but 0.92 at 1′, because INEGI reports large cities as one point. The 1′ "A" term in Mexico is largely an artefact of that benchmark, not light in genuinely unpopulated areas. Define r, and state this.
4. **Undisclosed dependence on the evaluation data.**
   - The backward splice (anchored on VIIRS) was chosen because it gave the best location errors against the US and European benchmarks (HARMONIZATION §10, "Förbehåll").
   - The 65° N mask was motivated by the Nordic benchmark errors (§11).
   - Only Brazil is out of sample, and Brazil tests no splice because it has only 2022.
   - The location stability across census years ("the 2000 and 2010 values come from DMSP while the 2020 value comes from VIIRS") is partly built in, because each unit's DMSP years are shifted onto its VIIRS 2013 position.
   - *Required:* state all of this in Methods and Results.
5. **Recommendation 3 is not backed by census results.** PCNL was never tested against census data. Fig. 6 shows only continuity at the join, which for PCNL and for this series holds by construction. The Li DN > 7 2016→17 step (2.43 km) is larger than its join step (2.28 km), so the join is not uniquely anomalous. Limit the recommendation to what was tested: the per-unit splice against Chiovelli, for location and change.
6. **Discussion claims without Results.**
   - "Continuous in global totals" (Recommendations) is not shown in Results. It comes from HARMONIZATION §12.
   - "Gridded products may inherit the same regional bias" is correctly flagged as untested.
   - "Sparsely populated northern and agricultural regions … pull the centre … away from dense metropolitan regions" generalises from five countries.

**Minor concerns.**

- Niedomysl et al. (2017) used local labour-market regions (median error about 4 km), not municipalities. "Where the error is about 1 km" does not describe that use.
- "Little more often than chance" undersells US 2000–2020 (66 %, 64–69). Use calibrated wording.
- An alternative source of non-demographic change in light per resident over 2010–2020 should be mentioned: the LED transition, to which the VIIRS DNB is partly insensitive.
- The Results do not report Chiovelli location results, although Methods says location tests were run on that series.

## Consistency

- **Mexico and Brazil national errors** are 133–134 / 94 in Results para. 1 and 131 / 95 in para. 2 and Figs. 1 and 3. Explain the difference (decomposition uses Aboufadel–Austin on a grid with raw VIIRS).
- **Centroid formula.** Methods says the Barmore centroid is reported. The decomposition (`check_identity.py`) and Fig. 6 (`fig6_seam.py`, `aa_lon` / `aa_lat`) use Aboufadel–Austin, and Fig. 6's comparison series use simplified polygons with centre-in-polygon rasterisation. Disclose this.
- **Figure order.** Figure 2 is cited before Figure 1 (Results, Location). Renumber, or cite Fig. 1 first.
- **Fig. 6.** The caption and text call 2013–14 "the join" for all series. For this series the join is 2012–13 (see table above).
- **Fig. 4 caption** does not give the periods (USA 2000–2020, Mexico 2010–2020, Europe 2011–2021). Only the panel titles do.
- **Fig. 2** labels Mexico level 2 as "1 km" (format `.2g` of 1.03). This is cosmetic but inconsistent with the "1.0–1.1" in the text.
- **Methods omissions:**
  - VNL v2.1 (2013–2021) and v2.2 (2022–2025) are mixed; Brazil 2022 uses v2.2.
  - The flare mask also uses EOG annual catalogues 2012–2024, including LNG plants and refineries.
  - Midland–Odessa loses 65 % of its light at 2 km.
  - The Chiovelli variant (bloomtopcode_fix) is not named.
  - INEGI records large cities as one point per municipality. This limits level-2 true shifts in urban municipalities and drives the 1′ A term.
- **"Harmonized series"** (abstract) versus "spliced". The release series is a per-unit splice of co-registered DMSP and raw VIIRS, not a radiometric harmonisation. Use one term consistently.
- **Repository inconsistency (not in the manuscript).** HARMONIZATION §14 gives Norrbotten 2.33 and Finnmark 3.00 "(rutnät 1′)" under a 10′ heading. The JSON 10′ values (3.07, 3.98) match the manuscript. Fix the documentation.

## Journal requirements

- **Length.** About 6,600 words including references (body about 5,500), within the 11,000 limit. The abstract is 183 words, within 200. Five key words.
- **Chicago author-date (Annals variant)** is followed in form: "Journal Vol (Issue): pages. doi:". Fix the reference entries listed above and add the data citations. Remove the bracketed placeholders and the internal note above the reference list.
- **Double-anonymous review. Remove or rephrase:**
  - the byline line "Nov 18, 2026 · @O";
  - the "Target: …" note;
  - "the second-order coefficients used for the first version of the dataset (Hall et al. 2019)" and "version 2.0.0" (Data availability; Fig. 6 legend "This series (v2.0.0)"). Together these present the authors as producers of the Hall et al. dataset. Write "the coefficients of Hall et al. (2019)" and "the series";
  - the repository URL (it contains the author's GitHub handle). Use an anonymised mirror. The code and docs comments are in Swedish, which narrows identification further.
  - The Zenodo DOI must be withheld or anonymised for review.
- **Self-citations** (Hall et al. 2019; Niedomysl et al. 2017; Archila Bustos et al. 2015) are written in the third person, which is acceptable. "Limits the use proposed by Hall et al. (2019)" is fine.
- **GADM licence.** It forbids redistribution without permission. Confirm that the deposited series contains no GADM geometry.

## Prioritised list of required changes

1. **Recalibrate the change claim.** Make it conditional on shift size. Report Mexico level 1 as positive evidence of detection. State the threshold per benchmark (US only; not reached in Mexico). Edit the abstract, Results heading and Conclusion accordingly.
2. **Fix the national-bias section:**
   - decompose the USA or drop it from the causal sentence;
   - state the dependence on the light product (35–134 km for Mexico across series);
   - correct the Nordic mask explanation (Norway reverses);
   - replace the "correction requires light per resident" claim and test population-weighted regional light centroids;
   - add a zoning sensitivity for Bb/Bw, or restrict the claim to first-level regions;
   - define r and note its grid and INEGI dependence.
3. **Correct the numerical errors** in the table above: 65–134; "61 against 58"; 34 composites; splice counts; European country counts and composition; NW not NNW; orders of magnitude; the Discussion's error-correlation argument; "in median" for Barmore versus Aboufadel–Austin.
4. **Fix the Fig. 6 metric** to use each series' own join year (2012–13 for this series). Remove "by construction" as evidence of quality. Restrict Recommendation 3 to what was tested.
5. **Disclose** that the splice direction and the mask were chosen on the US and European benchmarks, and that location stability across years is partly built in by the splice.
6. **Add robustness checks:**
   - naive baselines for location (geometric centroid, a gridded population product) and for change (zero shift);
   - multi-year averaging at each end;
   - a permutation null within blocks, and discussion of the significant below-chance European results;
   - US states without Puerto Rico;
   - block-bootstrap or design-effect intervals in Fig. 5.
7. **Anonymise** (byline, notes, version and "first version" wording, repository, DOI). Complete and correct the references (Wang: add Ren X.; Zhizhin DOI; X. Li initials; L. Li verification and et al.; data citations). Add the hedge "to our knowledge". Add Sutton et al. 2001 and Dong et al. 2021.
8. **Minor:**
   - figure order;
   - Fig. 4 caption periods;
   - methods details (VNL v2.1/v2.2, flare catalogues, Chiovelli variant, Aboufadel–Austin use in the decomposition and Fig. 6);
   - Niedomysl et al. (2017) error scale;
   - "harmonized" versus "spliced".
