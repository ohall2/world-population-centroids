# Can nighttime lights track where people move? – manuscript (Annals of the AAG)

Nov 18, 2026 · @O

## Title, abstract and key words

**Can Nighttime Lights Track Where People Move? Census Evidence on Light-Weighted Population Centroids, 1992–2025**

*Target: Annals of the American Association of Geographers, Geographic Methods. Limit 11,000 words including abstract, references, notes, tables and captions; abstract ≤ 200 words; Chicago author-date; double-anonymous review.*

**Abstract.** Light-weighted centroids are used to describe where populations live and how their centres move over decades, often across the 2013 change from the DMSP-OLS to the VIIRS sensor. We test both uses against census population centroids in the United States, Europe, Mexico and Brazil, with a spliced series of centroids for 1992–2025 at three administrative levels. Within counties and municipalities the light-weighted centroid lies a median 0.7–1.5 km from the census centroid. National centroids of the United States, Mexico and Brazil are displaced by 65–134 km; in Mexico, Brazil and three Nordic countries an exact covariance decomposition attributes most of the displacement to differences in light per resident between first-level regions. Change is a different matter. Over 10–20 years, municipal centres of population moved a median 0.13–0.54 km. The measured shift points the right way in 51–66 percent of municipalities, against 50 percent by chance, and in 72–92 percent assuming no movement is closer to the truth. Larger shifts are detected more often: the direction is right in 72 percent of Mexican states, whose centres moved a median 2 km. Nighttime lights locate populations well but cannot track shifts of the size observed in these benchmarks. \[195 words\]

**Key Words:** census validation, centre of population, change detection, nighttime lights, VIIRS

## Introduction

The centre of population is one of geography's oldest summary measures, and its movement is read as a record of redistribution. The U.S. Census Bureau has computed the national mean centre of population for every census since 1790, and its westward and southward drift is a familiar shorthand for two centuries of settlement history (Plane and Rogerson 2015). Centre points are also used to summarize economic activity at the global scale (Quah 2011; Grether and Mathys 2011). Computing them requires population counts at fine resolution and comparable units over time. Such data are missing or irregular in much of the world, and where they exist they are tied to decadal censuses with changing boundaries.

Nighttime lights seem to offer a way around this. Annual composites from the Defense Meteorological Satellite Program Operational Linescan System (DMSP-OLS, 1992–2013) and the Visible Infrared Imaging Radiometer Suite (VIIRS, from 2012) cover the world at about 1 km resolution, and light is strongly associated with where people live. Niedomysl et al. (2017) showed in this journal that light-weighted centres lie close to census centres of population and can be used to estimate migration distances. Hall et al. (2019) built on that result to publish annual light-weighted centroids for the world's administrative units, 1992–2013, and proposed that they be used to analyse changes in population distribution. A growing literature does exactly that with centres of gravity and standard deviational ellipses of light, typically over 20–30 years and across the change of sensor in 2013 (e.g., Dong et al. 2020; Wang et al. 2026).

That use rests on an assumption that, to our knowledge, has not been tested: that the movement of a light-weighted centre reflects the movement of the population. Validation of light-weighted centres has so far been cross-sectional, comparing locations in a single year. A sound location in each year does not imply a sound change between years. If the error of the light-weighted centre is of the same order as the true shift of the population centre over a decade, the measured movement is mostly noise. Evidence from neighbouring problems makes this a real possibility. Nighttime lights predict the level of economic activity far better than its change (Goldblatt, Heilmann, and Vaizman 2020; Asher et al. 2021), gridded population products reproduce the location of census populations but not their change (Zhang and Gibson 2025; L. Li et al. 2026), and population decline in Europe is not matched by declining light (Archila Bustos, Hall, and Andersson 2015).

This article tests both uses of light-weighted centroids — location in a given year and change over time — against census population centroids. We construct a series of centroids for 1992–2025 at three administrative levels for all countries, splicing DMSP and VIIRS at the overlap year, and compare it with census benchmarks for the United States (2000, 2010, 2020), 25–32 European countries and territories (2011, 2021), Mexico (2010, 2020) and Brazil (2022). We ask four questions. How accurately do light-weighted centroids locate population at each level? Why are national centroids biased? Can the shift of a centroid over 10–20 years be detected, and how large must it be? Does the answer depend on how the two sensors are joined?

The answers separate the two uses. Within counties and municipalities the light-weighted centroid lies about 1 km from the census centroid, and at that scale it is a usable proxy for where people live. National centroids are displaced by tens to more than a hundred kilometres; where we decompose the displacement, it is dominated by differences in light per resident between regions. Change is detected only where it is large. The shifts that dominate our benchmarks, a few hundred metres per municipality over a decade, are measured in the right direction little more often than chance, while shifts of several kilometres are detected more often. The finding holds for an alternative harmonized series and within the DMSP era alone. It limits the use proposed by Hall et al. (2019) and the practice that has grown around it, and it extends the level–change gap documented for economic activity and gridded population to the location of population.

## Light as a Weight for Population Location

### Centres of Population

A centre of population is a weighted mean of locations, and its value depends on the weights as much as on the geometry. On a sphere the mean can be defined in more than one way. Barmore (1991) proposed iterating the planar mean in an azimuthal equidistant projection centred on the current estimate; Aboufadel and Austin (2006) averaged three-dimensional unit vectors and projected the mean back to the surface. Plane and Rogerson (2015) compared five mean and median centres for the United States and showed how the centre can be decomposed into the contributions of births, deaths and migration. All these definitions assume that the weights are counts of people. When the weights are light, the centre of light coincides with the centre of population only if light per resident is the same everywhere in the unit.

### Light per Resident

Light per resident is not constant, and it varies with income and settlement form. Mirza et al. (2021) used the spatial distribution of light per person as a measure of income inequality, on the grounds that wealthier people live in larger, better lit dwellings and neighbourhoods. Gaughan et al. (2019) found that distributing emissions by nighttime lights concentrates them differently from distributing them by population. For centroids the consequence is direct. Any region whose residents produce more light than the national average pulls the light-weighted centre towards itself; any region producing less pushes it away. Niedomysl et al. (2017) found that the distance between light-weighted and population centres of Swedish local labour market regions was uncorrelated with the shape and population density of the regions and only very weakly correlated with their area, a correlation they attributed to six very large, sparsely populated regions in northern Sweden. That is what one would expect if the error is governed by light per resident rather than by geometry, and the northern outliers are where light per resident is highest.

### Joining Two Sensors

A long series requires joining DMSP-OLS and VIIRS, which differ in resolution, saturation, detection threshold and geolocation. Published harmonized series convert one sensor to the other's scale and are evaluated on aggregates: total light and lit area over time (X. Li et al. 2020) or pixel-level agreement with the target sensor (Chiovelli et al. 2026). The Pixel-scale Corrected Nighttime Light dataset (PCNL) instead anchors DMSP on VIIRS in the overlap year, adding a 2013 difference raster to every earlier DMSP year (S. Li et al. 2023). How the join behaves for a single administrative unit, and therefore for a centroid, is rarely reported. X. Li et al. (2026) note that the absence of a jump in aggregate totals does not by itself establish temporal stability at the county level.

### Level Versus Change

A proxy can be accurate in levels and uninformative about change, because the variance of change is small relative to the variance of levels while the proxy's error does not shrink in proportion. Studies of economic activity show this repeatedly. Village-level lights in India are significant predictors of population, employment and consumption, but with much lower elasticities in time series than in cross-section (Asher et al. 2021). Neither nighttime lights nor Landsat imagery predict changes in local economic activity in Vietnam well (Goldblatt, Heilmann, and Vaizman 2020), and the ability of lights to capture change depends on local conditions (Bluhm and McCord 2022). For population, four gridded products explain more than 95 percent of the variation in census population across Chinese counties but much less of its change (Zhang and Gibson 2025). Six gridded products assessed against Chinese township censuses get the direction of decadal change right in 53 percent of townships on average (L. Li et al. 2026). In Europe, nighttime lights did not decline where population did (Archila Bustos, Hall, and Andersson 2015). A review of multitemporal applications takes the more optimistic view that time series of lights can track socioeconomic dynamics reliably (Bennett and Smith 2017). None of these studies concerns the location of population, which is what centroid-based studies of redistribution measure.

## Data and Methods

### Lights

The DMSP-OLS series consists of the version 4 stable-lights composites (1992–2013, 30 arc-seconds); the 34 composites with intercalibration coefficients are used, and F15 2008 is not. Digital numbers are intercalibrated with the second-order coefficients of Hall et al. (2019), with satellite F18 in 2010 as reference; pixels with a value of zero remain zero, and years with two satellites are averaged. VIIRS lights are the annual VNL version 2 composites (v2.1 for 2013–2021, v2.2 for 2022–2025), band average\_masked, for 2013–2025 (15 arc-seconds) (Elvidge et al. 2021). Gas flares are removed with one static mask for all years and both sensors: the flare polygons of Elvidge et al. (2009) plus a 2-km buffer around every flare site in the VIIRS Nightfire catalogue for 2012–2025 (Zhizhin et al. 2026) and in the annual flare catalogues of the Earth Observation Group for 2012–2024, which include liquefied natural gas plants and refineries. The 2-km radius removes 83–99 percent of the light in the Bakken and Permian oil fields while removing 6–10 percent of the light in Houston, Edmonton and Oklahoma City, but 65 percent in Midland–Odessa, which lies inside the Permian field.

The DMSP composites are not geolocated consistently with each other or with VIIRS, a known problem (Tuttle et al. 2013; Ash and Mazur 2020). We co-register each composite to VIIRS 2013, aggregated and smoothed to the DMSP grid, by cross-correlating log-transformed light in 5° × 5° tiles and taking the weighted median shift. The estimated shifts reach 1.4 pixels (about 1.3 km), and almost all composites lie north or northwest of VIIRS, consistent with the northward offset found by Tuttle et al. (2013).

### Units and Centroids

Centroids are computed for the 263 countries, 3,662 first-level units and 47,217 second-level units of GADM version 4.1, with the same boundaries for all years (GADM n.d.). We use “region” for the first level and “municipality” for the second; in the United States these are states and counties. Each pixel enters a unit with the weight light × the share of the pixel inside the unit × the pixel's area on the sphere. The weight is linear in light intensity; no threshold or transformation is applied. We compute both the Barmore (1991) centroid, which we report, and the Aboufadel and Austin (2006) centroid; at the first and second levels the two differ by less than 0.01 km in median. At the national level, light north of 65° N is excluded, as in Hall et al. (2019), because high-latitude light is unreliable and its share of light far exceeds its share of population in the Nordic countries.

### Splicing the Sensors

The two sensors give different centroids for the same unit in the same year, because they differ in resolution, saturation and detection threshold. We join them per unit at the overlap year. Each unit's DMSP centroids for 1992–2012 are shifted by the difference between its VIIRS and its DMSP centroid in 2013, in kilometres east and north; VIIRS centroids from 2013 onwards are left as observed. The series is therefore anchored on VIIRS, as PCNL is, but at the level of units rather than pixels, and its join falls between 2012 and 2013. A unit is spliced when it has a centroid from both sensors in 2013: 250 of 263 countries, 3,496 of the 3,661 regions with centroids and 44,019 of 47,217 municipalities. The released series contains the centroids for 1992–2025 with flags for units whose centroid falls outside the unit or that could not be spliced.

Two design choices were made with the U.S. and European benchmarks in view: anchoring the splice on VIIRS rather than on DMSP, which gave smaller location errors, and masking light north of 65° N at the national level. Brazil, added afterwards, is the only benchmark held out from these choices, and it has a single census year. Because each unit's DMSP years are shifted onto its 2013 VIIRS position, the similarity of location errors across census years is partly built into the series.

### Census Benchmarks

Table 1 lists the four benchmarks. For each unit and census year, the census centroid is the Barmore centroid of the benchmark's population points, assigned to GADM units by point-in-polygon. A country is used only if at least 90 percent of its finest GADM units contain census population, and a unit only if it contains at least ten census points. For the 2011 GEOSTAT grid we exclude countries whose grid mixes disaggregated and aggregated sources (Bulgaria and Czechia), because their cell populations are modelled rather than counted. Puerto Rico is reported separately from the United States: its 78 municipios are first-level units in GADM and would otherwise dominate the U.S. first-level statistics. Brazil was added after the methods were fixed and serves as an out-of-sample test of the location results.

| Benchmark | Census years | Population points | Coverage |
| --- | --- | --- | --- |
| U.S. Census Bureau, block-group centres of population | 2000, 2010, 2020 | block-group centre of population | United States, Puerto Rico |
| Eurostat GEOSTAT 1-km grid | 2011, 2021 | populated 1-km cell | 25–26 (2011) and 29–32 (2021) countries and territories, by level |
| INEGI, locality register (ITER) | 2010, 2020 | locality (a large city is a single point) | Mexico |
| IBGE, Grade Estatística 1-km grid | 2022 | populated 1-km cell (1.50 million cells, 203.0 million people) | Brazil |

*Table 1. Census benchmarks.* Sources: U.S. Census Bureau (n.d.); Eurostat (n.d.); INEGI (2011, 2021); IBGE (2025).

### Tests

*Location.* The location error is the great-circle distance between the light-weighted and the census centroid of a unit in a census year, summarized by level and benchmark.

*National bias.* For a unit with population P and light L on a grid of cells, the difference between the light-weighted and the population-weighted mean position equals the population-weighted covariance between position x and light per resident r, divided by mean light per resident (Price 1970, 1972):

```latex
\bar{x}_L - \bar{x}_P = \frac{\mathrm{Cov}_P(x, r)}{\mathrm{E}_P(r)}
```

By the law of total covariance, the right-hand side splits exactly into a between-region term Bb, from differences in mean light per resident between first-level regions, and a within-region term Bw. Two further terms collect light in cells without population (A) and the difference between the light and population centres within cells (C), so that D = A + Bb + Bw + C holds as vectors. We compute the decomposition with the Aboufadel–Austin centroid on a common grid, with raw VIIRS for the year of each benchmark (Mexico 2020, Brazil 2022, Sweden, Norway and Finland 2021); its national errors therefore differ slightly from the Barmore errors of the released series (131 against 133 km for Mexico, 95 against 94 km for Brazil). Light per resident r is computed over populated cells, so its regional values depend on the grid and, for Mexico, on INEGI recording a large city as a single point. The United States is not decomposed. The terms are identities, not causal estimates; how D splits between A and B depends on the cell size, which we report as a scale analysis, and the split between Bb and Bw depends on how the regions are drawn.

*Change.* For each unit observed in two census years, the true shift t is the displacement of the census centroid and the measured shift e that of the light-weighted centroid, both in kilometres east and north. We report the share of units in which e points in the direction of t (positive cosine; 50 percent by chance), the Spearman rank correlation and the ordinary least squares slope between the components of e and t, and the share of units in which the true shift exceeds the error |e − t|. Intervals are 95 percent block-bootstrap intervals (1,000 draws), with the parent region as block for municipalities and the country as block for European regions. We also report the share in the right direction by size class of the true shift. The periods are 2000–2010, 2010–2020 and 2000–2020 for the United States, 2010–2020 for Mexico and 2011–2021 for Europe; Brazil has a single census year and enters only the location tests.

*Comparison series.* To check that the results do not depend on our splice, we run the location and change tests on the series of Chiovelli et al. (2026), which converts VIIRS to DMSP-like values pixel by pixel (variant bloomtopcode\_fix). To describe the sensor join itself, we compute the median year-to-year shift of first-level centroids for 2010–2017 in five series — X. Li et al. (2020) with all pixels and with the recommended threshold DN > 7, Chiovelli et al. (2026), PCNL (S. Li et al. 2023) and ours — on the 2,975 regions with at least ten lit pixels in every year and series, using Aboufadel–Austin centroids. Each series is evaluated at its own join: 2013–2014 for the four published series, which use DMSP through 2013, and 2012–2013 for ours.

## Results

### Location: Good Within Municipalities, Poor for Countries

The location error shrinks by one to two orders of magnitude from countries to municipalities (Figure 1). At the second level the light-weighted centroid lies a median 1.5 km from the census centroid in U.S. counties, 0.7–0.9 km in European municipalities, 1.0–1.1 km in Mexico and 1.1 km in Brazil; nine in ten units are within 2.4–4.3 km. Brazil, tested after the method was fixed, falls within the range of the other benchmarks. At the first level the error depends on the size of the units: 0.3 km for Puerto Rico's municipios, 2.5–3.0 km for European regions, 7–9 km for Mexican states, 11–13 km for U.S. states and 17 km for Brazilian states. At the national level the light-weighted centroid lies 65–88 km from the census centroid in the United States, 133–134 km in Mexico and 94 km in Brazil; the median across European countries is 11–13 km. The errors are stable across census years: in the United States the county median is 1.46–1.52 km in 2000, 2010 and 2020, and the 2000 and 2010 values come from DMSP while the 2020 value comes from VIIRS, although this stability is partly built in by the splice. The series of Chiovelli et al. (2026) gives larger second-level errors: 1.8–2.0 km in U.S. counties, 0.9–1.1 km in Europe, 1.7–1.8 km in Mexico and 1.5 km in Brazil.

### National Bias Follows Light per Resident

The national errors are not random: each points towards the regions where residents produce the most light (Figure 2). In Mexico the light-weighted centroid lies 131 km northwest of the census centroid. Light per resident is 1.7–1.8 times the national mean in the northern border states of Sonora and Chihuahua and 0.5–0.6 times the mean in Guerrero and Mexico City. In Brazil the centroid is displaced 95 km towards the west and south, where Mato Grosso produces 1.7 times the national light per resident and Amazonas 0.6 times. In the Nordic countries, without the 65° N mask, the displacements are 98 km (Sweden), 149 km (Norway) and 71 km (Finland), all northwards: on the 10′ grid, light per resident is 3.1 times the national mean in Norrbotten and 4.0 times in Finnmark, but 0.29 times in Oslo.

The decomposition attributes almost all of the national error to differences between regions (Figure 3a). On a 10′ grid the between-region term Bb has a length of 127 of 131 km in Mexico, 92 of 95 km in Brazil, 95 of 98 km in Sweden and 146 of 149 km in Norway. Finland is the exception: Bb is 38 km and the within-region term Bw 34 km of 71 km, because Finland's five first-level regions are large and internally heterogeneous. Light in unpopulated cells (A) and the offset within cells (C) contribute less than 2 km on this grid. The split is not invariant to the grid (Figure 3b). With 1′ cells, many cells contain light but no recorded population, and A takes up 115 of Mexico's 131 km; with 5′ cells and above, nearly all light falls in populated cells and Bb dominates. D itself does not depend on the grid. The robust statement is therefore that the national error equals the covariance between position and light per resident, and that at the scale of regions this covariance is dominated by between-region contrasts. Masking light north of 65° N reduces the Nordic errors (Sweden 98 → 33 km, Norway 149 → 61 km, Finland 71 → 42 km on the same grid) without removing them. In Sweden the remaining error still points north; in Norway it reverses and points south, because residents north of 65° N keep their weight while their light is removed; in Finland it lies mainly within regions. The size of the national bias also depends on the light product. For Mexico the national error is 133–134 km in our series, 95–119 km in the series of Chiovelli et al. (2026) and 35–42 km in a version of our series in which VIIRS is first converted to DMSP-like values; for Brazil the corresponding errors are 94, 111 and 73 km. The 65–134 km bias is therefore a property of raw VIIRS radiance as weights, and the dominance of between-region contrasts is a statement about GADM first-level regions.

### Small Shifts Are Not Detectable

Over 10–20 years, the measured shift of a municipal centroid is only weakly related to the true shift (Figure 4; Table 2). The true shifts are small: a median 0.13–0.54 km per municipality. The share of municipalities in which the measured shift points the right way is 66 percent for U.S. counties over 2000–2020, 58–61 percent over each decade, 56 percent in Mexico and 51 percent in Europe, against 50 percent by chance. Rank correlations between measured and true shifts are 0.02–0.34. In 72–92 percent of municipalities the error of the measured shift exceeds the true shift, so that assuming no movement at all comes closer to the census than the light-measured shift. Where a relation exists, the slope is 0.3–0.5: lights capture less than half of the true movement.

At the first level the result depends on how far the centres moved. Mexico's 32 states, whose centres moved a median 2.1 km, are measured in the right direction in 72 percent of cases (95 percent interval 56–88 percent), with a rank correlation of 0.42. European regions, which moved a median 0.42 km, are at 49 percent, and their rank correlation is slightly negative (−0.16 to −0.01). With the series of Chiovelli et al. (2026), European municipalities are also below chance (45 percent, 43–46). We read these below-chance results as a sign of differences between the 2011 and 2021 GEOSTAT grids rather than as a property of lights.

| Benchmark | Level | Period | Units | True shift, median (km) | Right direction | Rank correlation | Slope | True shift > error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| United States | 2 | 2000–2010 | 2,381 | 0.37 | 0.61 (0.58–0.64) | 0.24 (0.18–0.29) | 0.45 (0.30–0.60) | 26% |
| United States | 2 | 2010–2020 | 2,383 | 0.26 | 0.58 (0.56–0.60) | 0.18 (0.14–0.21) | 0.49 (0.19–0.77) | 16% |
| United States | 2 | 2000–2020 | 2,373 | 0.54 | 0.66 (0.64–0.69) | 0.34 (0.29–0.38) | 0.48 (0.36–0.59) | 28% |
| Mexico | 1 | 2010–2020 | 32 | 2.07 | 0.72 (0.56–0.88) | 0.42 (0.09–0.73) | 0.38 (−0.26 to 2.44) | 19% |
| Mexico | 2 | 2010–2020 | 1,915 | 0.25 | 0.56 (0.53–0.58) | 0.13 (0.08–0.18) | 0.26 (0.02–0.86) | 14% |
| Europe | 1 | 2011–2021 | 329 | 0.42 | 0.49 (0.43–0.53) | −0.08 (−0.16 to −0.01) | −0.14 (−0.70 to 0.53) | 5% |
| Europe | 2 | 2011–2021 | 6,161 | 0.13 | 0.51 (0.49–0.53) | 0.02 (−0.01 to 0.05) | 0.22 (−0.03 to 0.47) | 8% |

*Table 2. Detection of change by light-weighted centroids. 95 percent block-bootstrap intervals in parentheses. U.S. first-level units are omitted because the census sample at that level is dominated by Puerto Rico's municipios.*

Detection improves with the size of the shift but becomes reliable in only one benchmark (Figure 5). Among U.S. counties over 2000–2020, the share in the right direction rises from 54 percent for true shifts below 0.1 km to 69 percent for 0.5–1 km and 84 percent for 2–5 km; only 8 percent of counties moved more than 2 km. Over the single decades the 2–5 km class reaches 77 percent (2000–2010) and 70 percent (2010–2020). In Mexico the share stays below 70 percent in every size class with at least 20 units, including 58 percent for 2–5 km, and in Europe it stays at 50–52 percent in every class up to 2 km. The weak result is not specific to the sensor join: over 2000–2010, which lies entirely in the DMSP era, the share in the right direction for U.S. counties is 61 percent, against 58 percent for 2010–2020, which spans the join. The intervals in Figure 5 are binomial and do not account for spatial clustering; those in Table 2 do.

The result does not depend on our splice. With the series of Chiovelli et al. (2026), the share of municipalities in the right direction is 0.45–0.62 and the rank correlation −0.12 to 0.25; for 2010–2020 in the United States and in Mexico it is no better than chance.

### The Sensor Join Per Unit

Harmonized series that are continuous in aggregate are not necessarily continuous unit by unit (Figure 6). On 2,975 regions common to all series, the median shift of the centroid at the join is 3.3 times the median DMSP-era shift in X. Li et al. (2020) with all pixels, 1.6 times with the threshold DN > 7 and 1.5 times in Chiovelli et al. (2026). In PCNL and in our series, which anchor on VIIRS in the overlap year, the ratio is 0.8 and 0.9. For these two series, continuity at the join follows from the construction and says nothing about accuracy, which only the census tests address. All VIIRS-based series show a larger shift from 2016 to 2017, about 1.4 times a normal VIIRS year in our series and larger than the join itself in X. Li et al. (2020) with DN > 7; we report it but have not traced its source.

## Discussion

### A Good Snapshot and a Poor Film

Nighttime lights give a good snapshot of where people live within regions and a poor film of how that changes. The two results are consistent. A location error of about 1 km is small relative to the size of a municipality but large relative to how far a municipality's centre of population moves in a decade, which is 0.13–0.54 km in the median municipality of the three benchmarks with two census years. A unit's location errors are strongly correlated between census years: the median error of a measured shift is only 0.5–0.8 times the median location error, which implies a correlation of roughly 0.7–0.9 if the errors are normally distributed. The uncorrelated remainder is nonetheless larger than most true shifts. We interpret it as year-to-year noise in the light signal together with changes in light per resident that are unrelated to where people moved, among them the change to LED street lighting, to which the VIIRS day/night band is only partly sensitive \[REF needed\]. Separating these sources would require light and population at the same fine resolution in several years.

The implication for practice is direct. Studies that read 10–30-year movements of light-weighted centres or ellipses as shifts of population at regional scale measure noise as much as change, and the problem is not confined to the sensor join. Studies that track the centre of light itself, as a measure of urban expansion, are not tested here. Where population shifts are large the signal can exceed the noise, as for Mexican states; our results suggest that this requires shifts of kilometres per unit, and that the threshold must be shown for the region in question rather than assumed. The results also qualify the proposal of Hall et al. (2019) that light-weighted centroids can be used to analyse changes in population distribution. They do not affect cross-sectional uses such as the migration distances of Niedomysl et al. (2017), which are computed between centroids in a single year; the location error there was a few kilometres for Swedish labour-market regions, and it is about 1 km for the municipalities tested here.

### The Level–Change Gap Extends to Location

Our results extend a pattern documented for economic activity (Goldblatt, Heilmann, and Vaizman 2020; Asher et al. 2021) and for gridded population counts (Zhang and Gibson 2025; L. Li et al. 2026) to the location of population. The similarity with L. Li et al. (2026) is close, although a change in counts and a change in location are not directly comparable: six gridded products got the direction of township population change right in 53 percent of cases, and light-weighted centroids get the direction of municipal centroid shifts right in 51–66 percent. Many gridded population products use nighttime lights as a covariate, and in WorldPop's models lights are the most important predictor of population density (Nieves et al. 2021). If light per resident varies by region, as our decomposition shows, and with income, as Mirza et al. (2021) argue, gridded products may inherit the same regional bias. We have not tested this.

### Why National Centroids Fail

The national bias is a property of the weights, not of the centroid formula. The Barmore and Aboufadel–Austin centroids differ by less than 0.01 km in median at the first and second levels, and their national errors for Mexico and Brazil differ by less than 0.2 km. The bias is the covariance between position and light per resident, an identity known from Price (1970, 1972) that we apply to light. What our decomposition adds is that, for first-level regions, this covariance is dominated by between-region contrasts: in the countries decomposed, sparsely populated northern and agricultural regions with high light per resident pull the centre of light away from dense metropolitan regions with low light per resident. The identity also points to a correction. Weighting each region's light-weighted centroid by the region's census population, rather than by its light, removes the between-region term by construction and leaves only the within-region term. Regional population totals are available almost everywhere, even where fine census data are not. We have not tested this correction.

### Recommendations

For users of light-weighted centroids, three recommendations follow.

1. Use centroids at the first and second levels for location in a given year — distances, accessibility and gravity models. Use them for change only where the expected shift is kilometres rather than hundreds of metres, and check the detection threshold against census data in the region of interest.
2. Treat national centroids as biased towards regions with high light per resident; with raw VIIRS weights the bias is 65–134 km in the United States, Mexico and Brazil. Where regional population totals exist, reweighting regional centroids by population is a simple remedy to test.
3. For per-unit measures across the sensor change, prefer series anchored in the overlap year, per unit or per pixel. We have tested this choice for continuity at the join and against one alternative series, not against an independent census benchmark.

For producers of harmonized series, evaluation per unit should complement evaluation in aggregate. A series without a jump in aggregate totals can still have a join at which the centroid of a typical region moves 1.5–3.3 times more than in a normal year, a point also made by X. Li et al. (2026).

### Limitations

Five limitations bound the conclusions. First, change entirely within the VIIRS era is not tested, because no benchmark has two census years after 2013; the conclusions apply to DMSP-era and spliced series, and VIIRS alone may do better. Second, Europe is the weakest benchmark: its true shifts are small (median 0.13 km at the second level), and GEOSTAT 2011 and 2021 differ in method, so part of the “true” shift may be a difference between grids; the United States and Mexico, with the same type of census data in both years, carry the evidence. Third, there is no benchmark in Africa or Asia, where population redistribution is fastest and true shifts may exceed the detection threshold more often. Fourth, the VIIRS annual composites show a larger year-to-year change from 2016 to 2017 in all series, which we report but do not correct. Fifth, the splice direction and the national mask were chosen with the U.S. and European benchmarks in view, and only Brazil is held out; the location results for those benchmarks are therefore not fully independent of the design.

## Conclusion

Light-weighted centroids locate the population of a county or municipality to within about 1 km in four census systems, including one held out from the design of the series. National centroids are biased by tens to more than a hundred kilometres; the bias equals the covariance between position and light per resident, and in the countries decomposed it is dominated by differences between first-level regions. Change is detected only where it is large. The shifts of a few hundred metres that dominate municipal redistribution in the United States, Europe and Mexico are measured in the right direction little more often than chance, and in most municipalities assuming no movement comes closer to the census; the larger shifts of Mexican states are detected in about three cases out of four. The level–change gap known for economic activity and population counts therefore also holds for the location of population. Light-weighted centroids are a sound tool for describing where people live in a given year in places without fine census data, and a weak one for describing where they move.

## Data Availability and Figure Captions

**Data availability.** The centroid series (version 2.0.0, 1992–2025, three GADM levels, Barmore and Aboufadel–Austin centroids with flags) is deposited at Zenodo \[DOI to be added\] under CC BY-NC 4.0, which follows the licence of the flare-site catalogue used for masking. Input and comparison data: VIIRS Nightfire flare-site catalogue MYC25Q4 (doi:10.5281/zenodo.23091822), the series of Chiovelli et al. (doi:10.5281/zenodo.20558766) and PCNL (doi:10.5281/zenodo.7612389). Code, validation tables and the scripts that draw every figure are available at \[repository, anonymized for review\] under the MIT licence.

**Figure 1.** Distance between the light-weighted and the census population centroid per unit, by administrative level and benchmark, latest census year. Bars show medians, boxes the 25th–75th percentiles and lines the 10th–90th percentiles; groups with fewer than 20 units are shown as points.

**Figure 2.** Light per resident by first-level region, relative to the national mean (10′ grid), and the displacement of the national light-weighted centroid from the census centroid (Aboufadel–Austin; arrows drawn at three times their length; labels give the true distance). VIIRS 2020 (Mexico), 2022 (Brazil) and 2021 (Sweden, Norway, Finland); no masking north of 65° N.

**Figure 3.** (a) Decomposition of the national error D into light in unpopulated cells (A), between-region (Bb) and within-region (Bw) differences in light per resident, 10′ grid; component lengths in km. (b) Mexico: the components by grid cell size.

**Figure 4.** Light-measured shift of the centroid along the direction of the true shift, against the true shift of the census centroid, for counties in the United States (2000–2020) and municipalities in Mexico (2010–2020) and Europe (2011–2021). Points on the line y = x are measured perfectly; points below zero point the wrong way. Values beyond ±4 km are drawn at the edge.

**Figure 5.** Share of municipalities and counties in which the measured shift points in the right direction, by size class of the true shift, with 95 percent Wilson intervals. Classes with fewer than 20 units are omitted.

**Figure 6.** Median year-to-year shift of first-level centroids (Aboufadel–Austin), 2010–2017, in five series on 2,975 common regions. Circles mark each series' join: 2013–2014 for the published series, 2012–2013 for ours. The ratio in the legend is the join step divided by the median of the DMSP-era steps before it.

## References

*Every entry below has been checked against the publisher, PubMed Central, EconPapers or a library record in this project (docs/LITERATURE.md). Titles, journals, volumes and DOIs are checked; authors' initials are not yet checked against the sources, and the items marked in brackets must be completed before submission. Data citations (GADM, VNL v2, MYC25Q4, the four census benchmarks) are still to be added.*

Aboufadel, E., and D. Austin. 2006. A new method for computing the mean center of population of the United States. *The Professional Geographer* 58 (1): 65–69.

Archila Bustos, M. F., O. Hall, and M. Andersson. 2015. Nighttime lights and population changes in Europe 1992–2012. *Ambio* 44 (7): 653–65. doi:10.1007/s13280-015-0646-8.

Asher, S., T. Lunt, R. Matsuura, and P. Novosad. 2021. Development research at high geographic resolution. *The World Bank Economic Review* 35 (4): 845–71. doi:10.1093/wber/lhab003.

Ash, K., and E. Mazur. 2020. Identifying and correcting signal shift in DMSP-OLS data. *Remote Sensing* 12 (14): 2219. doi:10.3390/rs12142219.

Barmore, F. E. 1991. Where are we? Comments on the concept of the “center of population.” *The Wisconsin Geographer* 7:40–50.

Bennett, M. M., and L. C. Smith. 2017. Advances in using multitemporal night-time lights satellite imagery to detect, estimate, and monitor socioeconomic dynamics. *Remote Sensing of Environment* 192:176–97. doi:10.1016/j.rse.2017.01.005.

Bluhm, R., and G. C. McCord. 2022. What can we learn from nighttime lights for small geographies? Measurement errors and heterogeneous elasticities. *Remote Sensing* 14 (5): 1190. doi:10.3390/rs14051190.

Chiovelli, G., S. Michalopoulos, E. Papaioannou, and T. Regan. 2026. Illuminating the Global South. *The Economic Journal*, advance article \[volume and pages when issued\]. doi:10.1093/ej/ueaf134.

Dong, H., R. Li, J. Li, and S. Li. 2020. Study on urban spatiotemporal expansion pattern of three first-class urban agglomerations in China derived from integrated DMSP-OLS and NPP-VIIRS nighttime light data \[in Chinese\]. *Journal of Geo-information Science* 22 (5): 1161–74. doi:10.12082/dqxxkx.2020.190711.

Elvidge, C. D., D. Ziskin, K. E. Baugh, B. T. Tuttle, T. Ghosh, D. W. Pack, E. H. Erwin, and M. Zhizhin. 2009. A fifteen year record of global natural gas flaring derived from satellite data. *Energies* 2 (3): 595–622. doi:10.3390/en20300595.

Elvidge, C. D., M. Zhizhin, T. Ghosh, F.-C. Hsu, and J. Taneja. 2021. Annual time series of global VIIRS nighttime lights derived from monthly averages: 2012 to 2019. *Remote Sensing* 13 (5): 922. doi:10.3390/rs13050922.

Eurostat. n.d. GEOSTAT population grid 2011 (version 2.0.1) and Census 2021 population grid (version 3). GISCO, European Commission. https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/population-grids.

GADM. n.d. GADM database of global administrative areas, version 4.1. https://gadm.org.

Gaughan, A. E., T. Oda, A. Sorichetta, F. R. Stevens, M. Bondarenko, R. Bun, L. Krauser, G. Yetman, and S. V. Nghiem. 2019. Evaluating nighttime lights and population distribution as proxies for mapping anthropogenic CO2 emission in Vietnam, Cambodia and Laos. *Environmental Research Communications* 1:091006. doi:10.1088/2515-7620/ab3d91.

Goldblatt, R., K. Heilmann, and Y. Vaizman. 2020. Can medium-resolution satellite imagery measure economic activity at small geographies? Evidence from Landsat in Vietnam. *The World Bank Economic Review* 34 (3): 635–53. doi:10.1093/wber/lhz001.

Grether, J.-M., and N. A. Mathys. 2011. On the track of the world's economic center of gravity. In *Economic Growth and Development*, Frontiers of Economics and Globalization 11, 261–87. Bingley, UK: Emerald. doi:10.1108/S1574-8715(2011)0000011015.

Hall, O., M. F. Archila Bustos, N. Boke Olén, and T. Niedomysl. 2019. Population centroids of the world administrative units from nighttime lights 1992–2013. *Scientific Data* 6:235. doi:10.1038/s41597-019-0250-z.

IBGE (Instituto Brasileiro de Geografia e Estatística). 2025. Grade estatística, Censo Demográfico 2022, 1 km (BR1KM\_20251002). https://geoftp.ibge.gov.br/recortes\_para\_fins\_estatisticos/grade\_estatistica/censo\_2022/grade\_1km/.

INEGI (Instituto Nacional de Estadística y Geografía). 2011. Censo de Población y Vivienda 2010: Principales resultados por localidad (ITER). https://www.inegi.org.mx.

INEGI (Instituto Nacional de Estadística y Geografía). 2021. Censo de Población y Vivienda 2020: Principales resultados por localidad (ITER). https://www.inegi.org.mx/programas/ccpv/2020/.

Li, L., S. Fu, X. Zhou, K. Xiao, X. Cao, B. Zhang, F. Li, et al. 2026. Cross-sectional accuracy does not imply the reliability of population change in gridded population datasets of China. *Humanities and Social Sciences Communications* 13 (1) \[article number\]. doi:10.1057/s41599-026-07688-w.

Li, S., X. Cao, C. Zhao, N. Jie, L. Liu, X. Chen, and X. Cui. 2023. Developing a pixel-scale corrected nighttime light dataset (PCNL, 1992–2021) combining DMSP-OLS and NPP-VIIRS. *Remote Sensing* 15 (16): 3925. doi:10.3390/rs15163925.

Li, X., Y. Lei, J. Shen, R. Xin, and K. Qu. 2026. Spatiotemporal dynamics and relative growth of county-level CO2 emissions in China (1997–2021) derived from harmonized nighttime lights. *Environmental Monitoring and Assessment* 198 (12): 1287. doi:10.1007/s10661-026-16129-9.

Li, X., Y. Zhou, M. Zhao, and X. Zhao. 2020. A harmonized global nighttime light dataset 1992–2018. *Scientific Data* 7:168. doi:10.1038/s41597-020-0510-y.

Mirza, M. U., C. Xu, B. van Bavel, E. H. van Nes, and M. Scheffer. 2021. Global inequality remotely sensed. *Proceedings of the National Academy of Sciences* 118 (18): e1919913118. doi:10.1073/pnas.1919913118.

Niedomysl, T., O. Hall, M. F. Archila Bustos, and U. Ernstson. 2017. Using satellite data on nighttime lights intensity to estimate contemporary human migration distances. *Annals of the American Association of Geographers* 107 (3): 591–605. doi:10.1080/24694452.2016.1270191.

Nieves, J. J., M. Bondarenko, D. Kerr, N. Ves, G. Yetman, P. Sinha, D. J. Clarke, A. Sorichetta, F. R. Stevens, A. E. Gaughan, and A. J. Tatem. 2021. Measuring the contribution of built-settlement data to global population mapping. *Social Sciences & Humanities Open* 3:100102. doi:10.1016/j.ssaho.2020.100102.

Plane, D. A., and P. A. Rogerson. 2015. On tracking and disaggregating center points of population. *Annals of the Association of American Geographers* 105 (5): 968–86. doi:10.1080/00045608.2015.1066742.

Price, G. R. 1970. Selection and covariance. *Nature* 227:520–21. doi:10.1038/227520a0.

Price, G. R. 1972. Extension of covariance selection mathematics. *Annals of Human Genetics* 35 (4): 485–90.

Quah, D. 2011. The global economy's shifting centre of gravity. *Global Policy* 2 (1): 3–9.

Tuttle, B. T., S. J. Anderson, P. C. Sutton, C. D. Elvidge, and K. Baugh. 2013. It used to be dark here: Geolocation calibration of the Defense Meteorological Satellite Program Operational Linescan System. *Photogrammetric Engineering & Remote Sensing* 79 (3): 287–97. doi:10.14358/PERS.79.3.287.

U.S. Census Bureau. n.d. Centers of population: Block group, 2000, 2010 and 2020. https://www.census.gov/geographies/reference-files/time-series/geo/centers-population.html.

Wang, J., P. Hou, Y. Liu, J. Zhai, Y. Chen, Y. Liu, L. Xie, and X. Ren. 2026. Gradient characteristics and spatiotemporal differentiation of China's urbanization development based on nighttime light remote sensing data \[in Chinese\]. *National Remote Sensing Bulletin* 30 (1): 61–78. doi:10.11834/jrs.20254567.

Zhang, X., and J. Gibson. 2025. How well do gridded population estimates proxy for actual population changes? Evidence from China. *Population and Development Review* 51 (4): 1402–34. doi:10.1111/padr.70025.

Zhizhin, M., C. D. Elvidge, T. Ghosh, G. Gleason, and M. Bazilian. 2026. VIIRS Nightfire super-resolution method for multiyear cataloging of natural gas flaring sites: 2012–2025. *Remote Sensing* 18 (2): 314. doi:10.3390/rs18020314.
