#!/usr/bin/env python3
"""Bilagetabeller S1–S6 till manuset, på engelska, skrivna ur docs/extended/
(och docs/ANALYSIS_PLAN.md för S6). Utdata: docs/manuscript/supplement_tables.md.

S1 förändring per skattning (A2); S2 modellen per enhet (A3); S3 regionala totaler
(A1); S4 jämförelsekällor med bortfall (B); S5 Kanada och Japan med tolkningsregeln
(C); S6 avvikelserna ur planen.

Alla tal läses ur CSV-filerna. Den engelska texten i S6 är en översättning av
planens avvikelsetabell; skriptet stannar om planen har fått rader som inte är
översatta (datum och avsnitt jämförs rad för rad).

Kör från repots rot: python tools/write_supplement_tables.py"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

E = Path("docs/extended")
PLAN = Path("docs/ANALYSIS_PLAN.md")
OUT = Path("docs/manuscript/supplement_tables.md")

NAME = {"census": "United States", "geostat": "Europe", "inegi": "Mexico", "ibge": "Brazil",
        "canada": "Canada", "japan": "Japan"}
LEVEL = {1: "1", 2: "2"}
SOURCE = {"light": "Light", "ghsbuilts": "Built-up surface", "ghspop": "GHS-POP",
          "worldpop": "WorldPop", "geometric": "Geometric centre"}
EST_A2 = {"endpoint": "Endpoints", "three_year": "Three-year means", "trend": "Trend"}
EST_A1 = {"light": "Light", "hybrid_light": "Hybrid, light", "hybrid_geom": "Hybrid, geometric"}
PAIR_A1 = {"light_minus_hybrid_light": "Light − hybrid, light",
           "hybrid_geom_minus_hybrid_light": "Hybrid, geometric − hybrid, light"}
MODEL = {"snr": "log(|t|/σ)", "abs_t": "log |t|"}


# --- formatering ------------------------------------------------------------

def num(x, d: int = 2) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    s = f"{x:,.{d}f}"
    return s.replace("-", "−")


def ci(v, lo, hi, d: int = 2) -> str:
    """'0.61 (0.58–0.64)'; med 'to' när en gräns är negativ, som i manuset."""
    sep = " to " if (lo < 0 or hi < 0) else "–"
    return f"{num(v, d)} ({num(lo, d)}{sep}{num(hi, d)})"


def md_table(header, rows) -> str:
    esc = lambda c: str(c).replace("|", "\\|")
    out = ["| " + " | ".join(esc(h) for h in header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    out += ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def period(a, b) -> str:
    return f"{int(a)}–{int(b)}"


def read(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


# --- tolkningsregeln (planen §7; rad 4 tillagd i efterhand) -----------------

RULE = [
    (1, "Estimate ≥ 0.65 and lower bound > 0.50",
     "VIIRS measures direction better than the DMSP-era and spliced periods (0.51–0.66).", False),
    (2, "Estimate 0.50–0.65 and interval width ≤ 0.15",
     "The result from the DMSP era also holds for VIIRS.", False),
    (3, "Interval width > 0.15 or fewer than 150 units",
     "Inconclusive. A second benchmark within the VIIRS era is added under the same plan, chosen before its data are opened.", False),
    (4, "Estimate < 0.50 and interval width ≤ 0.15",
     "VIIRS does not measure direction in this benchmark. If the upper bound is below 0.50, the outcome is reported as "
     "significantly below chance and its cause as not investigated.", True),
]


def rule_outcome(est, lo, hi, n):
    """Raden i tolkningstabellen som utfallet hör till (1–4), annars None.
    Färre än 150 enheter är inkonklusivt (rad 3) oavsett punktskattning."""
    w = hi - lo
    if n < 150:
        return 3
    if est >= 0.65 and lo > 0.50:
        return 1
    if 0.50 <= est < 0.65 and w <= 0.15:
        return 2
    if w > 0.15 or n < 150:
        return 3
    if est < 0.50 and w <= 0.15:
        return 4
    return None


# --- avvikelser (S6) ----------------------------------------------------------

def parse_deviations(text: str):
    """(datum, avsnitt) för varje rad i planens tabell under '## 10. Avvikelser'."""
    part = text.split("## 10. Avvikelser", 1)[1]
    rows = []
    for line in part.splitlines():
        m = re.match(r"^\|\s*(\d{4}-\d{2}-\d{2})\s*\|\s*([^|]+?)\s*\|", line)
        if m:
            rows.append((m.group(1), m.group(2)))
    return rows


def check_deviations(rows, en):
    if [(d, s) for d, s, *_ in en] != list(rows):
        raise RuntimeError("S6: planens avvikelsetabell och den engelska översättningen stämmer inte "
                           f"rad för rad ({len(rows)} rader i planen, {len(en)} översatta). "
                           "Översätt de nya raderna i DEVIATIONS_EN.")


# Översättning av docs/ANALYSIS_PLAN.md §10, i planens ordning (datum, avsnitt, ändring, skäl).
DEVIATIONS_EN = [
    ("2026-11-21", "§2",
     "scipy added to requirements.txt and installed.",
     "The frozen measures in tools/change_detection.py use pandas' Spearman correlation, which requires scipy under pandas ≥ 3. No change to code or measures."),
    ("2026-11-21", "§6",
     "\"Same centroid code\" made precise: comparison sources whose cell value is a total per cell (GHS-POP population, GHS-BUILT-S built-up surface) are converted to density (value/cell area), so that the weight becomes total × share. wpc/centroids.py is not changed.",
     "The centroid code weights by value × share × area, which is right for radiance but not for totals per cell; the conversion gives the total-weighted centroid (tested in tests/test_ghspop.py)."),
    ("2026-11-21", "§6",
     "Caveat extended: WorldPop uses nighttime lights as a covariate and is not independent of light; GHS-POP and GHS-BUILT-S rest on census and building data respectively. WorldPop is reported separately from the other two.",
     "Independence from the benchmark differs between sources and must be stated in the text."),
    ("2026-11-21", "§6",
     "Exploratory addition (after the GHS-BUILT-S results had been seen): paired location comparison per unit between light and each source, levels 1 and 2 — share of units in which the source is closer to the census than light, and median paired difference, with block bootstrap.",
     "Median errors hide how often a source is closer in the individual unit; the paired comparison shows this and holds the units equal between source and light."),
    ("2026-11-21", "§7",
     "Freeze point: the scripts for A1–A3 and B were committed at 2e5d99c. Canadian benchmark data could be downloaded and opened only after this commit.",
     "Condition in §7: the same scripts are run unchanged on Canada; the commit hash is recorded here."),
    ("2026-11-21", "§7",
     "Level and point type for Canada: dissemination area (DA), DA representative point, DA population (sum of its blocks) as weight, in both years.",
     "Fallback rule in §7: the Geographic Attribute File has no dissemination-block points in either year, only DA points; hence DA in both years."),
    ("2026-11-21", "§7",
     "C is run with a separate driver, tools/canada_confirm.py, that imports the frozen functions and applies them to Canada. The frozen scripts are not changed. The driver reproduces the committed U.S. county results for 2010–2020 exactly (n = 2,383, right direction 0.582, median |e − t| 0.7947).",
     "The frozen scripts' lists of datasets and years do not include Canada; running them unchanged would otherwise require edits. The functions (the measures) are reused instead, which is the intent of \"the same scripts unchanged\"."),
    ("2026-11-21", "§7",
     "Benchmark data stored in data/validation_canada/ (not in the frozen data/validation/). About 0.9 percent of Canada's population lies in DA points outside the GADM polygons (coastal points).",
     "The freeze in §2 protects data/validation/. The points outside do not materially affect the selection but are reported."),
    ("2026-11-21", "§6/§7",
     "Epochs of the comparison sources for Canada: the epoch nearest 2016 is 2015, not 2020. GHS-POP and GHS-BUILT-S use epochs 2015 and 2020; WorldPop uses 2016 and 2020. Light uses exactly 2016 and 2021.",
     "§6: nearest epoch. The sources' periods (2015–2020 and 2016–2020) differ from the census period 2016–2021; stated in the results document."),
    ("2026-11-21", "§7",
     "The A3 secondary measure is not reported for Canada.",
     "The frozen A3 model's σ is not defined for a period entirely within VIIRS (2016–2021), and the model failed its own test (H3b, H3c). A σ for VIIRS would be an analysis beyond the plan."),
    ("2026-11-21", "§7",
     "Second benchmark within the VIIRS era (row 3 of the rule). Choice: Japan (censuses 2015 and 2020) if population is enumerated (not modelled) at a fine level with the same method and the same type of point or cell in both years; otherwise Australia (2016 and 2021) on the same terms. The second benchmark is tested alone against the same rule; if it falls in row 1 or 2 that conclusion holds and Canada is reported alongside; if it also falls in row 3 the VIIRS question is reported as unanswered. Stopping rule: no third benchmark.",
     "Canada was inconclusive (interval width 0.17 > 0.15); the choice and the stopping rule were recorded before any benchmark data were checked."),
    ("2026-11-21", "§7",
     "Second benchmark chosen: Japan, e-Stat grid statistics from the census (enumerated counts), 2015 and 2020. Blocks = prefectures. Australia was not checked.",
     "Checked at the source: enumerated counts, 1 km (and 500 m) mesh, the same mesh system in both years, licence compatible with CC BY 4.0."),
    ("2026-11-21", "§7",
     "Cell size: 500 m (half mesh) if available with enumerated counts in the same form in both years, otherwise 1 km. The cell centre (from the mesh code) as point, the cell population as weight.",
     "Finest comparable level in both years; point type and method the same in 2015 and 2020."),
    ("2026-11-21", "§7",
     "Disclosure control: cells that are suppressed or merged are documented and treated alike in both years (non-numeric population dropped; merged counts stay in the destination cell).",
     "Disclosure rules must not differ between years and bias the change."),
    ("2026-11-21", "§6/§7",
     "Comparison sources for Japan: GHS-POP and GHS-BUILT-S epochs 2015 and 2020; WorldPop 2015 and 2020. Light on the census years 2015 and 2020.",
     "The census years coincide with the GHS epochs and WorldPop's annual values, so there is no epoch offset (unlike Canada)."),
    ("2026-11-21", "§7",
     "Caveat in the results document: 2020 was a pandemic year, which may affect the census and the comparability of 2015 → 2020.",
     "Reported in docs/extended/C_japan.md."),
    ("2026-11-21", "§7",
     "The interpretation rule lacked the outcome of an estimate below 0.50 with a narrow interval. Row 4 added after the Japanese result had been seen: estimate < 0.50 and interval width ≤ 0.15 → VIIRS does not measure direction in this benchmark; if the upper bound is below 0.50, the outcome is reported as significantly below chance and its cause as not investigated.",
     "Japan 2015–2020 gave a right-direction share of 0.46 with a narrow interval (width 0.05), which none of rows 1–3 covered."),
    ("2026-11-22", "§6",
     "Error in B corrected (found after the review): tools/compare_benchmarks.py used the first source's (GHS-POP's) epochs for all sources. WorldPop was therefore tied to 2010 for GEOSTAT 2011, whose file lacks Europe, so the common set emptied and Europe 2011 dropped out of location and change. Corrected so that each source uses its own epochs; units lost are reported per benchmark and year, and a test requires every benchmark year and period to be present for every source with data. canada_confirm.py did not have the error. Europe's 2011 values for light, GHS-POP and GHS-BUILT-S are the same as before the error; other B values unchanged.",
     "The sources have different epochs (WorldPop has 2011 and 2016, which GHS lacks); the assumption of shared epochs failed when WorldPop was added."),
    ("2026-11-22", "§2 (manus)",
     "Sensitivity of the slope added as a descriptive analysis: OLS slope without the three most influential units and Theil–Sen slope, level 2 (tools/manuscript_numbers.py). Also scripted numbers for Sweden, Fisher's test Mexico–United States, FIPS against GADM assignment, and light against the geometric centre.",
     "Number audit of the manuscript: one sentence rested on a selection outside the repository and several numbers lacked scripts. Exploratory; changes no hypothesis outcome."),
]


# --- tabeller -----------------------------------------------------------------

def change_cols(r):
    return [num(r.n, 0), num(r.true_med), ci(r.err_med, r.err_med_lo, r.err_med_hi),
            ci(r.dir_right, r.dir_right_lo, r.dir_right_hi), ci(r.rho, r.rho_lo, r.rho_hi)]


def s1():
    rows = []
    for ds in ("census", "geostat", "inegi"):
        d = pd.read_csv(E / f"{ds}_trend_change.csv")
        for _, r in d.sort_values(["level", "from", "to"], kind="stable").iterrows():
            rows.append([NAME[ds], LEVEL[r.level], period(r["from"], r.to), EST_A2[r.series]] + change_cols(r))
    return ("## Table S1. Change by estimator\n\n"
            + md_table(["Benchmark", "Level", "Period", "Estimator", "Units", "True shift, median (km)",
                        "Error |e − t|, median (km)", "Right direction", "Rank correlation"], rows)
            + "\n\n*Endpoints: single-year centroids at both census years. Three-year means: centroids averaged over "
              "the census year and the years on either side. Trend: Theil–Sen slope within each sensor period "
              "(DMSP to 2012, VIIRS from 2013), at least four years in each part. Same units for the three estimators. "
              "95 percent block-bootstrap intervals in parentheses. Not computed for Brazil (one census year), "
              "Canada or Japan. Source: docs/extended/<benchmark>_trend_change.csv.*\n")


def s2():
    lab = lambda ds, a, b: f"{NAME[ds]} {period(a, b)}"
    rho = pd.read_csv(E / "detection_rho.csv")
    t1 = md_table(["Benchmark and period", "Units", "ρ(σ, |e − t|)"],
                  [[lab(r.dataset, r["from"], r.to), num(r.n, 0), ci(r.rho, r.rho_lo, r.rho_hi)] for _, r in rho.iterrows()])
    mo = pd.read_csv(E / "detection_models.csv")
    t2 = md_table(["Model", "Fitted on", "Units", "β0", "β1"],
                  [[MODEL[r.model], lab(r.fit_dataset, r["from"], r.to), num(r.n, 0), num(r.beta0, 3), num(r.beta1, 3)]
                   for _, r in mo.iterrows()])
    te = pd.read_csv(E / "detection_test.csv")
    t3 = md_table(["Tested on", "Model", "Units", "Base rate", "AUC", "Brier", "Brier, constant"],
                  [[lab(r.test_dataset, r["from"], r.to), MODEL[r.model], num(r.n, 0), num(r.base_rate),
                    num(r.auc, 3), num(r.brier, 3), num(r.brier_const, 3)] for _, r in te.iterrows()])
    ca = pd.read_csv(E / "detection_calibration.csv")
    t4 = md_table(["Tested on", "Class", "Units", "Predicted", "Observed", "Difference"],
                  [[lab(r.test_dataset, r["from"], r.to), int(r.bin), num(r.n, 0), num(r.pred_mean),
                    num(r.obs_prop), num(r.abs_gap)] for _, r in ca.iterrows()])
    mm = pd.read_csv(E / "detection_min_move.csv")
    t5 = md_table(["Benchmark and period", "Units", "σ, median (km)", "Minimum shift, median (km)", "25th pct.", "75th pct."],
                  [[lab(r.dataset, r["from"], r.to), num(r.n, 0), num(r.sigma_med), num(r.min_move_med),
                    num(r.min_move_q25), num(r.min_move_q75)] for _, r in mm.iterrows()])
    return ("## Table S2. Per-unit model for the right direction\n\n"
            "σ is the unit's noise from the release: `noise_dmsp_km` within the DMSP era and "
            "√(`noise_dmsp_km`² + `noise_viirs_km`²) across the join. Second-level units.\n\n"
            "**S2a. Rank correlation between σ and the error of the measured shift.**\n\n" + t1
            + "\n\n**S2b. Logistic models for the right direction.**\n\n" + t2
            + "\n\n**S2c. Out-of-sample test.**\n\n" + t3
            + "\n\n**S2d. Calibration of the σ model, ten equal-sized classes.**\n\n" + t4
            + "\n\n**S2e. Shift at which the σ model predicts 0.8 probability of the right direction.** "
              "The model failed its out-of-sample test (S2c, S2d), so these values must not be used as thresholds; "
              "they are shown for completeness.\n\n" + t5
            + "\n\n*95 percent block-bootstrap intervals in parentheses. Not computed for Canada and Japan: σ is not "
              "defined for a period entirely within VIIRS (Table S6). Source: docs/extended/detection_*.csv.*\n")


def s3():
    pos, pairs, chg = [], [], []
    for ds in ("census", "geostat", "inegi", "ibge"):
        d = read(E / f"{ds}_hybrid_position.csv")
        for _, r in d.iterrows():
            pos.append([NAME[ds], int(r.year), EST_A1[r.estimate], num(r.n, 0), ci(r.err_med, r.err_med_lo, r.err_med_hi)])
        for _, r in read(E / f"{ds}_hybrid_position_pairs.csv").iterrows():
            pairs.append([NAME[ds], int(r.year), PAIR_A1[r.pair], num(r.n, 0), ci(r.diff_med, r.diff_lo, r.diff_hi)])
        for _, r in read(E / f"{ds}_hybrid_change.csv").iterrows():
            chg.append([NAME[ds], period(r["from"], r.to), EST_A1[r.series]] + change_cols(r))
    for tag in ("canada", "japan"):
        for _, r in pd.read_csv(E / f"C_{tag}_hybrid.csv").iterrows():
            pos.append([NAME[tag], int(r.year), EST_A1[r.estimate], num(r.n, 0), num(r.err_med)])
    return ("## Table S3. Regional totals: three estimates of first-level centroids\n\n"
            "Light: the release's first-level light-weighted centroid. Hybrid, light: mean of the second-level "
            "light-weighted centroids, weighted by their census population. Hybrid, geometric: the same weights "
            "applied to the second-level geometric centres. A first-level unit is used if second-level units with a "
            "light-weighted centroid cover at least 99 percent of its census population. The hybrids use the "
            "benchmark's own population counts and are therefore not independent estimates.\n\n"
            "**S3a. Location: median distance to the census centroid (km).**\n\n"
            + md_table(["Benchmark", "Year", "Estimate", "Units", "Median error (km)"], pos)
            + "\n\n**S3b. Paired difference in location error per unit (km; positive: the second estimate is closer).**\n\n"
            + md_table(["Benchmark", "Year", "Difference", "Units", "Median difference (km)"], pairs)
            + "\n\n**S3c. Change.**\n\n"
            + md_table(["Benchmark", "Period", "Estimate", "Units", "True shift, median (km)", "Error |e − t|, median (km)",
                        "Right direction", "Rank correlation"], chg)
            + "\n\n*95 percent block-bootstrap intervals in parentheses; for Canada and Japan the hybrid is reported "
              "descriptively, without intervals, and its change was not computed. Brazil has one census year. "
              "Sources: docs/extended/<benchmark>_hybrid_*.csv, C_canada_hybrid.csv, C_japan_hybrid.csv.*\n")


def s4():
    sel, pos, pairs, chg = [], [], [], []
    miss_cols = ["missing_light", "missing_ghspop", "missing_ghsbuilts", "missing_worldpop", "missing_geometric"]
    for ds in ("census", "geostat", "inegi", "ibge"):
        for _, r in pd.read_csv(E / f"{ds}_benchmark_selection.csv").iterrows():
            sel.append([NAME[ds], LEVEL[r.level], int(r.year), num(r.n_facit, 0), num(r.n_used, 0), num(r.n_dropped, 0)]
                       + [num(r[c], 0) for c in miss_cols])
        for _, r in pd.read_csv(E / f"{ds}_benchmark_position.csv").iterrows():
            pos.append([NAME[ds], LEVEL[r.level], int(r.year), SOURCE[r.source], num(r.n, 0),
                        ci(r.err_med, r.err_med_lo, r.err_med_hi)])
        for _, r in pd.read_csv(E / f"{ds}_benchmark_pairs.csv").iterrows():
            pairs.append([NAME[ds], LEVEL[r.level], int(r.year), SOURCE[r.source], num(r.n, 0),
                          ci(r.closer_share, r.closer_lo, r.closer_hi), ci(r.diff_med, r.diff_lo, r.diff_hi)])
        for _, r in read(E / f"{ds}_benchmark_change.csv").iterrows():
            chg.append([NAME[ds], LEVEL[r.level], period(r["from"], r.to), SOURCE[r.source]] + change_cols(r))
    for tag in ("canada", "japan"):
        for _, r in pd.read_csv(E / f"C_{tag}_position.csv").iterrows():
            pos.append([NAME[tag], LEVEL[r.level], int(r.year), SOURCE[r.source], num(r.n, 0),
                        ci(r.err_med, r.err_med_lo, r.err_med_hi)])
        for _, r in pd.read_csv(E / f"C_{tag}_benchmark_pairs.csv").iterrows():
            pairs.append([NAME[tag], LEVEL[r.level], int(r.year), SOURCE[r.source], num(r.n, 0),
                          ci(r.closer_share, r.closer_lo, r.closer_hi), ci(r.diff_med, r.diff_lo, r.diff_hi)])
    return ("## Table S4. Comparison sources\n\n"
            "GHS-POP and WorldPop contain census counts and are not independent of the benchmarks; WorldPop also uses "
            "nighttime lights. Each source is taken at its own epoch nearest the census year (Table S6). Location "
            "measures are computed on the units with a value from every source.\n\n"
            "**S4a. Units used and lost.** Units with census data (complete, at least ten census points), units with a "
            "value from every source, and units lacking a value, by source.\n\n"
            + md_table(["Benchmark", "Level", "Year", "Census units", "Used", "Lost", "No light", "No GHS-POP",
                        "No built-up", "No WorldPop", "No geometric"], sel)
            + "\n\n*Canada and Japan: the confirmatory driver reports the units with a value per source (S4b); the "
              "number of units is the same for all sources in both benchmarks.*\n\n"
            + "**S4b. Location: median distance to the census centroid (km).**\n\n"
            + md_table(["Benchmark", "Level", "Year", "Source", "Units", "Median error (km)"], pos)
            + "\n\n**S4c. Paired comparison with light, per unit.** Share of units in which the source is closer to the "
              "census centroid than light, and median of light's error minus the source's error (km; positive: the "
              "source is closer).\n\n"
            + md_table(["Benchmark", "Level", "Year", "Source", "Units", "Share closer than light", "Median difference (km)"], pairs)
            + "\n\n**S4d. Change.** The sources' epochs can differ from the census years by up to a year.\n\n"
            + md_table(["Benchmark", "Level", "Period", "Source", "Units", "True shift, median (km)",
                        "Error |e − t|, median (km)", "Right direction", "Rank correlation"], chg)
            + "\n\n*95 percent block-bootstrap intervals in parentheses. Change was not computed for Brazil (one census "
              "year) or for the comparison sources in Canada and Japan. Sources: docs/extended/<benchmark>_benchmark_*.csv, "
              "C_canada_*.csv, C_japan_*.csv.*\n")


def s5():
    main, chg, outc = [], [], []
    for tag in ("canada", "japan"):
        m = pd.read_csv(E / f"C_{tag}_main.csv").iloc[0]
        y1, y2 = re.findall(r"\d{4}", m.measure)
        main.append([NAME[tag], f"{y1}–{y2}", num(m.n, 0), num(m.n_blocks, 0), ci(m.dir_right, m.dir_right_lo, m.dir_right_hi),
                     num(m.ci_width), ci(m.perm_null_mean, m.perm_null_lo, m.perm_null_hi),
                     num(m.perm_p_upper, 3), num(m.perm_p_lower, 3)])
        row = rule_outcome(m.dir_right, m.dir_right_lo, m.dir_right_hi, m.n)
        outc.append((NAME[tag], row))
        for _, r in pd.read_csv(E / f"C_{tag}_change.csv").iterrows():
            chg.append([NAME[tag], LEVEL[r.level], period(r["from"], r.to), num(r.n, 0), num(r.true_med),
                        num(r.err_med), num(r.ratio), ci(r.dir_right, r.dir_right_lo, r.dir_right_hi)])
    rule = [[str(i), cond + (" *(added after the Japanese result)*" if late else ""), concl,
             ", ".join(n for n, r in outc if r == i) or "–"] for i, cond, concl, late in RULE]
    return ("## Table S5. Confirmatory tests: Canada and Japan\n\n"
            "Primary measure: share of second-level units in which the measured shift points in the direction of the "
            "true shift. Blocks for the bootstrap: provinces and territories (Canada), prefectures (Japan). The "
            "permutation test shuffles measured shifts among units within the same block.\n\n"
            "**S5a. Primary measure.**\n\n"
            + md_table(["Benchmark", "Period", "Units", "Blocks", "Right direction", "Interval width",
                        "Permutation null", "p (upper tail)", "p (lower tail)"], main)
            + "\n\n**S5b. Interpretation rule, fixed in advance (plan §7), and the outcome.**\n\n"
            + md_table(["Row", "Outcome", "Conclusion", "Benchmark in this row"], rule)
            + "\n\n**S5c. Secondary measures of change.**\n\n"
            + md_table(["Benchmark", "Level", "Period", "Units", "True shift, median (km)", "Error |e − t|, median (km)",
                        "Ratio", "Right direction"], chg)
            + "\n\n*95 percent block-bootstrap intervals in parentheses. Row 3 for Canada triggered the second benchmark; "
              "under the stopping rule no third benchmark was added. Sources: docs/extended/C_canada_*.csv, "
              "C_japan_*.csv; the rule in docs/ANALYSIS_PLAN.md §7 and §10.*\n")


def s6():
    rows = parse_deviations(PLAN.read_text())
    check_deviations(rows, DEVIATIONS_EN)
    return ("## Table S6. Deviations from the analysis plan\n\n"
            + md_table(["Date", "Plan section", "Change", "Reason"], [list(r) for r in DEVIATIONS_EN])
            + "\n\n*Translated from the log in docs/ANALYSIS_PLAN.md §10, which is authoritative. Section numbers refer "
              "to the plan: §2 frozen measures, §6 comparison sources (B), §7 confirmatory tests (C).*\n")


def main():
    parts = ["# Supplementary Tables S1–S6\n",
             "*Generated by tools/write_supplement_tables.py from docs/extended/ and docs/ANALYSIS_PLAN.md. "
             "Do not edit by hand. Decimal points; e is the light-measured shift, t the true shift of the census "
             "centroid; right direction is the share of units with a positive cosine between e and t.*\n",
             s1(), s2(), s3(), s4(), s5(), s6()]
    OUT.write_text("\n".join(parts))
    print("skrivet:", OUT)


if __name__ == "__main__":
    main()
