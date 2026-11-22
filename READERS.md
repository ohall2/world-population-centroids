# Guide for readers of the manuscript

This guide is for readers of the draft *Can Nighttime Lights Track Where People Move? Census Evidence
on Light-Weighted Population Centroids, 1992–2025* (Ola Hall, Lund University). It says where
everything is, what can be rerun and what cannot. The repository is private and the draft is not for
circulation.

Everything refers to the branch **`utokade-kontroller`**. The branches `v2` and `master` are older.

**Language.** The manuscript, the supplement, the data set's README and this guide are in English. The
working documentation and the code comments are in Swedish. Table columns and file names are in English.

## 1. Where things are

| What | Where |
| --- | --- |
| Reading copy: manuscript, figures, supplement and supplementary tables in one PDF | [`docs/manuscript/manuscript_review.pdf`](docs/manuscript/manuscript_review.pdf) |
| Manuscript source, supplement, supplementary tables | [`docs/manuscript/`](docs/manuscript/) |
| Figures and the scripts that draw them | [`docs/figures/`](docs/figures/), [`tools/figures/`](tools/figures/) |
| The centroid series (data set v2.0.0) | Zenodo draft, reserved DOI 10.5281/zenodo.24336701, not yet public. The author sends a private preview link. |
| Census benchmarks and comparison-source centroids per unit | [`replication/`](replication/) |
| Raw data (about 40 GB, all from public sources) | Not redistributed. Sources, file names and download commands: [`docs/DATA.md`](docs/DATA.md) |
| Analysis plan, locked before the extended analyses were run (commit `dbbd88e`), with a log of deviations | [`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md) |
| Reviews of the manuscript and of the analyses, and the audit of every number in the manuscript | [`docs/review/`](docs/review/) |
| Method, and the decisions behind the series step by step | [`docs/METHODS.md`](docs/METHODS.md), [`docs/HARMONIZATION.md`](docs/HARMONIZATION.md) |
| Literature with verification status | [`docs/LITERATURE.md`](docs/LITERATURE.md) |

## 2. From the manuscript to the scripts

| In the manuscript | Result files | Script | Write-up |
| --- | --- | --- | --- |
| Table 1, census benchmarks | `replication/data/validation*/` | `tools/download_validation.py`, `tools/download_estat_jp.py`, `tools/validation_centroids.py` | `docs/DATA.md` |
| Location error (Table 2, Figure 1) | `docs/validation/*_position.csv`, `docs/figures/data/unit_errors.csv.gz` | `tools/validate_series.py`, `tools/figures/unit_errors.py` | `docs/HARMONIZATION.md` |
| Comparison sources (Table 2, Figure 6, Table S4) | `docs/extended/*_benchmark_*.csv` | `tools/compare_sources.py`, `tools/compare_benchmarks.py` | `docs/extended/B_benchmarks.md` |
| National bias and its decomposition (Figures 2–3) | `docs/check_identity/`, `docs/robustness/robust_national_correction.csv` | `tools/check_identity.py`, `tools/robustness_checks.py` | `docs/HARMONIZATION.md` |
| Detection of change (Table 3, Figures 4–5) | `docs/change_detection/`, `docs/robustness/`, `docs/figures/data/unit_shifts.csv.gz` | `tools/change_detection.py`, `tools/robustness_checks.py`, `tools/figures/unit_shifts.py` | `docs/HARMONIZATION.md` |
| Three-year means and trend estimates (Table S1) | `docs/extended/*_trend_change.csv` | `tools/trend_estimate.py` | `docs/extended/A2_trend.md` |
| Predicting which units can be trusted (Table S2) | `docs/extended/detection_*.csv` | `tools/detection_model.py` | `docs/extended/A3_detection.md` |
| Change within the VIIRS era, Canada and Japan (Figure 7, Table S5) | `docs/extended/C_canada_*.csv`, `docs/extended/C_japan_*.csv` | `tools/canada_confirm.py` | `docs/extended/C_canada.md`, `docs/extended/C_japan.md` |
| Regional totals (Figure 7b, Table S3) | `docs/extended/*_hybrid_*.csv` | `tools/hybrid_estimate.py` | `docs/extended/A1_hybrid.md` |
| The sensor join (Supplement S2, Figure S1) | `docs/check_seam/` | `tools/check_seam_measures.py`, `tools/check_seam_summary.py`, `tools/figures/figS1_seam.py` | `docs/HARMONIZATION.md` |
| Supplementary tables S1–S6 | `docs/manuscript/supplement_tables.md` | `tools/write_supplement_tables.py` | |
| Every number in the text | `docs/manuscript/numbers.csv` | `tools/manuscript_numbers.py` | `docs/review/talrevision-2026-11-22.md` |

## 3. What can be rerun from this repository and the data set

With the three `.csv.gz` files of the data set and the tables in `replication/`, the tests of change
can be rerun without any raw data. This was checked on 22 November 2026 in a fresh copy of the
repository: the three analyses below reproduced the committed tables, with identical point estimates.
Bootstrap interval bounds from `change_detection.py` can differ in the third decimal, depending on which
series and periods are run together.

```bash
python -m venv ~/wpc-venv && source ~/wpc-venv/bin/activate
pip install numpy pandas scipy matplotlib

# data: benchmarks from the repository, the series from the Zenodo preview link
mkdir -p data/release/v2.0.0 && cp -R replication/data/. data/
# put wpc_v2_gadm41_adm0.csv.gz, _adm1.csv.gz and _adm2.csv.gz in data/release/v2.0.0/

# detection of change (Table 3); one benchmark shown, a few seconds
python tools/change_detection.py --dataset inegi --pairs 2010-2020 --release data/release/v2.0.0 --out /tmp/cd
# three-year means and trends (Table S1), about 1.5 minutes
python -m tools.trend_estimate --out /tmp/trend
# model per unit (Table S2), about 20 seconds
python -m tools.detection_model --out /tmp/det

# figures (the scripts overwrite the files in docs/figures/)
python tools/figures/unit_errors.py && python tools/figures/fig1_location_error.py
python tools/figures/unit_shifts.py && python tools/figures/fig4_measured_vs_true.py && python tools/figures/fig5_direction_by_shift.py
python tools/figures/fig3_national_bias.py && python tools/figures/fig6_sources.py
python tools/figures/fig7_change_overview.py && python tools/figures/figS1_seam.py
```

Figures 3, 6, 7 and S1 are drawn from result tables in `docs/`; Figures 1, 4 and 5 are recomputed
from the series and the benchmarks.

Compare the output with the files of the same name in `docs/change_detection/` and `docs/extended/`.

## 4. What needs the raw data

The series itself, the comparison sources, the regional totals, the Canada and Japan tests, the
decomposition of the national bias, Figure 2 and the test suite read rasters or GADM boundaries. They need the full requirements
(`pip install -r requirements.txt`) and the downloads in `docs/DATA.md`, about 40 GB. The pipeline is
described in `README.md` under *Snabbstart*.

GADM 4.1 boundaries are free for academic use but may not be redistributed, so they are not in the
data set or in this repository. DMSP and VIIRS are downloaded through Google Earth Engine, which needs
an account; GADM and the Statistics Canada files are downloaded by hand.

## 5. Known gaps

- `docs/METHODS.md` section 9 describes the main series as of 17 November and has not been updated
  for the extended analyses; the manuscript and `docs/extended/` take precedence.
- The reference L. Li et al. (2026) has been checked against a summary only, not the full text.
- The reference list carries a working note on what is not yet checked: authors' initials in the
  older entries, and data citations that are still to be added.
- The manuscript is not yet anonymised and is about 10,900 words against a limit of 11,000.
- Questions and errors: ola.hall@keg.lu.se, or an issue in this repository.
