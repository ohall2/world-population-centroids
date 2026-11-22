# Replication data

Small derived tables that are not part of the published centroid series but are needed to rerun the
tests in the manuscript. They are copies of what the pipeline writes under `data/`, which is otherwise
kept outside git. Nothing here is raw data.

To use them, copy them into `data/` from the repository root:

```bash
mkdir -p data && cp -R replication/data/. data/
```

## Contents

| Folder | What | Written by |
| --- | --- | --- |
| `data/validation/census/{2000,2010,2020}` | Census population centroids per GADM 4.1 unit, United States and Puerto Rico | `tools/validation_centroids.py` |
| `data/validation/geostat/{2011,2021}` | Same, Europe (Eurostat GEOSTAT 1-km grid) | `tools/validation_centroids.py` |
| `data/validation/inegi/{2010,2020}` | Same, Mexico (INEGI locality register) | `tools/validation_centroids.py` |
| `data/validation/ibge/2022` | Same, Brazil (IBGE 1-km grid, Censo 2022) | `tools/validation_centroids.py` |
| `data/validation_canada/statcan/{2016,2021}` | Same, Canada (Statistics Canada Geographic Attribute File) | `tools/validation_centroids.py` |
| `data/validation_japan/estat_jp/{2015,2020}` | Same, Japan (e-Stat 500-m grid square statistics) | `tools/validation_centroids.py` |
| `data/results/gadm41/{ghsbuilts,ghspop,worldpop}` | Centroids of the benchmark units computed from GHS-BUILT-S, GHS-POP and WorldPop | `tools/compare_sources.py` |

Each benchmark folder holds `adm0.csv`, `adm1.csv`, `adm2.csv` and `summary.json`.

Benchmark columns: `unit_id` (GADM 4.1), `level`, `year`, `dataset`, `pop` (census population in the unit),
`n_points` (census points in the unit), `aa_lon`, `aa_lat` (Aboufadel–Austin centroid), `bm_lon`, `bm_lat`
(Barmore centroid), `missing_light_share`, `complete` (the unit enters the tests only if `True`).
`summary.json` records the source file, the national population total and coverage per level.

Comparison-source columns: `unit_id`, `year` (the source's epoch), `bm_lon`, `bm_lat`, `aa_lon`, `aa_lat`.

## Sources

How each source is downloaded, with file names, checksums and the licence notes we have recorded, is
in [`docs/DATA.md`](../docs/DATA.md).

- U.S. Census Bureau, block-group centres of population, 2000, 2010 and 2020.
- Eurostat, GEOSTAT 1-km population grid, 2011 (V2.0.1) and 2021 (V3).
- INEGI, Principales resultados por localidad (ITER), 2010 and 2020.
- IBGE, Grade Estatística 1 km, Censo Demográfico 2022.
- Statistics Canada, Geographic Attribute File, Census 2016 and 2021 (Open Government Licence – Canada).
- Statistics Bureau of Japan, grid square statistics of the Population Census 2015 and 2020. Source: e-Stat.
- JRC GHSL, GHS-BUILT-S and GHS-POP, release R2023A (CC BY 4.0).
- WorldPop, Global 2000–2020, unconstrained, 1 km (CC BY 4.0).

The tables here are aggregates computed by us from those sources. The terms for redistributing them
have not yet been reviewed source by source.
