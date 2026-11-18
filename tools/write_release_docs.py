#!/usr/bin/env python3
"""
Skriv validering, förändringstest och README för den publicerade datamängden.

Läser den metadata.json som tools/build_release.py skrev (filer, kontrollsummor,
skarvning) och ersätter avsnitten om validering och begränsningar med de
slutliga resultaten:

  docs/figures/data/unit_errors.csv.gz   lägesfel per enhet (tools/figures/unit_errors.py)
  docs/change_detection/*_change_detection.csv   test av förändring (tools/change_detection.py)

USA och Puerto Rico redovisas var för sig (HARMONIZATION.md §17). README skrivs
på engelska, eftersom datamängden är avsedd för internationella användare.
Bara metadata.json och README.md i utgåvans mapp skrivs om; datafilerna rörs inte.

Exempel
  python tools/write_release_docs.py                      # data/release/v<version>
  python tools/write_release_docs.py --release data/release/v2.0.0
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc import __version__  # noqa: E402

BENCH = {"census": "US Census Bureau, block-group centres of population",
         "geostat": "Eurostat GEOSTAT 1-km population grid",
         "inegi": "INEGI locality register (ITER)",
         "ibge": "IBGE Grade Estatística 1-km grid (Censo 2022)"}
LEVEL = {0: "country", 1: "region (GADM level 1)", 2: "municipality (GADM level 2)"}


def group_name(ds: str, gid0: str) -> str:
    if ds == "census":
        return "United States" if gid0 == "USA" else "Puerto Rico"
    return {"geostat": "Europe", "inegi": "Mexico", "ibge": "Brazil"}[ds]


def position_table(unit_errors: Path) -> pd.DataFrame:
    d = pd.read_csv(unit_errors)
    d["benchmark"] = [group_name(ds, g) for ds, g in zip(d.dataset, d.gid0)]
    t = (d.groupby(["benchmark", "dataset", "level", "year"]).err_km
         .agg(units="size", median_km="median", p90_km=lambda s: s.quantile(0.9)).reset_index())
    t[["median_km", "p90_km"]] = t[["median_km", "p90_km"]].round(2)
    order = {"United States": 0, "Puerto Rico": 1, "Europe": 2, "Mexico": 3, "Brazil": 4}
    return t.sort_values(["level", "benchmark", "year"], key=lambda s: s.map(order) if s.name == "benchmark" else s)


def r2(x: float, nd: int = 2) -> float:
    """Avrundning halvt uppåt (CSV-värdena har fyra decimaler; Pythons round ger 0,585 -> 0,58)."""
    q = 10 ** nd
    return float(np.floor(abs(x) * q + 0.5 + 1e-9) / q * np.sign(x))


def change_table(cd_dir: Path) -> pd.DataFrame:
    rows = []
    for ds, name in (("census", "United States"), ("inegi", "Mexico"), ("geostat", "Europe")):
        c = pd.read_csv(cd_dir / f"{ds}_change_detection.csv")
        c = c[c.series == "release"].rename(columns={"from": "y0", "to": "y1"})
        if ds == "census":
            c = c[c.level == 2]  # nivå 1 domineras av Puerto Rico (HARMONIZATION §17)
        for r in c.itertuples():
            rows.append({"benchmark": name, "level": int(r.level), "period": f"{r.y0}–{r.y1}",
                         "units": int(r.n), "true_shift_median_km": r2(r.true_med),
                         "right_direction": r2(r.dir_right),
                         "right_direction_95ci": [r2(r.dir_right_lo), r2(r.dir_right_hi)],
                         "rank_correlation": r2(r.rho), "slope": r2(r.slope),
                         "share_true_shift_exceeds_error": r2(r.snr_gt1)})
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--release", type=Path, default=Path("data/release") / f"v{__version__}")
    p.add_argument("--unit-errors", type=Path, default=Path("docs/figures/data/unit_errors.csv.gz"))
    p.add_argument("--change", type=Path, default=Path("docs/change_detection"))
    a = p.parse_args(argv)

    mp = a.release / "metadata.json"
    meta = json.loads(mp.read_text())
    pos = position_table(a.unit_errors)
    chg = change_table(a.change)
    meta.pop("validation_position_error", None)
    meta["licence"] = {"data": "CC BY-NC 4.0 (follows the MYC25Q4 flare-site catalogue used for masking)",
                       "code": "MIT"}
    meta["creator"] = {"name": "Ola Hall", "affiliation": "Lund University",
                       "orcid": "https://orcid.org/0000-0002-9231-4028"}
    meta["validation"] = {
        "benchmarks": BENCH,
        "note": ("Great-circle distance between the light-weighted and the census Barmore centroid per unit. "
                 "A country is used if >= 90 % of its finest GADM units have census population; a unit needs >= 10 "
                 "census points. Level 0 is measured with light north of 65 N masked, levels 1-2 without. "
                 "The United States and Puerto Rico are reported separately."),
        "position_error_km": pos.to_dict(orient="records"),
        "change_detection": chg.to_dict(orient="records"),
        "change_detection_note": ("Per unit, the measured shift (light) is compared with the true shift (census). "
                                  "right_direction: share with positive cosine (0.5 by chance); 95 % block-bootstrap "
                                  "intervals. U.S. level 1 omitted (dominated by Puerto Rico's municipios)."),
    }
    mp.write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n")

    def fmt_pos(lv):
        s = pos[pos.level == lv]
        return "\n".join(f"| {r.benchmark} | {int(r.year)} | {int(r.units)} | {r.median_km:g} | {r.p90_km:g} |"
                         for r in s.itertuples())
    ctab = "\n".join(
        f"| {r.benchmark} | {r.level} | {r.period} | {r.units} | {r.true_shift_median_km:g} | "
        f"{r.right_direction:.2f} ({r.right_direction_95ci[0]:.2f}–{r.right_direction_95ci[1]:.2f}) | "
        f"{r.rank_correlation:.2f} | {r.share_true_shift_exceeds_error * 100:.0f}% |" for r in chg.itertuples())
    files = "\n".join(f"- `{n}`" for lv in meta["files"].values() for n in lv["sha256"])
    readme = f"""# World Population Centroids v2 ({__version__})

Light-weighted population centroids for the countries (GADM 4.1 level 0), regions (level 1) and
municipalities (level 2) of the world, annually 1992–2025. DMSP-OLS (1992–2012, intercalibrated and
co-registered to VIIRS) is spliced per unit onto VIIRS VNL v2 (2013–2025) at the overlap year 2013.
Code and full documentation: https://github.com/ohall2/world-population-centroids (docs/METHODS.md,
docs/HARMONIZATION.md).

Creator: Ola Hall, Lund University (ORCID 0000-0002-9231-4028).
Licence: data CC BY-NC 4.0 (follows the VIIRS Nightfire flare-site catalogue MYC25Q4 used for masking);
code MIT. GADM boundaries and identifiers are free for academic and other non-commercial use.

## What the data are good for

- **Location in a given year.** Within municipalities the centroid lies about 1 km from the census
  centroid (median 0.7–1.5 km in four census systems). Use it for distances, accessibility and
  gravity models.
- **Not change over time at levels 1–2.** Over 10–20 years the measured shift of a municipal centroid
  points the right way in 51–66 % of units (50 % by chance). Direction becomes reliable only for true
  shifts of several kilometres.
- **National centroids are biased** towards regions with high light per resident: 65–135 km in the
  United States, Mexico and Brazil.

## Files

{files}

One row per unit and year:

| Column | Meaning |
| --- | --- |
| unit_id, GID_*, COUNTRY, NAME_* | GADM 4.1 identifiers and names |
| year, sensor | 1992–2025; DMSP (1992–2012) or VIIRS (2013–2025) |
| bm_lon, bm_lat | Barmore centroid (degrees, WGS84) – main |
| aa_lon, aa_lat | Aboufadel–Austin centroid |
| bm_flag, aa_flag | 0 inside its unit, 1 outside all units, 2 in another unit of the same level, -1 no centroid |
| splice_ok | 1 if the unit had a centroid from both sensors in 2013; 0 means its DMSP years are not comparable with its VIIRS years |
| splice_shift_km | how far the DMSP years were shifted (Barmore); 0 for VIIRS years |
| noise_dmsp_km, noise_viirs_km | the unit's median year-to-year shift in each era – a noise level, not a confidence interval |
| n_cells, area_km2 | grid cells covered and area |
| max_lat | northernmost latitude of light used (65 at level 0, empty at levels 1–2) |

The GeoPackages hold the same rows with the Barmore centroid as point geometry (EPSG:4326).

## Location error against censuses (Barmore, km)

Level 0 (countries; light north of 65° N masked):

| Benchmark | Year | Units | Median | 90th percentile |
| --- | --- | --- | --- | --- |
{fmt_pos(0)}

Level 1:

| Benchmark | Year | Units | Median | 90th percentile |
| --- | --- | --- | --- | --- |
{fmt_pos(1)}

Level 2:

| Benchmark | Year | Units | Median | 90th percentile |
| --- | --- | --- | --- | --- |
{fmt_pos(2)}

## Detection of change

| Benchmark | Level | Period | Units | True shift, median (km) | Right direction (95 % CI) | Rank correlation | True shift > error |
| --- | --- | --- | --- | --- | --- | --- | --- |
{ctab}

## Limitations

- Change entirely within the VIIRS era (2013–2025) has not been tested against censuses.
- No census benchmark outside the Americas and Europe.
- National centroids exclude light north of 65° N (as in v1); levels 1–2 use all light. Units
  entirely north of 65° N have no national centroid.
- All VIIRS-based series, this one included, show a larger change from 2016 to 2017 (about 1.4 times
  a normal VIIRS year); it is not corrected.
- Units with splice_ok = 0 have DMSP years that are not comparable with their VIIRS years.
- Identifiers follow GADM 4.1, including its own codes such as "?" for unknown units.

## Benchmarks

{chr(10).join(f"- {k}: {v}" for k, v in BENCH.items())}
"""
    (a.release / "README.md").write_text(readme)
    print(f"skrev {mp} och {a.release / 'README.md'}")


if __name__ == "__main__":
    main()
