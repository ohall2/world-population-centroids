#!/usr/bin/env python3
"""
World Population Centroids v2 – fackelmasker för VIIRS-perioden (2012–2025).

Polygonerna från 2009 (tools/prepare_flares.py) räcker inte för VIIRS: fält
som vuxit fram efter 2009 saknas helt (Bakken, stora delar av Permian och
Västsibirien). Masken byggs därför av fackelplatser från VIIRS Nightfire:

  1. MYC25Q4 (Zhizhin m.fl. 2026, doi:10.5281/zenodo.23091822, CC BY-NC 4.0):
     27 229 upstream-platser (olja, gas, blandat) 2012–2025. Platser som
     aldrig har en detektion (n_det = 0 alla år) tas inte med.
  2. EOG:s årskataloger 2012–2024 (Excel, eogdata.mines.edu/global_flare_data):
     alla blad med facklor – upstream, midstream och downstream (gas: LNG och
     gasbehandling; olja: raffinaderier). De ger de anläggningar som MYC25Q4
     inte omfattar, t.ex. Hammerfest LNG.
  3. Polygonerna från 2009 (Elvidge m.fl. 2009), som för DMSP.

Varje plats får en cirkel med radien r km (sfärisk jord, R = 6371 km, 32
hörn). Masken är statisk: samma platser alla år, som för DMSP. En pixel som
cirkeln eller en 2009-polygon berör sätts till 0 (all_touched).

Cirklarna och 2009-polygonerna (fem av dem är ogiltiga och lagas med
make_valid) slås ihop till en yta, som delas i sina sammanhängande delar.

Utdata: en GeoPackage per radie, data/processed/flares_viirs_r<r>km.gpkg,
listan över platser (flares_viirs_sites.csv) och en JSON-rapport. Radien 0
ger bara 2009-polygonerna (för jämförelse).

Exempel
  python tools/prepare_flares_viirs.py --radii 0 1 2 3 5
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.progress import Progress, log  # noqa: E402

R_KM = 6371.0
N_VERTICES = 32


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_myc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    det = df.groupby("flare_id")["n_det"].sum()
    keep = det[det > 0].index
    sites = (df[df.flare_id.isin(keep)].groupby("flare_id")
             .agg(lat=("lat", "first"), lon=("lon", "first"),
                  kind=("upstream_class", "first")).reset_index())
    sites["kind"] = "upstream " + sites["kind"]
    sites["source"] = "MYC25Q4"
    return sites[["lat", "lon", "kind", "source"]]


def sheet_kind(name: str) -> str:
    """Bladnamn -> typ av anläggning (namnen varierar mellan åren)."""
    n = name.lower().replace("_", " ")
    for pat, kind in (("midstream", "midstream"), ("upstream", "upstream"),
                      ("oil downstream", "downstream oil"), ("downstream oil", "downstream oil"),
                      ("gas downstream", "downstream gas"), ("downstream gas", "downstream gas"),
                      ("downstream", "downstream")):
        if pat in n:
            return kind
    return "other"


def read_legacy(folder: Path) -> tuple[pd.DataFrame, list[dict]]:
    """Alla blad med fackelpositioner (kolumnerna Latitude och Longitude).

    Blad med landsummor saknar koordinater och hoppas över. Bladet "all flares"
    (2019) upprepar de andra bladen och hoppas också över.
    """
    rows, files = [], []
    for f in sorted(folder.glob("VIIRS_Global_flaring_*.xlsx")):
        years = re.search(r"_(\d{4}(?:-\d{4})?)_", f.name).group(1)
        sheets = pd.read_excel(f, sheet_name=None, engine="openpyxl")
        n, used = 0, []
        for sheet, d in sheets.items():
            if not {"Latitude", "Longitude"} <= set(d.columns) or sheet.lower() == "all flares":
                continue
            d = d[["Latitude", "Longitude"]].apply(pd.to_numeric, errors="coerce").dropna()
            rows.append(pd.DataFrame({"lat": d.Latitude, "lon": d.Longitude,
                                      "kind": sheet_kind(sheet), "source": f"EOG {years}"}))
            n += len(d)
            used.append(sheet)
        files.append({"file": f.name, "years": years, "sheets": used, "sites": n,
                      "sha256": sha256(f)})
    return pd.concat(rows, ignore_index=True), files


def circles(lon: np.ndarray, lat: np.ndarray, r_km: float) -> np.ndarray:
    """Sfäriska cirklar med radien r_km kring varje punkt, som polygoner."""
    d = r_km / R_KM
    brg = np.linspace(0, 2 * np.pi, N_VERTICES, endpoint=False)
    p1, l1 = np.radians(lat)[:, None], np.radians(lon)[:, None]
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(brg))
    l2 = l1 + np.arctan2(np.sin(brg) * np.sin(d) * np.cos(p1),
                         np.cos(d) - np.sin(p1) * np.sin(p2))
    ring = np.stack([np.degrees(l2), np.degrees(p2)], axis=-1)
    ring = np.concatenate([ring, ring[:, :1]], axis=1)
    return shapely.polygons(ring)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--myc", type=Path,
                   default=Path("data/raw/flares_viirs/annual_flare_site_records.csv"))
    p.add_argument("--legacy", type=Path, default=Path("data/raw/flares_viirs/legacy"))
    p.add_argument("--polygons-2009", type=Path, default=Path("data/processed/flares_dmsp.gpkg"))
    p.add_argument("--out-dir", type=Path, default=Path("data/processed"))
    p.add_argument("--radii", type=float, nargs="+", default=[2.0])
    args = p.parse_args(argv)

    with Progress("Läser fackelkataloger"):
        myc = read_myc(args.myc)
        legacy, files = read_legacy(args.legacy)
        old = gpd.read_file(args.polygons_2009, engine="pyogrio")
    sites = pd.concat([myc, legacy], ignore_index=True)
    sites = sites[sites.lat.between(-90, 90) & sites.lon.between(-180, 180)]
    n_all = len(sites)
    # Samma plats i flera kataloger: en cirkel räcker (avrundat till ca 10 m).
    sites = sites.assign(_k=sites.lat.round(4).astype(str) + "," + sites.lon.round(4).astype(str))
    sites = sites.drop_duplicates("_k").drop(columns="_k").reset_index(drop=True)
    log(f"MYC25Q4: {len(myc)} platser, EOG årskataloger: {len(legacy)} rader; "
        f"{len(sites)} unika platser av {n_all}. 2009-polygoner: {len(old)}")
    log("Per typ: " + ", ".join(f"{k} {v}" for k, v in sites.kind.value_counts().items()))

    report = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
              "myc25q4": {"file": str(args.myc), "sha256": sha256(args.myc), "sites": len(myc)},
              "legacy_files": files, "polygons_2009": {"file": str(args.polygons_2009), "n": len(old)},
              "unique_sites": len(sites),
              "sites_by_kind": {k: int(v) for k, v in sites.kind.value_counts().items()},
              "outputs": {}}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(args.out_dir / "flares_viirs_sites.csv", index=False, float_format="%.5f")
    old_geoms = shapely.make_valid(old.geometry.values)
    for r in args.radii:
        tag = f"{r:g}".replace(".", "p")
        out = args.out_dir / f"flares_viirs_r{tag}km.gpkg"
        geoms = old_geoms
        if r > 0:
            geoms = np.concatenate([geoms, circles(sites.lon.values, sites.lat.values, r)])
        with Progress(f"Slår ihop masken för {r:g} km"):
            parts = shapely.get_parts(shapely.union_all(geoms))
        gdf = gpd.GeoDataFrame({"radius_km": np.full(len(parts), r)}, geometry=parts, crs=4326)
        tmp = out.with_suffix(".part.gpkg")
        tmp.unlink(missing_ok=True)
        with Progress(f"Skriver {out.name} ({len(gdf)} polygoner)"):
            gdf.to_file(tmp, driver="GPKG", engine="pyogrio")
        tmp.replace(out)
        report["outputs"][out.name] = {"radius_km": r, "polygons": len(gdf), "sha256": sha256(out)}
    rep = args.out_dir / "flares_viirs_report.json"
    rep.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    log(f"Klart. Rapport: {rep}")


if __name__ == "__main__":
    main()
