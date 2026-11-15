#!/usr/bin/env python3
"""
World Population Centroids v2 – gör GADM 2.8 nivå 0–2 klar för reproduktionen av v1.

Bakgrund
  GADM:s egen server för version 2.8 finns inte längre. Geometrierna
  (gadm28_adm{0,1,2}.shp/.shx/.prj) hämtas från NOAA GML:s spegel, men den
  servern lämnar inte ut .dbf-filerna. adm0-tabellen finns i en Zenodo-post
  (doi:10.5281/zenodo.5103739, md5 090efe30d3cb99f04a4814966832da6e).
  För adm1 och adm2 byggs attributtabellen upp ur v1:s egna utdata: varje
  enhet finns antingen bland de placerade centroiderna (*_final_flagged.shp)
  eller bland de oplacerade (*_novalue_final.csv), och OBJECTID är
  radnumret (1, 2, …) i GADM-filen.

Kontroller (skriptet avbryter om någon misslyckas)
  1. Antalet poster i varje shapefil är 256 / 3 609 / 46 311.
  2. Attributen ur v1 täcker OBJECTID 1..N exakt en gång.
  3. ID-fälten (OBJECTID, ID_*, ISO) är identiska i alla 22 år och för
     båda metoderna i v1. Namnen tas från Barmore 1992 (UTF-8); avvikande
     namn i andra år och i Aboufadel–Austins Windows-1252-filer redovisas
     men stoppar inte, eftersom namnen inte påverkar centroiderna.
     Blanktecken i början och slutet av namn tas bort (GADM 2.8 har t.ex.
     "Sark\\r\\n"), och saknade namn blir tomma strängar.
  4. adm0: OBJECTID i GADM-tabellen är 1..256 i radordning, och ID, ISO och
     namn stämmer med v1.
  5. Varje v1-centroid som ligger i sin egen enhet (OUT_FLAG = 0) ligger
     inuti polygonen med samma OBJECTID. Det visar att radordningen stämmer.

Utdata
  En GeoPackage med lagren adm0, adm1 och adm2 (EPSG:4326) samt en
  JSON-rapport med kontrollernas resultat.

Exempel
  python tools/prepare_gadm28.py --gadm-dir data/raw/gadm28 \
      --v1-dir legacy/v1 --out data/processed/gadm28.gpkg
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio

log = logging.getLogger("prepare_gadm28")

EXPECTED_N = {0: 256, 1: 3609, 2: 46311}
KEYS = {
    0: ["OBJECTID", "ID_0", "ISO", "NAME_ENGLI"],
    1: ["OBJECTID", "ID_0", "ISO", "NAME_0", "ID_1", "NAME_1"],
    2: ["OBJECTID", "ID_0", "ISO", "NAME_0", "ID_1", "ID_2", "NAME_1", "NAME_2"],
}
METHODS = {"Barmore": "Barmore", "AboufadelAustin": "AboufadelA"}
YEARS = range(1992, 2014)


class CheckFailed(RuntimeError):
    pass


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise CheckFailed(msg)
    log.info("OK  %s", msg)


# --------------------------------------------------------------------------
# v1-attribut
# --------------------------------------------------------------------------

def v1_paths(v1_dir: Path, method: str, year: int, level: int) -> tuple[Path, Path]:
    tag = METHODS[method]
    stem = f"DMSP_{year}_gadm28_adm{level}_gcs_{tag}"
    return (v1_dir / method / f"{stem}_final_flagged.shp",
            v1_dir / method / f"{stem}_novalue_final.csv")


def read_v1(v1_dir: Path, method: str, year: int, level: int) -> tuple[gpd.GeoDataFrame, pd.DataFrame]:
    shp, csv = v1_paths(v1_dir, method, year, level)
    placed = gpd.read_file(shp, engine="pyogrio")
    # Barmore-filerna är UTF-8; Aboufadel–Austins CSV-filer för adm1/adm2 är
    # Windows-1252. Prova UTF-8 först.
    try:
        novalue = pd.read_csv(csv, sep=";", encoding="utf-8-sig", dtype=str,
                              keep_default_na=False)
    except UnicodeDecodeError:
        novalue = pd.read_csv(csv, sep=";", encoding="cp1252", dtype=str,
                              keep_default_na=False)
    return placed, novalue


def normalise(df: pd.DataFrame, level: int) -> pd.DataFrame:
    """Samma kolumner och typer för placerade och oplacerade enheter."""
    out = df[KEYS[level]].copy()
    for col in KEYS[level]:
        if col == "OBJECTID" or col.startswith("ID_"):
            out[col] = pd.to_numeric(out[col]).astype("int64")
        else:
            # Namnfälten kan ha avslutande radbyten i GADM 2.8 (t.ex. "Sark\r\n"),
            # men bara i den ena av v1:s två källor. Blanktecken i början och
            # slutet tas därför bort; ID-fälten påverkas inte.
            out[col] = out[col].fillna("").astype(str).str.strip()
    return out.sort_values("OBJECTID").reset_index(drop=True)


def v1_attributes(v1_dir: Path, level: int) -> tuple[pd.DataFrame, dict]:
    """Attributtabell för alla enheter, kontrollerad över år och metoder.

    Barmore 1992 (UTF-8) är referens. Mot den kontrolleras alla 22 år för
    båda metoderna:
      - ID-fälten (OBJECTID, ID_*, ISO) ska vara identiska, annars avbryts.
      - Avvikande namn räknas och redovisas men stoppar inte. Namnen påverkar
        inte centroiderna, och GADM 2.8 har namn med radbyten inuti
        (t.ex. "Mongomeyen\\r\\nMongomeyen") som läses olika mellan år.
    """
    n = EXPECTED_N[level]
    id_cols = [c for c in KEYS[level] if c == "OBJECTID" or c.startswith("ID_") or c == "ISO"]
    name_cols = [c for c in KEYS[level] if c not in id_cols]
    reference = None
    name_diffs = {"Barmore": 0, "AboufadelAustin": 0}
    for method in ("Barmore", "AboufadelAustin"):
        for year in YEARS:
            placed, novalue = read_v1(v1_dir, method, year, level)
            table = normalise(pd.concat(
                [pd.DataFrame(placed.drop(columns="geometry")), novalue],
                ignore_index=True), level)
            if reference is None:
                reference = table
                check(len(table) == n and table["OBJECTID"].tolist() == list(range(1, n + 1)),
                      f"adm{level}: v1-attributen täcker OBJECTID 1..{n} exakt en gång")
                continue
            if len(table) != n or not table[id_cols].equals(reference[id_cols]):
                raise CheckFailed(f"adm{level}: ID-fälten i v1 skiljer sig för {method} {year}")
            diff = int((table[name_cols] != reference[name_cols]).any(axis=1).sum())
            name_diffs[method] = max(name_diffs[method], diff)
    log.info("OK  adm%d: ID-fälten är identiska i 22 år x 2 metoder", level)
    log.info("    adm%d: enheter med avvikande namn mot Barmore 1992, högst per år: "
             "Barmore %d, Aboufadel–Austin %d", level,
             name_diffs["Barmore"], name_diffs["AboufadelAustin"])
    return reference, {"name_reference": "Barmore 1992",
                       "max_name_differences_per_year": name_diffs}


# --------------------------------------------------------------------------
# Kontroll av radordning med v1-centroider
# --------------------------------------------------------------------------

def check_order(polys: gpd.GeoDataFrame, v1_dir: Path, level: int) -> dict:
    placed, _ = read_v1(v1_dir, "Barmore", 2013, level)
    # v1-filernas .prj anger en sfär (radie 6 371 000 m) i stället för WGS84.
    # Koordinaterna är grader i båda fallen, så punkterna kan jämföras direkt.
    placed = placed.set_crs(polys.crs, allow_override=True)
    inside = placed[placed["OUT_FLAG"] == 0]
    geoms = polys.set_index("OBJECTID").geometry
    own = geoms.loc[inside["OBJECTID"].to_numpy()].reset_index(drop=True)
    hit = own.covers(inside.geometry.reset_index(drop=True)).to_numpy()
    n_miss = int((~hit).sum())
    result = {"checked": int(len(inside)), "outside_own_polygon": n_miss,
              "missing_objectids": inside.loc[~hit, "OBJECTID"].astype(int).tolist()[:20]}
    check(n_miss == 0,
          f"adm{level}: alla {len(inside)} v1-centroider med OUT_FLAG=0 ligger i "
          f"polygonen med samma OBJECTID")
    return result


# --------------------------------------------------------------------------
# Huvudflöde
# --------------------------------------------------------------------------

def prepare(gadm_dir: Path, v1_dir: Path, out: Path) -> dict:
    report: dict = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    "gadm_dir": str(gadm_dir), "v1_dir": str(v1_dir), "levels": {}}
    tmp = out.with_suffix(".part.gpkg")
    tmp.unlink(missing_ok=True)

    for level in (0, 1, 2):
        shp = gadm_dir / f"gadm28_adm{level}.shp"
        n = pyogrio.read_info(shp)["features"]
        check(n == EXPECTED_N[level], f"adm{level}: {n} poster i {shp.name}")

        if level == 0:
            polys = gpd.read_file(shp, engine="pyogrio", columns=KEYS[0])
            polys = polys.astype({"OBJECTID": "int64", "ID_0": "int64"})
            for col in ("ISO", "NAME_ENGLI"):
                polys[col] = polys[col].astype(str).str.strip()
            check(polys["OBJECTID"].tolist() == list(range(1, n + 1)),
                  "adm0: OBJECTID i GADM-tabellen är 1..256 i radordning")
            attrs, attr_info = v1_attributes(v1_dir, 0)
            check(pd.DataFrame(polys[KEYS[0]]).reset_index(drop=True).equals(attrs),
                  "adm0: GADM-tabellen och v1-attributen är identiska")
        else:
            polys = gpd.read_file(shp, engine="pyogrio")
            polys = polys[["geometry"]].copy()
            polys.insert(0, "OBJECTID", np.arange(1, n + 1, dtype="int64"))
            attrs, attr_info = v1_attributes(v1_dir, level)
            polys = polys.merge(attrs, on="OBJECTID", how="left", validate="one_to_one")
            polys = gpd.GeoDataFrame(polys[KEYS[level] + ["geometry"]],
                                     geometry="geometry", crs=polys.crs)

        check(str(polys.crs.to_epsg()) == "4326", f"adm{level}: koordinatsystem EPSG:4326")
        order = check_order(polys, v1_dir, level)
        n_invalid = int((~polys.geometry.is_valid).sum())
        n_empty = int(polys.geometry.is_empty.sum())
        log.info("    adm%d: %d ogiltiga och %d tomma geometrier (lämnas orörda)",
                 level, n_invalid, n_empty)

        polys.to_file(tmp, layer=f"adm{level}", driver="GPKG", engine="pyogrio")
        report["levels"][f"adm{level}"] = {
            "features": n, "fields": KEYS[level], "order_check": order,
            "attributes": attr_info,
            "invalid_geometries": n_invalid, "empty_geometries": n_empty}
        log.info("    adm%d skrivet till GeoPackage", level)

    tmp.replace(out)
    return report


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Gör GADM 2.8 nivå 0–2 klar för reproduktionen av v1.")
    p.add_argument("--gadm-dir", default="data/raw/gadm28", type=Path)
    p.add_argument("--v1-dir", default="legacy/v1", type=Path,
                   help="Klon av ohall2/world-population-centroids (v1).")
    p.add_argument("--out", default="data/processed/gadm28.gpkg", type=Path)
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = prepare(args.gadm_dir, args.v1_dir, args.out)
    except CheckFailed as exc:
        sys.exit(f"KONTROLL MISSLYCKADES: {exc}")
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    log.info("Klart: %s och %s", args.out, report_path)


if __name__ == "__main__":
    main()
