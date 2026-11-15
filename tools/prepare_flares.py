#!/usr/bin/env python3
"""
World Population Centroids v2 – hämta och slå ihop NGDC:s fackelpolygoner för DMSP.

Källa
  Elvidge m.fl. (2009), "A Fifteen Year Record of Global Natural Gas Flaring
  Derived from Satellite Data", Energies 2:595–622. Polygonerna ritades för
  hand kring gasfacklor i DMSP-data, en fil per land eller region
  (version 2009-06-18). Samma polygoner användes i v1 (Hall m.fl. 2019).
  https://eogdata.mines.edu/interest/flare_docs/gas_flares_countries_shapefiles.html

Vad skriptet gör
  1. Hämtar de 65 arkiven (kräver inget konto) till --raw-dir.
  2. Packar upp dem och läser alla polygoner.
  3. Kontrollerar att koordinaterna är giltiga grader. Tre filer
     (Elfenbenskusten, Ghana, Mauretanien) har koordinater som inte är
     grader – troligen bildkoordinater – och utesluts. Alla tre gäller
     facklor till havs, som ligger utanför GADM:s landpolygoner och därför
     inte påverkar centroiderna.
  4. Skriver en GeoPackage (EPSG:4326) med en rad per polygon och
     källfilens namn, samt en JSON-rapport.

Exempel
  python tools/prepare_flares.py --raw-dir data/raw/flares \
      --out data/processed/flares_dmsp.gpkg
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import tarfile
import threading
import time
import urllib.request
from pathlib import Path

import geopandas as gpd
import pandas as pd

BASE_URL = ("https://eogdata.mines.edu/interest/flare_docs/gasflares_v2/"
            "country_vectors_20090618/")
NAMES = [
    "Algeria", "Angola", "Argentina", "Australia", "Azerbaijan", "Belgium",
    "Bolivia", "Brazil", "Brunei", "Cameroon", "Canada", "Chad", "Chile",
    "China", "Colombia", "Congo", "Cote_dIvoire", "Dem_Rep_Congo", "Denmark",
    "East_Timor", "Ecuador", "Egypt", "Eq_Guinea", "Gabon", "Germany", "Ghana",
    "India", "Indonesia", "Iran", "Iraq", "Kazakhstan", "Kuwait", "Libya",
    "Malaysia", "Mauritania", "Mexico", "Myanmar", "Netherlands", "Nigeria",
    "Norway", "Oman", "Peru", "Philippines", "PNG", "Qatar", "Romania",
    "Russia_KM", "Russia_not_KM", "Saudi_Arabia", "South_Africa", "Sudan",
    "Syria", "Thailand", "Trinidad", "Tunisia", "Turkmenistan", "UAE", "UK",
    "USA_Alaska", "USA_CONUS", "Uzbekistan", "Venezuela", "Vietnam", "Yemen",
    "New_Zealand",
]
SPECIAL = {"Venezuela": "Flares_Venezuela_1_point.tgz",
           "New_Zealand": "Flares_New_Zealand.tgz"}
EXCLUDED = {"Cote_dIvoire", "Ghana", "Mauritania"}  # ogiltiga koordinater, se ovan


def archive_name(name: str) -> str:
    return SPECIAL.get(name, f"Flares_{name}_1.tgz")


class Ticker:
    """Visar förfluten tid och en snurra så att det syns att skriptet arbetar."""

    def __init__(self, text: str):
        self.text, self.t0 = text, time.time()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self) -> None:
        spin = "|/-\\"
        i = 0
        while not self._stop.wait(0.25):
            sys.stdout.write(f"\r  {spin[i % 4]} {self.text}  {time.time() - self.t0:6.1f} s")
            sys.stdout.flush()
            i += 1

    def update(self, text: str) -> None:
        self.text = text

    def __enter__(self) -> "Ticker":
        self._thread.start()
        return self

    def __exit__(self, *exc) -> None:
        self._stop.set()
        self._thread.join()
        sys.stdout.write(f"\r  ✓ {self.text}  {time.time() - self.t0:6.1f} s\n")
        sys.stdout.flush()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download(raw_dir: Path) -> dict:
    raw_dir.mkdir(parents=True, exist_ok=True)
    checksums = {}
    with Ticker("Hämtar arkiv") as t:
        for i, name in enumerate(NAMES, start=1):
            fname = archive_name(name)
            target = raw_dir / fname
            t.update(f"Hämtar arkiv {i}/{len(NAMES)}: {fname}")
            if not target.exists():
                for attempt in range(1, 4):
                    try:
                        urllib.request.urlretrieve(BASE_URL + fname, target.with_suffix(".part"))
                        target.with_suffix(".part").replace(target)
                        break
                    except Exception:  # noqa: BLE001
                        if attempt == 3:
                            raise
                        time.sleep(2 * attempt)
            checksums[fname] = sha256(target)
        t.update(f"Hämtade {len(NAMES)} arkiv")
    return checksums


def read_all(raw_dir: Path) -> tuple[gpd.GeoDataFrame, list[dict]]:
    work = raw_dir / "extracted"
    work.mkdir(exist_ok=True)
    parts, log = [], []
    with Ticker("Läser polygoner") as t:
        for i, name in enumerate(NAMES, start=1):
            fname = archive_name(name)
            t.update(f"Läser {i}/{len(NAMES)}: {fname}")
            with tarfile.open(raw_dir / fname) as tar:
                members = [m for m in tar.getmembers() if m.isfile() and "/" not in m.name]
                tar.extractall(work, members=members, filter="data")
                shp = next(m.name for m in members if m.name.endswith(".shp"))
            g = gpd.read_file(work / shp, engine="pyogrio")
            minx, miny, maxx, maxy = g.total_bounds
            valid = -180 <= minx and maxx <= 180 and -90 <= miny and maxy <= 90
            entry = {"archive": fname, "shapefile": shp, "polygons": int(len(g)),
                     "crs": str(g.crs), "bounds": [float(minx), float(miny), float(maxx), float(maxy)],
                     "used": bool(valid and name not in EXCLUDED)}
            if name in EXCLUDED and valid:
                raise RuntimeError(f"{fname} har giltiga koordinater men står som utesluten – kontrollera")
            if not valid and name not in EXCLUDED:
                raise RuntimeError(f"{fname} har ogiltiga koordinater {g.total_bounds}")
            log.append(entry)
            if entry["used"]:
                parts.append(g.to_crs(4326).assign(source=fname)[["source", "geometry"]])
        t.update("Läste alla polygoner")
    flares = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=4326)
    return flares, log


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Hämta och slå ihop NGDC:s fackelpolygoner.")
    p.add_argument("--raw-dir", default="data/raw/flares", type=Path)
    p.add_argument("--out", default="data/processed/flares_dmsp.gpkg", type=Path)
    args = p.parse_args(argv)

    checksums = download(args.raw_dir)
    flares, log = read_all(args.raw_dir)
    if set(flares.geom_type) != {"Polygon"}:
        sys.exit(f"Oväntade geometrityper: {set(flares.geom_type)}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    flares.to_file(args.out, layer="flares", driver="GPKG", engine="pyogrio")
    n_used = sum(e["used"] for e in log)
    report = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "source": BASE_URL,
        "reference": "Elvidge et al. 2009, Energies 2:595-622",
        "archives": len(log), "archives_used": n_used,
        "polygons_used": int(len(flares)),
        "excluded": [e["archive"] for e in log if not e["used"]],
        "excluded_reason": "koordinaterna är inte grader (troligen bildkoordinater); "
                           "gäller facklor till havs utanför landpolygonerna",
        "files": log, "sha256": checksums,
    }
    args.out.with_suffix(".report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"Klart: {n_used} av {len(log)} filer, {len(flares)} polygoner -> {args.out}")


if __name__ == "__main__":
    main()
