#!/usr/bin/env python3
"""
World Population Centroids v2 – hämta valideringsdata (befolkning).

USA (U.S. Census Bureau, public domain): befolkning och befolkningens
mittpunkt per kvartersgrupp (block group), 2000, 2010 och 2020.
  https://www.census.gov/geographies/reference-files/time-series/geo/centers-population.html
  2020: cenpop2020/blkgrp/CenPop2020_Mean_BG.txt (en fil)
  2010: cenpop2010/blkgrp/CenPop2010_Mean_BG.txt (en fil)
  2000: cenpop2000/blkgrp/bg_<fips>_<st>.txt (en fil per delstat, utan rubrikrad)
  Kolumner: STATEFP, COUNTYFP, TRACTCE, BLKGRPCE, POPULATION, LATITUDE, LONGITUDE.

Europa (Eurostat GISCO, folkräkningar på rutnät 1 km², ETRS89-LAEA):
  2011: GEOSTAT-grid-POP-1K-2011-V2-0-1.zip
  2021: Eurostat_Census-GRID_2021_V3.zip
  https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/population-grids
  (JRC-rutnätet 2018 tas inte med: det är modellerat, inte en folkräkning.)

Filerna hamnar i data/raw/validation/<källa>/. Avbrutna nedladdningar
fortsätter där de slutade. En manifestfil (manifest.json) anger URL, storlek,
SHA-256 och tidpunkt för varje fil.

Kräver bara Pythons standardbibliotek.

Exempel
  python tools/download_validation.py
  python tools/download_validation.py --only census
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

CENSUS = "https://www2.census.gov/geo/docs/reference"
FILES = {
    "census": [
        ("2020", f"{CENSUS}/cenpop2020/blkgrp/CenPop2020_Mean_BG.txt"),
        ("2010", f"{CENSUS}/cenpop2010/blkgrp/CenPop2010_Mean_BG.txt"),
    ],
    "geostat": [
        ("2011", "https://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/GEOSTAT-grid-POP-1K-2011-V2-0-1.zip"),
        ("2021", "https://gisco-services.ec.europa.eu/census/2021/Eurostat_Census-GRID_2021_V3.zip"),
    ],
}
CENSUS_2000_INDEX = f"{CENSUS}/cenpop2000/blkgrp/"
TIMEOUT = 120
CHUNK = 1 << 20


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def census_2000_files() -> list[tuple[str, str]]:
    with urllib.request.urlopen(CENSUS_2000_INDEX, timeout=TIMEOUT) as r:
        html = r.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r'href="(bg_\d\d_[a-z]{2}\.txt)"', html)))
    if not names:
        sys.exit(f"Hittade inga filer i {CENSUS_2000_INDEX}")
    return [("2000", CENSUS_2000_INDEX + n) for n in names]


def fetch(url: str, dst: Path, retries: int = 8) -> None:
    part = dst.with_name(dst.name + ".part")
    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                total = r.headers.get("Content-Range", "").split("/")[-1] or r.headers.get("Content-Length")
                if have and r.status != 206:
                    have = 0
                    part.unlink(missing_ok=True)
                t0, last = time.time(), 0.0
                with open(part, "ab") as out:
                    while chunk := r.read(CHUNK):
                        out.write(chunk)
                        have += len(chunk)
                        if time.time() - last > 5:
                            last = time.time()
                            tot = f" av {int(total) / 1e6:.0f}" if total and total.isdigit() else ""
                            print(f"\r  {dst.name}: {have / 1e6:.0f}{tot} MB, {time.time() - t0:.0f} s   ",
                                  end="", flush=True)
            print()
            part.replace(dst)
            return
        except Exception as exc:  # noqa: BLE001
            if attempt == retries:
                raise
            wait = min(60, 2 ** attempt)
            log(f"  {dst.name}: {str(exc)[:100]} – fortsätter om {wait} s")
            time.sleep(wait)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(description="Hämta valideringsdata (befolkning).")
    p.add_argument("--out", type=Path, default=Path("data/raw/validation"))
    p.add_argument("--only", choices=sorted(FILES), default=None)
    args = p.parse_args(argv)
    t0 = time.time()
    for src, files in FILES.items():
        if args.only and src != args.only:
            continue
        if src == "census":
            files = files + census_2000_files()
        d = args.out / src
        d.mkdir(parents=True, exist_ok=True)
        man_path = d / "manifest.json"
        manifest = json.loads(man_path.read_text()) if man_path.exists() else {}
        for year, url in files:
            name = url.rsplit("/", 1)[1]
            dst = d / year / name
            dst.parent.mkdir(exist_ok=True)
            if dst.exists() and name in manifest:
                continue
            log(f"{src} {year}: {name}")
            fetch(url, dst)
            manifest[name] = {"year": year, "url": url, "bytes": dst.stat().st_size,
                              "sha256": sha256(dst),
                              "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
            man_path.write_text(json.dumps(manifest, indent=1) + "\n")
        log(f"{src}: {len(manifest)} filer i {d}")
    log(f"Klart på {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
