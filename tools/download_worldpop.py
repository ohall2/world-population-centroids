#!/usr/bin/env python3
"""Hämta WorldPop Global (unconstrained, 1 km) direkt från WorldPop.

Jämförelsekälla för analys B (docs/ANALYSIS_PLAN.md §6). Den globala mosaiken
ppp_<år>_1km_Aggregated (befolkning per cell, 1 km, EPSG:4326, unconstrained)
laddas ner i sitt eget rutnät, ingen pixel resamplas. För varje år hämtas
GeoTIFF:en (återupptagbart) och en JSON-fil skrivs bredvid med ursprung,
rutnät, datatyp och SHA-256 (samma mönster som tools/download_ghsl.py;
pipelinen hittar filerna via JSON-filernas "year"). Finns .tif och .json redan
hoppas året över.

WorldPop använder nattljus som en av förklaringsvariablerna och är därför inte
oberoende av ljuset; redovisas skilt (se docs/DATA.md och planen §6).

Källa och licens: WorldPop (www.worldpop.org), CC BY 4.0.

Exempel
  python tools/download_worldpop.py --years 2000 2010 2011 2020
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import rasterio

BASE = "https://data.worldpop.org/GIS/Population/Global_2000_2020"


def stem(year: int) -> str:
    return f"ppp_{year}_1km_Aggregated"


def url(year: int) -> str:
    return f"{BASE}/{year}/0_Mosaicked/{stem(year)}.tif"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(year: int, out: Path) -> None:
    tif = out / f"{stem(year)}.tif"
    js = out / f"{stem(year)}.json"
    if tif.exists() and js.exists():
        print(f"{year}: finns redan, hoppar över")
        return
    print(f"{year}: hämtar {url(year)}")
    subprocess.run(["curl", "-fSL", "-C", "-", "-o", str(tif), url(year)], check=True)
    with rasterio.open(tif) as r:
        meta = {"source": "WorldPop (www.worldpop.org)",
                "product": "WorldPop Global 2000-2020, unconstrained, 1 km, population per cell",
                "version": "Global_2000_2020 (unconstrained, 1km Aggregated)", "license": "CC BY 4.0",
                "url": url(year), "year": year,
                "crs": str(r.crs), "transform": list(r.transform)[:6],
                "width": r.width, "height": r.height, "dtype": r.dtypes[0],
                "nodata": r.nodata, "sha256": sha256(tif)}
    js.write_text(json.dumps(meta, indent=2))
    print(f"{year}: klar ({tif.name}, {tif.stat().st_size/1e9:.2f} GB)")


def main(argv=None):
    p = argparse.ArgumentParser(description="Hämta WorldPop Global unconstrained 1 km.")
    p.add_argument("--years", type=int, nargs="+", default=[2000, 2010, 2011, 2020])
    p.add_argument("--out", type=Path, default=Path("data/raw/worldpop"))
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    for y in a.years:
        download(y, a.out)


if __name__ == "__main__":
    main()
