#!/usr/bin/env python3
"""Hämta GHS-POP R2023A, WGS84 30 bågsekunder, direkt från JRC GHSL.

Jämförelsekälla för analys B (docs/ANALYSIS_PLAN.md §6). Produkten laddas ner i
sitt eget rutnät (EPSG:4326, 30"), ingen pixel resamplas. För varje epok:
  1. zip hämtas (återupptagbart) från jeodpp.jrc.ec.europa.eu,
  2. GeoTIFF:en packas upp,
  3. en JSON-fil skrivs bredvid med ursprung, rutnät, datatyp och SHA-256
     (samma mönster som tools/download_ee.py; pipelinen hittar filerna via
     JSON-filernas "year").
Finns .tif och .json redan hoppas epoken över. Zip-filen tas bort efter
uppackning (dess SHA-256 sparas i JSON-filen); behåll den med --keep-zip.

Källa och licens: JRC GHSL, release R2023A, licens CC BY 4.0 (se docs/DATA.md).

Exempel
  python tools/download_ghspop.py --epochs 2000 2010 2020
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

import rasterio

BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A"
LICENSE = "CC BY 4.0"
PRODUCT = "GHS-POP R2023A, WGS84 30 arcsec, residential population count per cell"


def stem(epoch: int) -> str:
    return f"GHS_POP_E{epoch}_GLOBE_R2023A_4326_30ss_V1_0"


def url(epoch: int) -> str:
    return f"{BASE}/GHS_POP_E{epoch}_GLOBE_R2023A_4326_30ss/V1-0/{stem(epoch)}.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(epoch: int, out: Path, keep_zip: bool) -> None:
    tif = out / f"{stem(epoch)}.tif"
    js = out / f"{stem(epoch)}.json"
    if tif.exists() and js.exists():
        print(f"{epoch}: finns redan, hoppar över")
        return
    zp = out / f"{stem(epoch)}.zip"
    print(f"{epoch}: hämtar {url(epoch)}")
    subprocess.run(["curl", "-fSL", "-C", "-", "-o", str(zp), url(epoch)], check=True)
    sha_zip = sha256(zp)
    with zipfile.ZipFile(zp) as z:
        name = next(n for n in z.namelist() if n.endswith(".tif"))
        with z.open(name) as src, tif.open("wb") as dst:
            while chunk := src.read(1 << 20):
                dst.write(chunk)
    with rasterio.open(tif) as r:
        meta = {"source": "JRC GHSL (direct, jeodpp.jrc.ec.europa.eu)",
                "product": PRODUCT, "version": "R2023A V1-0", "license": LICENSE,
                "url": url(epoch), "epoch": epoch, "year": epoch,
                "crs": str(r.crs), "transform": list(r.transform)[:6],
                "width": r.width, "height": r.height, "dtype": r.dtypes[0],
                "nodata": r.nodata, "sha256_zip": sha_zip, "sha256": sha256(tif)}
    js.write_text(json.dumps(meta, indent=2))
    if not keep_zip:
        zp.unlink()
    print(f"{epoch}: klar ({tif.name}, {tif.stat().st_size/1e9:.2f} GB)")


def main(argv=None):
    p = argparse.ArgumentParser(description="Hämta GHS-POP R2023A 30\" (WGS84).")
    p.add_argument("--epochs", type=int, nargs="+", default=[2000, 2010, 2020])
    p.add_argument("--out", type=Path, default=Path("data/raw/ghspop"))
    p.add_argument("--keep-zip", action="store_true")
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    for e in a.epochs:
        download(e, a.out, a.keep_zip)


if __name__ == "__main__":
    main()
