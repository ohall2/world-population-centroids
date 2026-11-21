#!/usr/bin/env python3
"""Hämta GHSL-produkter (R2023A, WGS84 30 bågsekunder) direkt från JRC GHSL.

Jämförelsekällor för analys B (docs/ANALYSIS_PLAN.md §6). Produkterna laddas ner
i sitt eget rutnät (EPSG:4326, 30"), ingen pixel resamplas. Produkter:
  pop      GHS-POP, befolkning per cell
  built_s  GHS-BUILT-S, bebyggd yta (m²) per cell

För varje epok hämtas zip:en (återupptagbart), GeoTIFF:en packas upp, och en
JSON-fil skrivs bredvid med ursprung, rutnät, datatyp och SHA-256 (samma mönster
som tools/download_ee.py; pipelinen hittar filerna via JSON-filernas "year").
Finns .tif och .json redan hoppas epoken över. Zip-filen tas bort efter
uppackning (dess SHA-256 sparas i JSON-filen); behåll den med --keep-zip.

Källa och licens: JRC GHSL, release R2023A, licens CC BY 4.0 (se docs/DATA.md).

Exempel
  python tools/download_ghsl.py --product pop --epochs 2000 2010 2020
  python tools/download_ghsl.py --product built_s --epochs 2000 2010 2020
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

import rasterio

BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL"
PRODUCTS = {
    "pop": {"code": "GHS_POP", "out": "ghspop",
            "desc": "GHS-POP R2023A, WGS84 30 arcsec, residential population count per cell"},
    "built_s": {"code": "GHS_BUILT_S", "out": "ghsbuilts",
                "desc": "GHS-BUILT-S R2023A, WGS84 30 arcsec, built-up surface (m2) per cell"},
}


def stem(code: str, epoch: int) -> str:
    return f"{code}_E{epoch}_GLOBE_R2023A_4326_30ss_V1_0"


def url(code: str, epoch: int) -> str:
    return f"{BASE}/{code}_GLOBE_R2023A/{code}_E{epoch}_GLOBE_R2023A_4326_30ss/V1-0/{stem(code, epoch)}.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(product: str, epoch: int, out: Path, keep_zip: bool) -> None:
    code, desc = PRODUCTS[product]["code"], PRODUCTS[product]["desc"]
    tif = out / f"{stem(code, epoch)}.tif"
    js = out / f"{stem(code, epoch)}.json"
    if tif.exists() and js.exists():
        print(f"{product} {epoch}: finns redan, hoppar över")
        return
    zp = out / f"{stem(code, epoch)}.zip"
    print(f"{product} {epoch}: hämtar {url(code, epoch)}")
    subprocess.run(["curl", "-fSL", "-C", "-", "-o", str(zp), url(code, epoch)], check=True)
    sha_zip = sha256(zp)
    with zipfile.ZipFile(zp) as z:
        name = next(n for n in z.namelist() if n.endswith(".tif"))
        with z.open(name) as src, tif.open("wb") as dst:
            while chunk := src.read(1 << 20):
                dst.write(chunk)
    with rasterio.open(tif) as r:
        meta = {"source": "JRC GHSL (direct, jeodpp.jrc.ec.europa.eu)",
                "product": desc, "version": "R2023A V1-0", "license": "CC BY 4.0",
                "url": url(code, epoch), "epoch": epoch, "year": epoch,
                "crs": str(r.crs), "transform": list(r.transform)[:6],
                "width": r.width, "height": r.height, "dtype": r.dtypes[0],
                "nodata": r.nodata, "sha256_zip": sha_zip, "sha256": sha256(tif)}
    js.write_text(json.dumps(meta, indent=2))
    if not keep_zip:
        zp.unlink()
    print(f"{product} {epoch}: klar ({tif.name}, {tif.stat().st_size/1e9:.2f} GB)")


def main(argv=None):
    p = argparse.ArgumentParser(description="Hämta GHSL-produkter R2023A 30\" (WGS84).")
    p.add_argument("--product", choices=sorted(PRODUCTS), default="pop")
    p.add_argument("--epochs", type=int, nargs="+", default=[2000, 2010, 2020])
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--keep-zip", action="store_true")
    a = p.parse_args(argv)
    out = a.out or Path("data/raw") / PRODUCTS[a.product]["out"]
    out.mkdir(parents=True, exist_ok=True)
    for e in a.epochs:
        download(a.product, e, out, a.keep_zip)


if __name__ == "__main__":
    main()
