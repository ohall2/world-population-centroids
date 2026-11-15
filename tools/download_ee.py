#!/usr/bin/env python3
"""
World Population Centroids v2 – hämta DMSP-OLS v4 från Google Earth Engine.

Skriptet laddar ner bandet `stable_lights` ur Earth Engine-katalogen
NOAA/DMSP-OLS/NIGHTTIME_LIGHTS direkt till disk, utan Google Drive och utan
kodeditorn. Varje satellit-år blir en GeoTIFF i källans eget rutnät
(samma CRS, samma affina transform, samma antal rader och kolumner), så att
inga pixlar resamplas.

Hur det fungerar
  1. Kompositerna för de valda åren listas ur katalogen.
  2. För varje komposit läses källans rutnät (crs, crs_transform, dimensions).
  3. Rutnätet delas i rutor (standard 4096 x 4096 pixlar) som hämtas med
     ee.data.computePixels i formatet NUMPY_NDARRAY, flera samtidigt.
  4. Rutorna skrivs på rätt plats i en GeoTIFF (DEFLATE-komprimerad).
  5. Bredvid varje GeoTIFF skrivs en JSON-fil med ursprung, rutnät och
     SHA-256-checksumma. Finns båda filerna redan hoppas kompositen över.

Maskade pixlar sätts till 0 (unmask(0)), vilket motsvarar EOG:s
originalfiler där pixlar utan stabilt ljus har värdet 0.

Exempel
  # De sju kompositerna som Fas 1 behöver (1992, 1994, 2004, 2010, 2013):
  python download_ee.py --project mitt-projekt-id --out data/raw/dmsp

  # Alla år 1992–2013:
  python download_ee.py --project mitt-projekt-id --out data/raw/dmsp --years all

  # Bara lista vad som skulle hämtas, utan att hämta:
  python download_ee.py --project mitt-projekt-id --dry-run
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

import numpy as np
import rasterio
from rasterio.transform import Affine
from rasterio.windows import Window

COLLECTION = "NOAA/DMSP-OLS/NIGHTTIME_LIGHTS"
DEFAULT_BAND = "stable_lights"
PHASE1_YEARS = [1992, 1994, 2004, 2010, 2013]
HIGHVOLUME_URL = "https://earthengine-highvolume.googleapis.com"

log = logging.getLogger("download_ee")


# --------------------------------------------------------------------------
# Rutnät och rutor
# --------------------------------------------------------------------------

def make_tiles(width: int, height: int, tile: int) -> list[tuple[int, int, int, int]]:
    """Dela ett rutnät i rutor. Returnerar (col0, row0, w, h) per ruta."""
    tiles = []
    for row0 in range(0, height, tile):
        for col0 in range(0, width, tile):
            tiles.append((col0, row0, min(tile, width - col0), min(tile, height - row0)))
    return tiles


def tile_grid(crs: str, transform: list[float], col0: int, row0: int,
              w: int, h: int) -> dict:
    """Rutnät för en ruta i källans pixelsystem, i computePixels-format.

    transform = [scaleX, shearX, translateX, shearY, scaleY, translateY].
    Rutans övre vänstra hörn flyttas exakt ett heltal pixlar, så rutorna
    ligger i samma rutnät som källan.
    """
    sx, shx, tx, shy, sy, ty = transform
    if shx != 0 or shy != 0:
        raise ValueError("Rutnät med skjuvning stöds inte.")
    return {
        "dimensions": {"width": w, "height": h},
        "affineTransform": {
            "scaleX": sx, "shearX": 0, "translateX": tx + col0 * sx,
            "shearY": 0, "scaleY": sy, "translateY": ty + row0 * sy,
        },
        "crsCode": crs,
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# Sammansättning (oberoende av Earth Engine, så att den kan testas lokalt)
# --------------------------------------------------------------------------

def assemble(out_tif: Path, crs: str, transform: list[float], width: int,
             height: int, dtype: str, tile: int, workers: int,
             fetch: Callable[[int, int, int, int], np.ndarray]) -> int:
    """Hämta alla rutor med `fetch` och skriv dem till en GeoTIFF.

    Skrivningen sker i huvudtråden; bara hämtningen är parallell.
    Returnerar antalet rutor.
    """
    tiles = make_tiles(width, height, tile)
    profile = {
        "driver": "GTiff", "dtype": dtype, "count": 1,
        "width": width, "height": height,
        "crs": crs, "transform": Affine(*transform),
        "tiled": True, "blockxsize": 512, "blockysize": 512,
        "compress": "deflate", "predictor": 2, "BIGTIFF": "IF_SAFER",
    }
    tmp = out_tif.with_suffix(out_tif.suffix + ".part")
    with rasterio.open(tmp, "w", **profile) as dst:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(fetch, *t): t for t in tiles}
            for i, fut in enumerate(as_completed(futures), start=1):
                col0, row0, w, h = futures[fut]
                arr = fut.result()
                if arr.shape != (h, w):
                    raise RuntimeError(
                        f"Ruta col0={col0} row0={row0}: fick {arr.shape}, väntade {(h, w)}")
                dst.write(arr.astype(dtype, copy=False), 1,
                          window=Window(col0, row0, w, h))
                if i % 10 == 0 or i == len(tiles):
                    log.info("  %d/%d rutor klara", i, len(tiles))
    tmp.replace(out_tif)
    return len(tiles)


# --------------------------------------------------------------------------
# Earth Engine
# --------------------------------------------------------------------------

def init_ee(project: str, highvolume: bool) -> None:
    import ee
    try:
        if highvolume:
            ee.Initialize(project=project, url=HIGHVOLUME_URL)
        else:
            ee.Initialize(project=project)
    except Exception as exc:  # noqa: BLE001
        sys.exit(
            "Kunde inte ansluta till Earth Engine.\n"
            f"Fel: {exc}\n\n"
            "Kontrollera att du har kört 'earthengine authenticate' och att\n"
            "projekt-ID:t är rätt och registrerat för Earth Engine."
        )


def list_images(years: list[int] | None) -> list[dict]:
    """Lista kompositerna i katalogen, eventuellt filtrerade på år."""
    import ee
    coll = ee.ImageCollection(COLLECTION)
    ids = coll.aggregate_array("system:index").getInfo()
    starts = coll.aggregate_array("system:time_start").getInfo()
    images = []
    for img_id, t in zip(ids, starts):
        year = dt.datetime.fromtimestamp(t / 1000, tz=dt.timezone.utc).year
        if years is None or year in years:
            images.append({"id": img_id, "year": year})
    images.sort(key=lambda d: (d["year"], d["id"]))
    return images


def fetch_tile_factory(asset_id: str, band: str, crs: str,
                       transform: list[float], ee_dtype: str,
                       retries: int = 6) -> Callable[[int, int, int, int], np.ndarray]:
    import ee
    caster = {"uint8": "toUint8", "uint16": "toUint16",
              "float32": "toFloat"}[ee_dtype]
    image = getattr(ee.Image(asset_id).select(band).unmask(0), caster)()

    def fetch(col0: int, row0: int, w: int, h: int) -> np.ndarray:
        request = {
            "expression": image,
            "fileFormat": "NUMPY_NDARRAY",
            "grid": tile_grid(crs, transform, col0, row0, w, h),
        }
        for attempt in range(1, retries + 1):
            try:
                arr = ee.data.computePixels(request)
                return np.asarray(arr[band])
            except Exception as exc:  # noqa: BLE001
                if attempt == retries:
                    raise RuntimeError(
                        f"Ruta col0={col0} row0={row0} misslyckades efter "
                        f"{retries} försök: {exc}") from exc
                wait = min(60, 2 ** attempt)
                log.warning("  ruta col0=%d row0=%d: %s – nytt försök om %d s",
                            col0, row0, exc, wait)
                time.sleep(wait)
        raise AssertionError("unreachable")

    return fetch


def download_image(img: dict, band: str, out_dir: Path, tile: int,
                   workers: int, project: str) -> None:
    import ee
    asset_id = f"{COLLECTION}/{img['id']}"
    out_tif = out_dir / f"DMSP_{img['id']}_{band}.tif"
    out_json = out_tif.with_suffix(".json")
    if out_tif.exists() and out_json.exists():
        log.info("%s finns redan – hoppar över", out_tif.name)
        return

    info = ee.Image(asset_id).select(band).getInfo()["bands"][0]
    if "dimensions" not in info:
        raise RuntimeError(f"{asset_id}: saknar fasta dimensioner i katalogen.")
    crs = info["crs"]
    transform = [float(v) for v in info["crs_transform"]]
    width, height = (int(v) for v in info["dimensions"])
    dtype = "uint8"  # stable_lights är DN 0–63

    log.info("%s: %d x %d pixlar, %s, transform %s",
             img["id"], width, height, crs, transform)
    fetch = fetch_tile_factory(asset_id, band, crs, transform, dtype)
    t0 = time.time()
    n_tiles = assemble(out_tif, crs, transform, width, height, dtype,
                       tile, workers, fetch)

    meta = {
        "source": "Google Earth Engine",
        "asset_id": asset_id,
        "band": band,
        "year": img["year"],
        "crs": crs,
        "crs_transform": transform,
        "width": width,
        "height": height,
        "dtype": dtype,
        "masked_pixels_set_to": 0,
        "tile_size": tile,
        "n_tiles": n_tiles,
        "ee_project": project,
        "earthengine_api_version": ee.__version__,
        "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "seconds": round(time.time() - t0, 1),
        "sha256": sha256(out_tif),
    }
    out_json.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    log.info("%s klar (%.0f s)", out_tif.name, meta["seconds"])


# --------------------------------------------------------------------------
# Kommandorad
# --------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Hämta DMSP-OLS v4 stable_lights från Earth Engine till disk.")
    p.add_argument("--project", required=True,
                   help="Ditt Google Cloud-projekt-ID för Earth Engine.")
    p.add_argument("--out", default="data/raw/dmsp",
                   help="Mapp för GeoTIFF-filerna (skapas vid behov).")
    p.add_argument("--years", nargs="+", default=[str(y) for y in PHASE1_YEARS],
                   help="År att hämta, eller 'all'. Standard: Fas 1-åren "
                        "1992 1994 2004 2010 2013.")
    p.add_argument("--band", default=DEFAULT_BAND)
    p.add_argument("--tile", type=int, default=4096,
                   help="Rutstorlek i pixlar (standard 4096).")
    p.add_argument("--workers", type=int, default=8,
                   help="Antal rutor som hämtas samtidigt (standard 8).")
    p.add_argument("--no-highvolume", action="store_true",
                   help="Använd standardslutpunkten i stället för high-volume.")
    p.add_argument("--dry-run", action="store_true",
                   help="Lista kompositerna utan att hämta något.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    years = None if args.years == ["all"] else [int(y) for y in args.years]

    init_ee(args.project, highvolume=not args.no_highvolume)
    images = list_images(years)
    if not images:
        sys.exit("Inga kompositer hittades för de valda åren.")
    log.info("%d kompositer: %s", len(images),
             ", ".join(f"{d['id']} ({d['year']})" for d in images))
    if args.dry_run:
        return

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    for img in images:
        download_image(img, args.band, out_dir, args.tile, args.workers,
                       args.project)
    log.info("Klart. Filerna ligger i %s", out_dir.resolve())


if __name__ == "__main__":
    main()
