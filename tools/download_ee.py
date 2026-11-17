#!/usr/bin/env python3
"""
World Population Centroids v2 – hämta nattljusdata från Google Earth Engine.

Produkter (--product):
  dmsp   NOAA/DMSP-OLS/NIGHTTIME_LIGHTS, bandet stable_lights (1992–2013)
  viirs  NOAA/VIIRS/DNB/ANNUAL_V22, banden average_masked och cf_cvg

CCNL hämtas inte härifrån: BNU/FGS/CCNL/v1 i Earth Engine är omsamplad till
0,008983° (1 km vid ekvatorn). Originalet i 30" hämtas med
tools/download_zenodo.py --product ccnl.

Skriptet laddar ner banden direkt till disk, utan Google Drive och utan
kodeditorn. Varje bild och band blir en GeoTIFF i källans eget rutnät
(samma CRS, samma affina transform, samma antal rader och kolumner), så att
inga pixlar resamplas.

Hur det fungerar
  1. Bilderna för de valda åren listas ur katalogen.
  2. För varje bild och band läses källans rutnät och datatyp.
  3. Rutnätet delas i rutor (ca 16 MB per ruta) som hämtas med
     ee.data.computePixels i formatet NUMPY_NDARRAY, flera samtidigt.
     Varje anrop har en tidsgräns på 2 minuter och görs om vid fel.
  4. Rutorna skrivs på rätt plats i en GeoTIFF (DEFLATE-komprimerad).
  5. Bredvid varje GeoTIFF skrivs en JSON-fil med ursprung, rutnät och
     SHA-256-checksumma. Finns båda filerna redan hoppas filen över, så ett
     avbrutet körning kan startas om med samma kommando.

Maskade pixlar sätts till 0 (unmask(0)), vilket motsvarar EOG:s
originalfiler där pixlar utan ljus har värdet 0.

Under nedladdningen visas en rad med snurra, antal rutor, tid och
beräknad återstående tid.

Exempel
  # Visa vad som finns (år, band, rutnät, datatyp) utan att hämta:
  python tools/download_ee.py --product viirs --project mitt-projekt --dry-run

  # Hämta alla år av VIIRS till data/raw/viirs:
  python tools/download_ee.py --product viirs --project mitt-projekt
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

import numpy as np
import rasterio
from rasterio.transform import Affine
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.progress import Progress, fmt_seconds, log  # noqa: E402

# Produkter som kan hämtas. "prefix" ger filnamnet <prefix>_<id>_<band>.tif;
# för DMSP behålls namnen från första nedladdningen.
PRODUCTS = {
    "dmsp": {"collection": "NOAA/DMSP-OLS/NIGHTTIME_LIGHTS",
             "bands": ["stable_lights"], "prefix": "DMSP"},
    "viirs": {"collection": "NOAA/VIIRS/DNB/ANNUAL_V22",
              "bands": ["average_masked", "cf_cvg"], "prefix": "VIIRS"},
    # Publicerade harmoniserade serier, för jämförelse (docs/HARMONIZATION.md).
    # bands None = alla band i bilden (namnen läses från katalogen).
    "li2020-dmsp": {"collection": "projects/sat-io/open-datasets/Harmonized_NTL/dmsp",
                    "bands": None, "prefix": "LI2020"},
    "li2020-viirs": {"collection": "projects/sat-io/open-datasets/Harmonized_NTL/viirs",
                     "bands": None, "prefix": "LI2020"},
    "svnl": {"collection": "projects/sat-io/open-datasets/srunet-npp-viirs-ntl",
             "bands": None, "prefix": "SVNL"},
}
HIGHVOLUME_URL = "https://earthengine-highvolume.googleapis.com"
REQUEST_TIMEOUT_MS = 120_000     # ett anrop som inte svarat på 2 min avbryts
TILE_BYTES = 16 * 2**20          # rutstorlek styrs av ca 16 MB per ruta

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
# Datatyper
# --------------------------------------------------------------------------

_CASTERS = {"uint8": "toUint8", "uint16": "toUint16", "int16": "toInt16",
            "int32": "toInt32", "float32": "toFloat"}


def numpy_dtype(data_type: dict) -> str:
    """Minsta numpy-typ som rymmer bandets värden enligt Earth Engine."""
    if data_type.get("precision") != "int":
        return "float32"
    lo, hi = data_type.get("min", -2**31), data_type.get("max", 2**31 - 1)
    if lo >= 0 and hi <= 255:
        return "uint8"
    if lo >= 0 and hi <= 65535:
        return "uint16"
    if lo >= -32768 and hi <= 32767:
        return "int16"
    if lo >= -2**31 and hi <= 2**31 - 1:
        return "int32"
    raise ValueError(f"Heltalsintervall {lo}..{hi} ryms inte i int32")


def tile_size(dtype: str) -> int:
    """Rutsida i pixlar, multipel av 256, så att en ruta blir ca TILE_BYTES."""
    side = int(math.sqrt(TILE_BYTES / np.dtype(dtype).itemsize))
    return max(256, side // 256 * 256)


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
        "compress": "deflate", "predictor": 3 if dtype == "float32" else 2,
        "BIGTIFF": "IF_SAFER",
    }
    tmp = out_tif.with_suffix(out_tif.suffix + ".part")
    with Progress(out_tif.name, total=len(tiles)) as prog, \
            rasterio.open(tmp, "w", **profile) as dst, \
            ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch, *t): t for t in tiles}
        for fut in as_completed(futures):
            col0, row0, w, h = futures[fut]
            arr = fut.result()
            if arr.shape != (h, w):
                raise RuntimeError(
                    f"Ruta col0={col0} row0={row0}: fick {arr.shape}, väntade {(h, w)}")
            dst.write(arr.astype(dtype, copy=False), 1,
                      window=Window(col0, row0, w, h))
            prog.advance()
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
            "projekt-ID:t är rätt och registrerat för Earth Engine.")
    # Utan tidsgräns kan ett anrop hänga i det oändliga när uppkopplingen
    # tappas; då avbryts det i stället och görs om.
    if hasattr(ee.data, "setDeadline"):
        ee.data.setDeadline(REQUEST_TIMEOUT_MS)


def list_images(collection: str, years: list[int] | None) -> list[dict]:
    """Lista bilderna i katalogen, eventuellt filtrerade på år."""
    import ee
    coll = ee.ImageCollection(collection)
    ids = coll.aggregate_array("system:index").getInfo()
    starts = coll.aggregate_array("system:time_start").getInfo()
    images = []
    for img_id, t in zip(ids, starts):
        year = dt.datetime.fromtimestamp(t / 1000, tz=dt.timezone.utc).year
        if years is None or year in years:
            images.append({"id": img_id, "year": year})
    images.sort(key=lambda d: (d["year"], d["id"]))
    return images


def band_info(asset_id: str, band: str) -> dict:
    """Rutnät och datatyp för ett band, ur katalogen."""
    import ee
    info = ee.Image(asset_id).select(band).getInfo()["bands"][0]
    if "dimensions" not in info:
        raise RuntimeError(f"{asset_id}/{band}: saknar fasta dimensioner i katalogen.")
    return {"crs": info["crs"],
            "transform": [float(v) for v in info["crs_transform"]],
            "width": int(info["dimensions"][0]),
            "height": int(info["dimensions"][1]),
            "data_type": info.get("data_type", {}),
            "dtype": numpy_dtype(info.get("data_type", {}))}


def fetch_tile_factory(asset_id: str, band: str, crs: str,
                       transform: list[float], dtype: str,
                       retries: int = 6) -> Callable[[int, int, int, int], np.ndarray]:
    import ee
    image = getattr(ee.Image(asset_id).select(band).unmask(0), _CASTERS[dtype])()

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
                log(f"  ruta col0={col0} row0={row0}: {str(exc)[:120]} "
                    f"– nytt försök om {wait} s")
                time.sleep(wait)
        raise AssertionError("unreachable")

    return fetch


def out_name(prefix: str, img_id: str, band: str) -> str:
    return f"{prefix}_{img_id}_{band}.tif"


def download_band(img: dict, band: str, collection: str, prefix: str,
                  out_dir: Path, tile: int | None, workers: int,
                  project: str) -> None:
    import ee
    asset_id = f"{collection}/{img['id']}"
    out_tif = out_dir / out_name(prefix, img["id"], band)
    out_json = out_tif.with_suffix(".json")
    if out_tif.exists() and out_json.exists():
        log(f"{out_tif.name} finns redan – hoppar över")
        return

    g = band_info(asset_id, band)
    tile = tile or tile_size(g["dtype"])
    log(f"{out_tif.name}: {g['width']} x {g['height']} pixlar, {g['dtype']}, "
        f"rutor {tile} px")
    fetch = fetch_tile_factory(asset_id, band, g["crs"], g["transform"], g["dtype"])
    t0 = time.time()
    n_tiles = assemble(out_tif, g["crs"], g["transform"], g["width"], g["height"],
                       g["dtype"], tile, workers, fetch)

    meta = {
        "source": "Google Earth Engine",
        "collection": collection,
        "asset_id": asset_id,
        "band": band,
        "year": img["year"],
        "crs": g["crs"],
        "crs_transform": g["transform"],
        "width": g["width"],
        "height": g["height"],
        "ee_data_type": g["data_type"],
        "dtype": g["dtype"],
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


def dry_run(images: list[dict], collection: str, bands: list[str]) -> None:
    """Skriv ut år, band, rutnät och datatyp, och uppskattad storlek."""
    seen = set()
    for img in images:
        for band in bands:
            g = band_info(f"{collection}/{img['id']}", band)
            key = (band, g["width"], g["height"], tuple(g["transform"]), g["dtype"])
            mark = "" if key in seen else "   <- nytt rutnät/typ"
            seen.add(key)
            gb = g["width"] * g["height"] * np.dtype(g["dtype"]).itemsize / 1e9
            print(f"  {img['id']:>12} {img['year']}  {band:<15} "
                  f"{g['width']} x {g['height']}  {g['dtype']:<7} "
                  f"({g['data_type'].get('precision')}, "
                  f"{g['data_type'].get('min')}..{g['data_type'].get('max')})  "
                  f"okomprimerat {gb:.1f} GB  transform {g['transform']}{mark}")


# --------------------------------------------------------------------------
# Kommandorad
# --------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Hämta nattljusdata från Earth Engine till disk.")
    p.add_argument("--product", choices=sorted(PRODUCTS), default="dmsp")
    p.add_argument("--project", required=True,
                   help="Ditt Google Cloud-projekt-ID för Earth Engine.")
    p.add_argument("--out", default=None,
                   help="Mapp för GeoTIFF-filerna (standard data/raw/<produkt>; li2020-* -> data/raw/li2020).")
    p.add_argument("--years", nargs="+", default=["all"],
                   help="År att hämta, t.ex. 1992 1993 eller 1992-2000, eller "
                        "'all' (standard).")
    p.add_argument("--bands", nargs="+", default=None,
                   help="Band att hämta (standard: produktens band).")
    p.add_argument("--collection", default=None,
                   help="Annan katalogsamling än produktens standard.")
    p.add_argument("--tile", type=int, default=None,
                   help="Rutstorlek i pixlar (standard: ca 16 MB per ruta).")
    p.add_argument("--workers", type=int, default=8,
                   help="Antal rutor som hämtas samtidigt (standard 8).")
    p.add_argument("--no-highvolume", action="store_true",
                   help="Använd standardslutpunkten i stället för high-volume.")
    p.add_argument("--dry-run", action="store_true",
                   help="Visa vad som finns, utan att hämta något.")
    return p.parse_args(argv)


def parse_years(tokens: list[str]) -> list[int] | None:
    if tokens == ["all"]:
        return None
    years: list[int] = []
    for t in tokens:
        if "-" in t:
            a, b = (int(x) for x in t.split("-"))
            years.extend(range(a, b + 1))
        else:
            years.append(int(t))
    return sorted(set(years))


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    prod = PRODUCTS[args.product]
    collection = args.collection or prod["collection"]
    bands = args.bands or prod["bands"]
    years = parse_years(args.years)

    init_ee(args.project, highvolume=not args.no_highvolume)
    images = list_images(collection, years)
    if not images:
        sys.exit(f"Inga bilder i {collection} för de valda åren.")
    if bands is None:
        import ee
        bands = ee.Image(f"{collection}/{images[0]['id']}").bandNames().getInfo()
    log(f"{collection}: {len(images)} bilder, år "
        f"{images[0]['year']}–{images[-1]['year']}, band {', '.join(bands)}")
    if args.dry_run:
        dry_run(images, collection, bands)
        return

    out_dir = Path(args.out or f"data/raw/{args.product.split('-')[0]}")
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    for img in images:
        for band in bands:
            download_band(img, band, collection, prod["prefix"], out_dir,
                          args.tile, args.workers, args.project)
    log(f"Klart på {fmt_seconds(time.time() - t0)}. Filerna ligger i {out_dir.resolve()}")


if __name__ == "__main__":
    main()
