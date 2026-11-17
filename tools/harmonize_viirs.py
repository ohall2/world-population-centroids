#!/usr/bin/env python3
"""
World Population Centroids v2 – harmonisera VIIRS till DMSP-liknande värden.

Metoden beskrivs i wpc/harmonize.py och docs/METHODS.md. I korthet: VIIRS
aggregeras till DMSP:s 30"-rutnät, smetas ut med en gaussisk kärna och förs
över till DMSP-IRQR-skalan med en monoton funktion som anpassas på
överlappsåret 2013.

Steg
  1. Anpassning (2013). För varje kandidat-sigma samlas klasstatistik över
     hela jorden. Utsmetningen väljs med tvåfaldig rumslig korsvalidering:
     jorden delas i 10°-band i longitud, varannan i varje vikning; funktionen
     anpassas på den ena halvan och prövas på den andra. Vald sigma ger
     minst kvadratfel på testhalvorna. Funktionen anpassas sedan på alla
     pixlar med vald sigma.
  2. Tillämpning. Varje VIIRS-år (också 2013, för kontroll) omvandlas och
     skrivs som GeoTIFF i DMSP:s rutnät till data/derived/viirs_dmsp/
     VDMSP_<år>.tif, med en JSON-fil bredvid.

Fackelmasken (standard data/processed/flares_viirs_r2km.gpkg, som också
innehåller 2009-polygonerna) sätter VIIRS-pixlar till 0 före aggregeringen,
och DMSP-pixlar inom masken tas inte med i anpassningen.

Exempel
  python tools/harmonize_viirs.py              # anpassa och tillämpa alla år
  python tools/harmonize_viirs.py --fit-only
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc import harmonize as H  # noqa: E402
from wpc.lights import (DMSP_HEIGHT, DMSP_TRANSFORM, DMSP_WIDTH, DmspIrqr,  # noqa: E402
                        Window2D, annual_files, mask_flares, read_offsets, snap_transform)
from wpc.progress import Progress, fmt_seconds, log  # noqa: E402

SIGMAS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
BLOCK = 256                       # DMSP-rader per block
BAND_DEG = 10                     # bredd på longitudbanden i korsvalideringen


class ViirsReader:
    """Läser VIIRS-rader 2a-1 ... 2b+1 (fyller utanför bilden med 0) och
    maskar facklor i VIIRS-upplösning."""

    def __init__(self, path: Path, flares, index):
        self.src = rasterio.open(path)
        t = snap_transform(self.src.transform)
        if (self.src.width, self.src.height) != (2 * DMSP_WIDTH - 1, 2 * DMSP_HEIGHT - 1):
            raise ValueError(f"{path.name}: oväntad storlek {self.src.width} x {self.src.height}")
        if abs(t.c - (DMSP_TRANSFORM.c + 1 / 480)) > 1e-9 or abs(t.f - (DMSP_TRANSFORM.f - 1 / 480)) > 1e-9:
            raise ValueError(f"{path.name}: rutnätet ligger inte som väntat i förhållande till DMSP")
        self.t, self.flares, self.index = t, flares, index

    def rows(self, a: int, b: int) -> np.ndarray:
        r0, r1 = 2 * a - 1, 2 * b + 2            # [r0, r1)
        out = np.zeros((r1 - r0, self.src.width), "float32")
        lo, hi = max(r0, 0), min(r1, self.src.height)
        if hi > lo:
            v = self.src.read(1, window=Window(0, lo, self.src.width, hi - lo)).astype("float32")
            v[~np.isfinite(v) | (v < 0)] = 0.0
            w = Window2D(data=v, transform=self.t @ rasterio.Affine.translation(0, lo))
            mask_flares(w, self.flares, self.index)
            out[lo - r0:hi - r0] = w.data
        return out

    def close(self):
        self.src.close()


ROWS = (0, DMSP_HEIGHT)            # kan begränsas med --rows (för test)


def blocks(halo: int):
    for j0 in range(ROWS[0], ROWS[1], BLOCK):
        j1 = min(j0 + BLOCK, ROWS[1])
        yield j0, j1, j0 - halo, j1 - 1 + halo


def transformed_block(reader: ViirsReader, j0, j1, a, b, sigmas):
    """Aggregerad och utsmetad VIIRS för DMSP-raderna j0..j1-1, per sigma."""
    agg = H.aggregate_to_dmsp(reader.rows(a, b))
    return {s: H.blur(agg, s)[j0 - a:j1 - a] for s in sigmas}


def band_fold(width: int) -> np.ndarray:
    lon = DMSP_TRANSFORM.c + (np.arange(width) + 0.5) * DMSP_TRANSFORM.a
    return (np.floor((lon + 180) / BAND_DEG).astype(int) % 2).astype(bool)


def fit(viirs_2013: Path, dmsp: DmspIrqr, flares, index) -> dict:
    halo = math.ceil(3 * max(SIGMAS))
    stats = {s: (H.BinStats(), H.BinStats()) for s in SIGMAS}
    fold = band_fold(DMSP_WIDTH)
    reader = ViirsReader(viirs_2013, flares, index)
    blist = list(blocks(halo))
    with Progress("Anpassning 2013 (block)", total=len(blist)) as prog:
        for j0, j1, a, b in blist:
            y = dmsp.read_window(Window(0, j0, DMSP_WIDTH, j1 - j0))   # maskad
            keep = Window2D(data=np.ones_like(y.data), transform=y.transform)
            mask_flares(keep, flares, index)
            keep = keep.data > 0
            xs = transformed_block(reader, j0, j1, a, b, SIGMAS)
            for s, x in xs.items():
                idx = H.bin_index(x)
                for k, f in enumerate((fold, ~fold)):
                    sel = keep & f[None, :]
                    stats[s][k].add(idx[sel], y.data[sel])
            prog.advance()
    reader.close()

    cv = {}
    for s, (A, B) in stats.items():
        e = H.sse(B, H.fit_table(A)["value"]) + H.sse(A, H.fit_table(B)["value"])
        allst = A.merged(B)
        sst = float(np.sum(allst.ss) - np.sum(allst.s) ** 2 / np.sum(allst.n))
        cv[s] = {"sse": e, "r2": 1 - e / sst}
        log(f"  sigma {s:.1f} DMSP-pixlar: korsvaliderat R² = {cv[s]['r2']:.4f}")
    best = min(cv, key=lambda s: cv[s]["sse"])
    allst = stats[best][0].merged(stats[best][1])
    table = H.fit_table(allst)
    log(f"Vald sigma: {best} DMSP-pixlar (ca {best * 0.93:.1f} km vid ekvatorn)")
    return {"sigma": best, "cv": {str(k): v for k, v in cv.items()},
            "bin": {"log_min": H.LOG_MIN, "log_max": H.LOG_MAX, "log_step": H.LOG_STEP},
            "table": {k: v.tolist() for k, v in table.items()},
            "stats": allst.to_dict()}


def apply_year(year: int, path: Path, model: dict, out_dir: Path, flares, index) -> None:
    out = out_dir / f"VDMSP_{year}.tif"
    js = out.with_suffix(".json")
    if out.exists() and js.exists():
        log(f"{out.name} finns redan – hoppar över")
        return
    sigma, value = model["sigma"], np.asarray(model["table"]["value"])
    halo = math.ceil(3 * sigma)
    reader = ViirsReader(path, flares, index)
    prof = {"driver": "GTiff", "dtype": "float32", "count": 1, "width": DMSP_WIDTH,
            "height": DMSP_HEIGHT, "crs": "EPSG:4326", "transform": DMSP_TRANSFORM,
            "tiled": True, "blockxsize": 512, "blockysize": 512, "compress": "deflate",
            "predictor": 3, "BIGTIFF": "IF_SAFER"}
    tmp = out.with_name(out.name + ".part")
    t0 = time.time()
    blist = list(blocks(halo))
    with Progress(f"{out.name} (block)", total=len(blist)) as prog, rasterio.open(tmp, "w", **prof) as dst:
        for j0, j1, a, b in blist:
            x = transformed_block(reader, j0, j1, a, b, [sigma])[sigma]
            dst.write(H.apply_table(x, value), 1, window=Window(0, j0, DMSP_WIDTH, j1 - j0))
            prog.advance()
    reader.close()
    tmp.replace(out)
    js.write_text(json.dumps({"year": year, "source_file": path.name, "sigma": sigma,
                              "model": "data/derived/viirs_dmsp/model.json",
                              "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                              "seconds": round(time.time() - t0, 1)}, indent=2) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description="Harmonisera VIIRS till DMSP-liknande värden.")
    p.add_argument("--viirs-dir", type=Path, default=Path("data/raw/viirs"))
    p.add_argument("--dmsp-dir", type=Path, default=Path("data/raw/dmsp"))
    p.add_argument("--calibration", type=Path, default=Path("config/calibration_v1.csv"))
    p.add_argument("--flares", type=Path, default=Path("data/processed/flares_viirs_r2km.gpkg"))
    p.add_argument("--out-dir", type=Path, default=Path("data/derived/viirs_dmsp"))
    p.add_argument("--dmsp-offsets", type=Path, default=None,
                   help="Samregistrera DMSP 2013 med förskjutningarna i filen (config/dmsp_offsets.csv) "
                        "före anpassningen.")
    p.add_argument("--fit-only", action="store_true")
    p.add_argument("--refit", action="store_true", help="Anpassa om även om model.json finns.")
    p.add_argument("--rows", type=int, nargs=2, default=None, help=argparse.SUPPRESS)
    p.add_argument("--years", type=int, nargs="+", default=None, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    global ROWS
    if args.rows:
        ROWS = tuple(args.rows)

    files = annual_files("viirs", args.viirs_dir)
    if 2013 not in files:
        sys.exit(f"VIIRS 2013 saknas i {args.viirs_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    flares = gpd.read_file(args.flares, engine="pyogrio")
    index = flares.sindex
    t0 = time.time()
    model_path = args.out_dir / "model.json"
    if model_path.exists() and not args.refit:
        model = json.loads(model_path.read_text())
        log(f"Använder befintlig modell {model_path} (sigma {model['sigma']})")
    else:
        offsets = read_offsets(args.dmsp_offsets) if args.dmsp_offsets else None
        dmsp = DmspIrqr(2013, args.dmsp_dir, args.calibration, args.flares, offsets=offsets)
        model = fit(files[2013], dmsp, flares, index)
        model.update({"overlap_year": 2013, "dmsp": dmsp.describe(), "viirs_2013": files[2013].name,
                      "flares": str(args.flares),
                      "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
        dmsp.close()
        model_path.write_text(json.dumps(model, indent=1, ensure_ascii=False) + "\n")
    if not args.fit_only:
        for year in sorted(args.years or files):
            apply_year(year, files[year], model, args.out_dir, flares, index)
    log(f"Klart på {fmt_seconds(time.time() - t0)}.")


if __name__ == "__main__":
    main()
