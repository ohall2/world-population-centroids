#!/usr/bin/env python3
"""
World Population Centroids v2 – uppskatta subpixel-förskjutningar för alla
DMSP-kompositer mot VIIRS.

Metoden beskrivs i wpc/coregister.py. Referensen är den harmoniserade
VIIRS-bilden 2013 (data/derived/viirs_dmsp/VDMSP_2013.tif); dess läge
kommer enbart från VIIRS, eftersom harmoniseringen bara ändrar ljusvärden.

Utdata
  config/dmsp_offsets.csv          en rad per komposit: dx, dy (pixlar),
                                   spridning (viktad MAD), antal rutor,
                                   median av toppkorrelationen, samt dx/dy
                                   för Amerika (väster om 30° V) och övriga
                                   världen var för sig.
  data/derived/coreg_tiles.csv     förskjutning per ruta och komposit, för
                                   kontroll av regionala skillnader.

Exempel
  python tools/coregister_dmsp.py
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc import coregister as C  # noqa: E402
from wpc.lights import DMSP_HEIGHT, DMSP_WIDTH  # noqa: E402
from wpc.progress import Progress, fmt_seconds, log  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description="Subpixel-förskjutningar för DMSP-kompositer.")
    p.add_argument("--reference", type=Path, default=Path("data/derived/viirs_dmsp/VDMSP_2013.tif"))
    p.add_argument("--dmsp-dir", type=Path, default=Path("data/raw/dmsp"))
    p.add_argument("--out", type=Path, default=Path("config/dmsp_offsets.csv"))
    p.add_argument("--tiles-out", type=Path, default=Path("data/derived/coreg_tiles.csv"))
    args = p.parse_args(argv)

    files = sorted(args.dmsp_dir.glob("DMSP_F*_stable_lights.tif"))
    if not files:
        sys.exit(f"Inga DMSP-filer i {args.dmsp_dir}")
    t0 = time.time()
    tiles = [(r, c) for r in range(0, DMSP_HEIGHT - C.TILE + 1, C.TILE)
             for c in range(0, DMSP_WIDTH - C.TILE + 1, C.TILE)]
    use = []
    with rasterio.open(args.reference) as ref, Progress("Referensrutor", total=len(tiles)) as prog:
        for r, c in tiles:
            a = ref.read(1, window=Window(c, r, C.TILE, C.TILE))
            n = int((a > 0).sum())
            if n >= C.MIN_LIT:
                use.append((r, c, a, n))
            prog.advance()
    log(f"{len(use)} av {len(tiles)} rutor har tillräckligt med ljus i referensen")

    rows, tile_rows = [], []
    for f in files:
        sat, year = re.match(r"DMSP_(F\d\d)(\d{4})_", f.name).groups()
        res = []
        with rasterio.open(f) as src, Progress(f"{sat} {year}", total=len(use)) as prog:
            for r, c, a, n in use:
                b = src.read(1, window=Window(c, r, C.TILE, C.TILE))
                nb = int((b > 0).sum())
                if nb >= C.MIN_LIT:
                    o = C.tile_offset(a, b)
                    if o is not None:
                        res.append((o[0], o[1], o[2], min(n, nb), r, c))
                prog.advance()
        R = np.array(res)
        w = R[:, 3]
        lon = -180 + (R[:, 5] + C.TILE / 2) / 120
        am = lon < -30
        dx, dy = C.weighted_median(R[:, 0], w), C.weighted_median(R[:, 1], w)
        row = {"satellite": sat, "year": int(year), "dx": round(dx, 4), "dy": round(dy, 4),
               "mad_dx": round(C.weighted_median(np.abs(R[:, 0] - dx), w), 4),
               "mad_dy": round(C.weighted_median(np.abs(R[:, 1] - dy), w), 4),
               "n_tiles": len(R), "peak_median": round(float(np.median(R[:, 2])), 4),
               "dx_americas": round(C.weighted_median(R[am, 0], w[am]), 4),
               "dy_americas": round(C.weighted_median(R[am, 1], w[am]), 4),
               "dx_rest": round(C.weighted_median(R[~am, 0], w[~am]), 4),
               "dy_rest": round(C.weighted_median(R[~am, 1], w[~am]), 4)}
        rows.append(row)
        log(f"{sat} {year}: dx {dx:+.3f} dy {dy:+.3f} px, {len(R)} rutor")
        tile_rows += [{"satellite": sat, "year": int(year), "row": int(r_), "col": int(c_),
                       "dx": dx_, "dy": dy_, "peak": pk, "weight": int(w_)}
                      for dx_, dy_, pk, w_, r_, c_ in res]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("# Subpixel-förskjutning per DMSP-komposit mot VIIRS 2013 (tools/coregister_dmsp.py).\n"
                 "# referens(x, y) ≈ komposit(x + dx, y + dy); x = kolumn (öster), y = rad (söder).\n"
                 f"# Referens: {args.reference.name}; ruta {C.TILE} px, sökområde ±{C.MAX_SHIFT} px.\n")
        pd.DataFrame(rows).to_csv(fh, index=False)
    args.tiles_out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(tile_rows).to_csv(args.tiles_out, index=False, float_format="%.4f")
    log(f"Klart på {fmt_seconds(time.time() - t0)}. {args.out}")


if __name__ == "__main__":
    main()
