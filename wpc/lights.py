"""Ljuskällor för centroidberäkningen.

En ljuskälla levererar, för ett givet geografiskt område, ett fönster ur ett
globalt rutnät med ljusvärden (float32) och fönstrets affina transform.
Allt annat i pipelinen är oberoende av vilken källa som används.

Implementerat:
  DmspIrqr  DMSP-OLS v4 Stable Lights, kalibrerade med IRQR-koefficienterna i
            Hall m.fl. (2019, tabell 2), medelvärde per år, gasfacklor maskade.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio.features import rasterize
from rasterio.windows import Window

# DMSP-OLS v4: 43 201 x 16 801 pixlar à 30 bågsekunder. Pixelmittpunkterna
# ligger på hela 30"-steg från -180 till 180 och från 75 till -65, så hörnet
# ligger en halv pixel utanför. Earth Engine anger hörnet till (-180, 75) och
# skalan avrundad till 0.0083333333; vi använder EOG:s rutnät.
DMSP_WIDTH, DMSP_HEIGHT = 43201, 16801
DMSP_TRANSFORM = Affine(1 / 120, 0.0, -180.0 - 1 / 240, 0.0, -1 / 120, 75.0 + 1 / 240)


@dataclass
class Window2D:
    """Ett utsnitt ur ett globalt rutnät."""
    data: np.ndarray          # float32, rader x kolumner
    transform: Affine         # fönstrets transform

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        h, w = self.data.shape
        xmin, ymax = self.transform.c, self.transform.f
        return xmin, ymax + self.transform.e * h, xmin + self.transform.a * w, ymax


def grid_window(bounds, transform: Affine, width: int, height: int) -> Window | None:
    """Alla pixlar som rektangeln berör, begränsat till rutnätet."""
    minx, miny, maxx, maxy = bounds
    inv = ~transform
    c0, r0 = inv @ (minx, maxy)
    c1, r1 = inv @ (maxx, miny)
    col0 = max(int(np.floor(min(c0, c1))), 0)
    row0 = max(int(np.floor(min(r0, r1))), 0)
    col1 = min(int(np.ceil(max(c0, c1))), width)
    row1 = min(int(np.ceil(max(r0, r1))), height)
    if col1 <= col0 or row1 <= row0:
        return None
    return Window(col0, row0, col1 - col0, row1 - row0)


def calibrate(dn: np.ndarray, c0: float, c1: float, c2: float) -> np.ndarray:
    """IRQR-kalibrering c0 + c1*DN + c2*DN^2. Pixlar med DN = 0 förblir 0."""
    dn = dn.astype("float32")
    return np.where(dn > 0, c0 + c1 * dn + c2 * dn * dn, 0.0).astype("float32")


def read_calibration(path: Path) -> dict[tuple[str, int], tuple[float, float, float]]:
    rows = [r for r in csv.reader(l for l in path.open(encoding="utf-8") if not l.startswith("#"))]
    header, body = rows[0], rows[1:]
    idx = {k: header.index(k) for k in ("satellite", "year", "c0", "c1", "c2")}
    return {(r[idx["satellite"]], int(r[idx["year"]])):
            (float(r[idx["c0"]]), float(r[idx["c1"]]), float(r[idx["c2"]])) for r in body}


class DmspIrqr:
    """DMSP-OLS v4 Stable Lights, IRQR-kalibrerade, ett år.

    - Varje satellit kalibreras med c0 + c1*DN + c2*DN^2 (tabell 2 i Hall m.fl. 2019).
      Pixlar med DN = 0 förblir 0: utan kravet skulle varje mörk pixel få värdet
      c0 och få vikt i centroiden.
    - År med två satelliter: medelvärdet av de kalibrerade värdena.
    - Pixlar som berörs av en fackelpolygon sätts till 0.
    """

    name = "dmsp-irqr"

    def __init__(self, year: int, raw_dir: Path, calibration: Path, flares: Path | None):
        coeffs = read_calibration(calibration)
        self.year = year
        self.satellites = sorted(sat for (sat, y) in coeffs if y == year)
        if not self.satellites:
            raise ValueError(f"Inga kalibreringskoefficienter för {year}")
        self.coeffs = [coeffs[(s, year)] for s in self.satellites]
        self.paths = [raw_dir / f"DMSP_{s}{year}_stable_lights.tif" for s in self.satellites]
        for p in self.paths:
            if not p.exists():
                raise FileNotFoundError(p)
        self._srcs = [rasterio.open(p) for p in self.paths]
        for s in self._srcs:
            if (s.width, s.height) != (DMSP_WIDTH, DMSP_HEIGHT):
                raise ValueError(f"{s.name}: {s.width}x{s.height}, väntade {DMSP_WIDTH}x{DMSP_HEIGHT}")
        self.flares = gpd.read_file(flares, engine="pyogrio") if flares else None
        self._flare_index = self.flares.sindex if self.flares is not None else None

    def describe(self) -> dict:
        return {"source": self.name, "year": self.year, "satellites": self.satellites,
                "coefficients": dict(zip(self.satellites, self.coeffs)),
                "files": [str(p) for p in self.paths],
                "grid": "EOG DMSP v4, 30 arcsec, hörn (-180.0041667, 75.0041667)"}

    # Rutnätet som källan levererar fönster ur.
    transform, width, height = DMSP_TRANSFORM, DMSP_WIDTH, DMSP_HEIGHT

    def read(self, bounds) -> Window2D | None:
        win = grid_window(bounds, DMSP_TRANSFORM, DMSP_WIDTH, DMSP_HEIGHT)
        if win is None:
            return None
        return self.read_window(win)

    def read_window(self, win: Window) -> Window2D:
        acc = None
        for src, (c0, c1, c2) in zip(self._srcs, self.coeffs):
            cal = calibrate(src.read(1, window=win), c0, c1, c2)
            acc = cal if acc is None else acc + cal
        data = acc / len(self._srcs)
        transform = DMSP_TRANSFORM * Affine.translation(win.col_off, win.row_off)
        out = Window2D(data=data, transform=transform)
        self._mask_flares(out)
        return out

    def _mask_flares(self, w: Window2D) -> None:
        if self._flare_index is None:
            return
        hits = self._flare_index.query(_box(w.bounds))
        if len(hits) == 0:
            return
        mask = rasterize(((g, 1) for g in self.flares.geometry.iloc[hits]),
                         out_shape=w.data.shape, transform=w.transform,
                         fill=0, all_touched=True, dtype="uint8")
        w.data[mask.astype(bool)] = 0.0

    def close(self) -> None:
        for s in self._srcs:
            s.close()


def _box(bounds):
    from shapely.geometry import box
    return box(*bounds)
