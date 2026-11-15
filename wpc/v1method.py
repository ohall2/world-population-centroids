"""
Exakt återskapande av v1:s centroidberäkning (Hall m.fl. 2019).

v1 är skriven i R (paketen raster, sp, rgdal). Den här modulen gör samma sak
i Python, steg för steg, så att v1-utdata kan återskapas och avvikelser
förklaras. Förbättringar hör hemma i v2-koden, inte här.

v1:s steg per administrativ enhet (scripts/Barmore.R, scripts/AboufadelAustin.R):
  1. crop(light, extent(area), snap = 'out')
       Beskär rastret till enhetens omslutande rektangel, utåt till hela pixlar.
  2. mask(..., area)
       Pixlar vars mittpunkt ligger utanför polygonen blir NA.
  3. Om summan av ljuset i enheten är 0 får enheten ingen centroid (NA).
  4. rasterToPoints(...)
       Varje kvarvarande pixel blir en punkt i pixelns mittpunkt med ljusvärdet
       som vikt. Pixlar med värdet 0 följer med men väger 0.
  5a. Aboufadel–Austin: ljusviktat medel av 3D-enhetsvektorer på en sfär,
      tillbaka till latitud/longitud.
  5b. Barmore: start i ljusviktat medel av longitud och latitud; därefter
      ljusviktat medel av x och y i en azimutal ekvidistant projektion (sfär,
      radie 6 371 000 m) centrerad i föregående skattning, tills
      abs(x) < 1 m ELLER abs(y) < 1 m (v1 använder `&&` i while-villkoret).

Ingen ytviktning av pixlarna görs i v1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from affine import Affine
from pyproj import CRS, Transformer
from rasterio.features import rasterize
from rasterio.windows import Window
from shapely.geometry.base import BaseGeometry

SPHERE_RADIUS = 6371000.0
GCS_SPHERE = CRS.from_proj4(f"+proj=longlat +a={SPHERE_RADIUS} +b={SPHERE_RADIUS} +no_defs")


@dataclass
class Pixels:
    """Pixlar i en enhet: mittpunkternas koordinater (grader) och ljusvärden."""
    lon: np.ndarray
    lat: np.ndarray
    w: np.ndarray

    @property
    def total(self) -> float:
        return float(self.w.sum())


# --------------------------------------------------------------------------
# Steg 1–4: pixlar i en enhet
# --------------------------------------------------------------------------

def crop_window(bounds: tuple[float, float, float, float], transform: Affine,
                width: int, height: int) -> Window | None:
    """Fönster som motsvarar raster::crop(..., snap = 'out').

    Alla pixlar som rektangeln berör tas med. Returnerar None om rektangeln
    inte överlappar rastret (v1: intersect(extent(area), extent(light)) är NULL).
    """
    minx, miny, maxx, maxy = bounds
    inv = ~transform
    c0, r0 = inv @ (minx, maxy)
    c1, r1 = inv @ (maxx, miny)
    col0, col1 = int(np.floor(min(c0, c1))), int(np.ceil(max(c0, c1)))
    row0, row1 = int(np.floor(min(r0, r1))), int(np.ceil(max(r0, r1)))
    col0, row0 = max(col0, 0), max(row0, 0)
    col1, row1 = min(col1, width), min(row1, height)
    if col1 <= col0 or row1 <= row0:
        return None
    return Window(col0, row0, col1 - col0, row1 - row0)


def unit_pixels(geom: BaseGeometry, light: np.ndarray, transform: Affine) -> Pixels | None:
    """Pixlar vars mittpunkt ligger i polygonen, med sina ljusvärden.

    Returnerar None om enheten ligger utanför rastret.
    """
    height, width = light.shape
    win = crop_window(geom.bounds, transform, width, height)
    if win is None:
        return None
    r0, c0 = int(win.row_off), int(win.col_off)
    h, w = int(win.height), int(win.width)
    win_transform = transform * Affine.translation(c0, r0)
    # GDAL:s regel (all_touched=False) är densamma som raster::mask:
    # pixeln räknas om dess mittpunkt ligger i polygonen.
    inside = rasterize([(geom, 1)], out_shape=(h, w), transform=win_transform,
                       fill=0, all_touched=False, dtype="uint8").astype(bool)
    rows, cols = np.nonzero(inside)
    vals = light[r0 + rows, c0 + cols].astype("float64")
    xs, ys = win_transform @ (cols + 0.5, rows + 0.5)
    return Pixels(lon=np.asarray(xs, dtype="float64"), lat=np.asarray(ys, dtype="float64"), w=vals)


# --------------------------------------------------------------------------
# Steg 5a: Aboufadel–Austin
# --------------------------------------------------------------------------

def aboufadel_austin(px: Pixels) -> tuple[float, float]:
    lon, lat = np.deg2rad(px.lon), np.deg2rad(px.lat)
    wi = px.w.sum()
    x = (px.w * np.cos(lat) * np.cos(lon)).sum() / wi
    y = (px.w * np.cos(lat) * np.sin(lon)).sum() / wi
    z = (px.w * np.sin(lat)).sum() / wi
    a = np.sqrt(x * x + y * y + z * z)
    return float(np.rad2deg(np.arctan2(y / a, x / a))), float(np.rad2deg(np.arcsin(z / a)))


# --------------------------------------------------------------------------
# Steg 5b: Barmore
# --------------------------------------------------------------------------

def _aeqd(lon0: float, lat0: float) -> CRS:
    return CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +x_0=0 +y_0=0 "
        f"+a={SPHERE_RADIUS} +b={SPHERE_RADIUS} +units=m +no_defs")


def barmore(px: Pixels, max_iter: int = 1000, v1_stop_rule: bool = True) -> tuple[float, float, int]:
    """Returnerar (lon, lat, antal iterationer).

    v1_stop_rule=True följer v1 (`&&`): iterationen avbryts så snart en av
    komponenterna är under 1 m. Med False krävs att avståndet är under 1 m.
    """
    w = px.w
    lon0 = float((px.lon * w).sum() / w.sum())
    lat0 = float((px.lat * w).sum() / w.sum())
    crs = _aeqd(lon0, lat0)
    x, y = _project_mean(px, crs)
    n = 0

    def keep_going(x: float, y: float) -> bool:
        if v1_stop_rule:
            return abs(x) >= 1 and abs(y) >= 1
        return np.hypot(x, y) >= 1

    while keep_going(x, y):
        n += 1
        if n > max_iter:
            raise RuntimeError("Barmore konvergerade inte")
        back = Transformer.from_crs(crs, GCS_SPHERE, always_xy=True)
        lon0, lat0 = back.transform(x, y)
        crs = _aeqd(lon0, lat0)
        x, y = _project_mean(px, crs)

    back = Transformer.from_crs(crs, GCS_SPHERE, always_xy=True)
    lon, lat = back.transform(x, y)
    return float(lon), float(lat), n


def _project_mean(px: Pixels, crs: CRS) -> tuple[float, float]:
    fwd = Transformer.from_crs(GCS_SPHERE, crs, always_xy=True)
    xs, ys = fwd.transform(px.lon, px.lat)
    w = px.w
    return float((xs * w).sum() / w.sum()), float((ys * w).sum() / w.sum())


# --------------------------------------------------------------------------
# Hjälpfunktioner
# --------------------------------------------------------------------------

def great_circle_m(lon1, lat1, lon2, lat2) -> np.ndarray:
    """Storcirkelavstånd i meter på v1:s sfär."""
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * SPHERE_RADIUS * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
