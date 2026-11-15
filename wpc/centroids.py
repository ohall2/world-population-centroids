"""Centroidberäkning för v2.

Per administrativ enhet:
  1. Ljusfönstret över enhetens omslutande rektangel hämtas från ljuskällan.
  2. exactextract ger, för varje pixel som enheten berör, hur stor andel av
     pixeln som ligger i enheten (0–1) samt pixelns mittpunkt.
  3. Vikten för en pixel är  ljus x andel x pixelns yta på sfären.
     - Andelen ersätter v1:s regel att bara räkna pixlar vars mittpunkt ligger
       i polygonen. Små enheter får därmed alltid sina pixlar med.
     - Ytan behövs eftersom en pixel på 30" är mindre längre från ekvatorn;
       ljusvärdet är en intensitet per yta.
  4. Aboufadel–Austin: viktat medel av 3D-enhetsvektorer på en sfär.
     Barmore: viktat medel av x, y i en azimutal ekvidistant projektion (sfär),
     upprepat tills förflyttningen är under 1 m. Startpunkten är
     Aboufadel–Austin-centroiden (v1 startade i medel av longitud och latitud,
     vilket blir fel för enheter som korsar datumlinjen).

Pixlar som ligger delvis i flera enheter delas proportionellt mellan dem.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import shapely
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource
from pyproj import CRS, Transformer
from rasterio.windows import Window
from shapely.geometry.base import BaseGeometry

from .lights import grid_window

SPHERE_RADIUS = 6371000.0
GCS_SPHERE = CRS.from_proj4(f"+proj=longlat +a={SPHERE_RADIUS} +b={SPHERE_RADIUS} +no_defs")


@dataclass
class UnitPixels:
    lon: np.ndarray       # pixelns mittpunkt, grader
    lat: np.ndarray
    light: np.ndarray     # ljusvärde
    frac: np.ndarray      # andel av pixeln i enheten
    area: np.ndarray      # pixelns yta, km²

    @property
    def weight(self) -> np.ndarray:
        return self.light * self.frac * self.area


@dataclass
class UnitResult:
    aa_lon: float = np.nan
    aa_lat: float = np.nan
    bm_lon: float = np.nan
    bm_lat: float = np.nan
    bm_iter: int = 0
    aa_bm_m: float = np.nan
    sum_light: float = 0.0     # Σ ljus x andel (jämförbart med v1:s summa)
    sum_weight: float = 0.0    # Σ ljus x andel x yta
    n_cells: int = 0           # pixlar som enheten berör
    n_lit: int = 0             # varav med ljus > 0
    area_km2: float = 0.0      # enhetens yta enligt pixelandelarna


def cell_area_km2(lat_center: np.ndarray, dlat: float, dlon: float) -> np.ndarray:
    """Pixelns yta på sfären (km²) för pixlar med mittpunkt på latitud lat_center."""
    r = SPHERE_RADIUS / 1000
    north = np.deg2rad(np.clip(lat_center + dlat / 2, -90, 90))
    south = np.deg2rad(np.clip(lat_center - dlat / 2, -90, 90))
    return r * r * np.deg2rad(dlon) * (np.sin(north) - np.sin(south))


def unit_pixels(geom: BaseGeometry, data: np.ndarray, bounds, res: tuple[float, float]) -> UnitPixels:
    """Pixlar som enheten berör, med andel, mittpunkt, ljus och yta."""
    xmin, ymin, xmax, ymax = bounds
    rast = NumPyRasterSource(data, xmin, ymin, xmax, ymax)
    feat = {"type": "Feature", "properties": {}, "geometry": geom.__geo_interface__}
    out = exact_extract(rast, [feat], ["coverage", "values", "center_x", "center_y"], output="geojson")
    p = out[0]["properties"]
    lon = np.asarray(p["center_x"], dtype="float64")
    lat = np.asarray(p["center_y"], dtype="float64")
    return UnitPixels(lon=lon, lat=lat,
                      light=np.asarray(p["values"], dtype="float64"),
                      frac=np.asarray(p["coverage"], dtype="float64"),
                      area=cell_area_km2(lat, res[1], res[0]))


def aboufadel_austin(lon, lat, w) -> tuple[float, float]:
    lo, la = np.deg2rad(lon), np.deg2rad(lat)
    cl = np.cos(la)
    x, y, z = (w * cl * np.cos(lo)).sum(), (w * cl * np.sin(lo)).sum(), (w * np.sin(la)).sum()
    n = np.sqrt(x * x + y * y + z * z)
    return float(np.rad2deg(np.arctan2(y, x))), float(np.rad2deg(np.arcsin(z / n)))


def _aeqd(lon0: float, lat0: float) -> CRS:
    return CRS.from_proj4(f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +x_0=0 +y_0=0 "
                          f"+a={SPHERE_RADIUS} +b={SPHERE_RADIUS} +units=m +no_defs")


def barmore(lon, lat, w, start: tuple[float, float], tol_m: float = 1.0,
            max_iter: int = 100) -> tuple[float, float, int]:
    lon0, lat0 = start
    for n in range(1, max_iter + 1):
        crs = _aeqd(lon0, lat0)
        xs, ys = Transformer.from_crs(GCS_SPHERE, crs, always_xy=True).transform(lon, lat)
        x, y = float((xs * w).sum() / w.sum()), float((ys * w).sum() / w.sum())
        lon0, lat0 = Transformer.from_crs(crs, GCS_SPHERE, always_xy=True).transform(x, y)
        if np.hypot(x, y) < tol_m:
            return float(lon0), float(lat0), n
    raise RuntimeError("Barmore konvergerade inte inom max_iter")


def great_circle_m(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * SPHERE_RADIUS * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


MAX_CELLS = 4_000_000   # större fönster delas i rutor
TILE = 2048              # rutstorlek i pixlar
FRAC_EPS = 1e-9          # andelar under detta är avrundningsfel längs rutgränser


def _tiles(win: Window) -> list[Window]:
    """Dela ett fönster i rutor som följer rutnätets pixelgränser."""
    if win.width * win.height <= MAX_CELLS:
        return [win]
    out = []
    for r in range(int(win.row_off), int(win.row_off + win.height), TILE):
        for c in range(int(win.col_off), int(win.col_off + win.width), TILE):
            out.append(Window(c, r, min(TILE, win.col_off + win.width - c),
                              min(TILE, win.row_off + win.height - r)))
    return out


def compute_unit(geom: BaseGeometry, source) -> UnitResult:
    """Centroider för en enhet. Enheter utan ljus får NaN-koordinater.

    Stora enheter räknas ruta för ruta. Rutorna följer pixelgränserna, så varje
    pixel hör till exakt en ruta och andelen i enheten blir densamma som om
    hela enheten räknats på en gång. Bara pixlar med ljus sparas mellan rutorna.
    """
    res = UnitResult()
    if geom is None or geom.is_empty:
        return res
    full = grid_window(geom.bounds, source.transform, source.width, source.height)
    if full is None:
        return res
    lons, lats, ws = [], [], []
    for win in _tiles(full):
        w2 = source.read_window(win)
        xmin, ymin, xmax, ymax = w2.bounds
        part = shapely.clip_by_rect(geom, xmin, ymin, xmax, ymax) if win is not full else geom
        if part.is_empty:
            continue
        t = w2.transform
        px = unit_pixels(part, w2.data, w2.bounds, (t.a, -t.e))
        px.frac = np.where(px.frac > FRAC_EPS, px.frac, 0.0)
        keep = px.frac > 0
        w = px.weight
        lit = keep & (px.light > 0)
        res.n_cells += int(keep.sum())
        res.n_lit += int(lit.sum())
        res.area_km2 += float((px.frac * px.area).sum())
        res.sum_light += float((px.light * px.frac).sum())
        res.sum_weight += float(w.sum())
        lons.append(px.lon[lit]); lats.append(px.lat[lit]); ws.append(w[lit])
    if res.sum_weight <= 0:
        return res
    lon, lat, w = np.concatenate(lons), np.concatenate(lats), np.concatenate(ws)
    res.aa_lon, res.aa_lat = aboufadel_austin(lon, lat, w)
    res.bm_lon, res.bm_lat, res.bm_iter = barmore(lon, lat, w, (res.aa_lon, res.aa_lat))
    res.aa_bm_m = float(great_circle_m(res.aa_lon, res.aa_lat, res.bm_lon, res.bm_lat))
    return res
