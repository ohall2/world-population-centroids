"""Tester av centroidberäkningen med konstgjorda raster där svaret är känt."""

import numpy as np
import pytest
from affine import Affine
from shapely.geometry import Polygon, box

from wpc import centroids as C
from wpc.lights import Window2D, calibrate, grid_window
from wpc.run import flags


class FakeSource:
    """Globalt rutnät med 1-gradspixlar, (-180, 90) i övre vänstra hörnet."""

    transform = Affine(1.0, 0.0, -180.0, 0.0, -1.0, 90.0)
    width, height = 360, 180

    def __init__(self):
        self.grid = np.zeros((self.height, self.width), dtype="float32")

    def light(self, lon_center, lat_center, value):
        col, row = ~self.transform @ (lon_center, lat_center)
        self.grid[int(row), int(col)] = value

    def read_window(self, win):
        r0, c0 = int(win.row_off), int(win.col_off)
        data = self.grid[r0:r0 + int(win.height), c0:c0 + int(win.width)].copy()
        return Window2D(data=data, transform=self.transform @ Affine.translation(c0, r0))

    def read(self, bounds):
        win = grid_window(bounds, self.transform, self.width, self.height)
        return None if win is None else self.read_window(win)


def test_single_pixel_gives_pixel_centre():
    s = FakeSource()
    s.light(10.5, 45.5, 20)
    r = C.compute_unit(box(5, 40, 15, 50), s)
    assert r.aa_lon == pytest.approx(10.5) and r.aa_lat == pytest.approx(45.5)
    assert r.bm_lon == pytest.approx(10.5) and r.bm_lat == pytest.approx(45.5)
    assert r.n_lit == 1 and r.n_cells == 100


def test_symmetric_pixels_on_equator():
    s = FakeSource()
    s.light(-1.5, 0.5, 10)
    s.light(1.5, 0.5, 10)
    r = C.compute_unit(box(-5, 0, 5, 1), s)
    assert r.aa_lon == pytest.approx(0, abs=1e-9)
    assert r.bm_lon == pytest.approx(0, abs=1e-6)


def test_antimeridian_unit_stays_at_180():
    """En enhet med ljus på båda sidor om datumlinjen får centroiden vid 180°, inte 0°."""
    s = FakeSource()
    s.light(179.5, -17.5, 10)
    s.light(-179.5, -17.5, 10)
    geom = box(178, -18, 180, -17).union(box(-180, -18, -178, -17))
    r = C.compute_unit(geom, s)
    assert abs(abs(r.aa_lon) - 180) < 1e-6
    assert abs(abs(r.bm_lon) - 180) < 1e-3
    assert r.aa_lat == pytest.approx(-17.5, abs=1e-3)


def test_area_weighting_pulls_towards_larger_pixels():
    """Lika ljusa pixlar vid 0° och 60°: pixeln vid 60° är ungefär hälften så stor."""
    s = FakeSource()
    s.light(0.5, 0.5, 10)
    s.light(0.5, 60.5, 10)
    r = C.compute_unit(box(0, 0, 1, 61), s)
    a0 = C.cell_area_km2(np.array([0.5]), 1, 1)[0]
    a60 = C.cell_area_km2(np.array([60.5]), 1, 1)[0]
    assert a60 / a0 == pytest.approx(np.cos(np.deg2rad(60.5)), rel=1e-3)
    expected = C.aboufadel_austin(np.array([0.5, 0.5]), np.array([0.5, 60.5]), np.array([a0, a60]))
    assert r.aa_lat == pytest.approx(expected[1])
    assert r.aa_lat < 30.5  # mitten utan ytviktning


def test_partial_pixel_coverage():
    """En polygon som täcker halva den enda ljusa pixeln får hälften av ljuset."""
    s = FakeSource()
    s.light(10.5, 45.5, 20)
    r = C.compute_unit(box(10, 45, 10.5, 46), s)
    assert r.sum_light == pytest.approx(10.0)
    assert r.aa_lon == pytest.approx(10.5)  # pixelns mittpunkt används som läge


def test_tiling_gives_identical_result(monkeypatch):
    rng = np.random.default_rng(1)
    s = FakeSource()
    s.grid[30:80, 100:200] = rng.random((50, 100)).astype("float32") * 63
    geom = Polygon([(-80, 60), (20, 55), (15, 10), (-75, 0)])
    whole = C.compute_unit(geom, s)
    monkeypatch.setattr(C, "MAX_CELLS", 100)
    monkeypatch.setattr(C, "TILE", 7)
    tiled = C.compute_unit(geom, s)
    assert tiled.n_cells == whole.n_cells
    assert tiled.sum_weight == pytest.approx(whole.sum_weight, rel=1e-9)
    assert C.great_circle_m(whole.aa_lon, whole.aa_lat, tiled.aa_lon, tiled.aa_lat) < 1e-3


def test_unit_without_light_has_no_centroid():
    r = C.compute_unit(box(0, 0, 5, 5), FakeSource())
    assert np.isnan(r.aa_lon) and r.n_lit == 0


def test_unit_outside_grid():
    s = FakeSource()
    s.transform = Affine(1.0, 0.0, -180.0, 0.0, -1.0, 75.0)
    s.height = 140
    r = C.compute_unit(box(0, -89, 5, -80), s)
    assert np.isnan(r.aa_lon) and r.n_cells == 0


def test_calibration_keeps_zero_and_follows_formula():
    dn = np.array([0, 1, 10, 63], dtype="uint8")
    out = calibrate(dn, 0.149, 1.8118, -0.015)
    assert out[0] == 0
    assert out[2] == pytest.approx(0.149 + 18.118 - 1.5, rel=1e-6)


def test_flags():
    import geopandas as gpd
    units = gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1), box(1, 0, 2, 1), box(5, 5, 6, 6)])
    lon = np.array([0.5, 0.5, 10.0])     # i egen, i annan enhet, utanför alla
    lat = np.array([0.5, 0.5, 10.0])
    assert flags(units, lon, lat).tolist() == [0, 2, 1]
    assert flags(units, np.array([np.nan, 1.5, np.nan]), np.array([np.nan, 0.5, np.nan])).tolist() == [-1, 0, -1]
