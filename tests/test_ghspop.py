"""Test av GhsPop-källan (analys B): totaler per cell -> täthet, så att
centroidkoden ger den befolkningsviktade centroiden.

Visar att centroiden ur GhsPop (täthet x andel x yta) är exakt densamma som en
direkt beräkning med vikten antal x andel, och att den skiljer sig från den
felaktiga vikten antal x andel x yta (som en vanlig radianskälla skulle ge).
"""

import json

import numpy as np
import pytest
import rasterio
from affine import Affine
from rasterio.crs import CRS
from shapely.geometry import box

from wpc import centroids as C
from wpc.lights import make_source

STEM = "GHS_POP_E2000_GLOBE_R2023A_4326_30ss_V1_0"


def _write_raster(tmp_path, counts, transform):
    tif = tmp_path / f"{STEM}.tif"
    with rasterio.open(tif, "w", driver="GTiff", height=counts.shape[0], width=counts.shape[1],
                       count=1, dtype="float32", crs=CRS.from_epsg(4326), transform=transform,
                       nodata=-200.0) as dst:
        dst.write(counts.astype("float32"), 1)
    (tmp_path / f"{STEM}.json").write_text(json.dumps({"year": 2000}))
    return tif


def _counts_and_grid():
    # 0,5-graders celler, 6x6, övre vänstra hörnet (10, 63): latituder 60–63,
    # så cellytan varierar tydligt mellan raderna.
    transform = Affine(0.5, 0.0, 10.0, 0.0, -0.5, 63.0)
    rng = np.random.default_rng(0)
    counts = (rng.random((6, 6)) * 100).astype("float32")
    counts[0, 0] = 0.0          # några tomma celler
    counts[5, 5] = 0.0
    return counts, transform


def _direct_centroid(counts, transform, geom, weight_with_area: bool):
    """Centroid direkt ur cellerna med vikten antal x andel (eller x yta)."""
    h, w = counts.shape
    bounds = (transform.c, transform.f + transform.e * h, transform.c + transform.a * w, transform.f)
    px = C.unit_pixels(geom, counts, bounds, (transform.a, -transform.e))
    weight = px.light * px.frac
    if weight_with_area:
        weight = weight * px.area
    lit = (px.frac > 0) & (px.light > 0)
    aa = C.aboufadel_austin(px.lon[lit], px.lat[lit], weight[lit])
    bm = C.barmore(px.lon[lit], px.lat[lit], weight[lit], aa)
    return aa, bm[:2]


def test_ghspop_density_equals_count_times_fraction(tmp_path):
    counts, transform = _counts_and_grid()
    _write_raster(tmp_path, counts, transform)
    src = make_source("ghspop", 2000, tmp_path, calibration=None, flares=None)
    geom = box(10.75, 60.25, 12.75, 62.25)      # skär genom celler -> andelar < 1

    r = C.compute_unit(geom, src)
    aa_direct, bm_direct = _direct_centroid(counts, transform, geom, weight_with_area=False)

    # lika så långt float32-rutnätet räcker (täthet = antal/yta lagras i float32;
    # 1e-6 grader är omkring 0,1 m)
    assert r.aa_lon == pytest.approx(aa_direct[0], abs=1e-6)
    assert r.aa_lat == pytest.approx(aa_direct[1], abs=1e-6)
    assert r.bm_lon == pytest.approx(bm_direct[0], abs=1e-5)
    assert r.bm_lat == pytest.approx(bm_direct[1], abs=1e-5)


def test_ghspop_differs_from_area_weighted(tmp_path):
    """Täthetsomräkningen spelar roll: vikten antal x andel ger en annan centroid
    än antal x andel x yta (som en råtotal genom centroidkoden skulle ge)."""
    counts, transform = _counts_and_grid()
    _write_raster(tmp_path, counts, transform)
    src = make_source("ghspop", 2000, tmp_path, calibration=None, flares=None)
    geom = box(10.0, 60.0, 13.0, 63.0)

    r = C.compute_unit(geom, src)
    _, bm_area = _direct_centroid(counts, transform, geom, weight_with_area=True)
    assert abs(r.bm_lat - bm_area[1]) > 1e-3     # ytviktningen drar centroiden söderut
