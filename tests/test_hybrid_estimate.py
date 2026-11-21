"""Tester av A1-hybridskattningen med konstgjorda data där svaret är känt.

Skrivna före analysskriptet (TDD). Alla förväntade värden är uträknade för hand
eller följer av symmetri. Måtten för förändring (metrics, en, block_ids,
bootstrap) återanvänds oförändrade ur tools/change_detection.py och prövas
därför inte här; dessa tester gäller det som är nytt i A1:

  1. nivå 2 -> nivå 1-förälder ur GADM-id
  2. befolkningsviktat 3D-medel (hybrid-ljus och hybrid-geometri)
  3. geometrisk ytcentroid i ytriktig projektion (USA cos-phi; övriga LAEA)
  4. täckningsregeln 99 %
  5. aggregeringen per nivå 1-enhet: hybridcentroider, täckning och urval
  6. great circle-avstånd (lägesfelet)
"""

import numpy as np
import pytest
from shapely.geometry import box

import geopandas as gpd

from tools import hybrid_estimate as H


# 1. Förälder på nivå 1 ur GADM-id på nivå 2 -------------------------------

def test_parent_adm1():
    assert H.parent_adm1("USA.5.13_1") == "USA.5_1"
    assert H.parent_adm1("USA.5.13_2") == "USA.5_1"   # delnummer ignoreras
    assert H.parent_adm1("BRA.11.102_1") == "BRA.11_1"
    assert H.parent_adm1("MEX.9.3_1") == "MEX.9_1"


# 2. Befolkningsviktat 3D-medel --------------------------------------------

def test_mean_3d_identical_points():
    """Alla punkter i samma läge ger det läget."""
    lon = np.array([12.3, 12.3, 12.3])
    lat = np.array([-4.5, -4.5, -4.5])
    w = np.array([1.0, 7.0, 2.0])
    glon, glat = H.mean_centroid_3d(lon, lat, w)
    assert glon == pytest.approx(12.3, abs=1e-9)
    assert glat == pytest.approx(-4.5, abs=1e-9)


def test_mean_3d_symmetric_gives_zero():
    """Två lika tunga punkter speglade i (0,0) ger (0,0)."""
    glon, glat = H.mean_centroid_3d(
        np.array([-1.0, 1.0]), np.array([30.0, -30.0]), np.array([5.0, 5.0]))
    assert glon == pytest.approx(0.0, abs=1e-9)
    assert glat == pytest.approx(0.0, abs=1e-9)


def test_mean_3d_weighted_known_angle():
    """(0,0) med vikt 3 och (90,0) med vikt 1: vektorsumma (3,1,0) på ekvatorn,
    longitud = atan2(1, 3)."""
    glon, glat = H.mean_centroid_3d(
        np.array([0.0, 90.0]), np.array([0.0, 0.0]), np.array([3.0, 1.0]))
    assert glon == pytest.approx(np.degrees(np.arctan2(1.0, 3.0)), abs=1e-9)
    assert glat == pytest.approx(0.0, abs=1e-9)


# 3. Geometrisk ytcentroid --------------------------------------------------

def test_geometric_centroid_usa_cos_scaling_is_box_centre():
    """USA-grenen (longitud skalad med cos phi per enhet): en rektangel får sin
    mittpunkt som ytcentroid, exakt."""
    gdf = gpd.GeoDataFrame(
        {"unit_id": ["USA.1.1_1"]},
        geometry=[box(-100.5, 39.0, -99.5, 41.0)], crs=4326)
    glon, glat = H.geometric_centroids(gdf, "USA")
    assert glon[0] == pytest.approx(-100.0, abs=1e-9)
    assert glat[0] == pytest.approx(40.0, abs=1e-9)


def test_geometric_centroid_non_usa_laea_is_box_centre():
    """Övriga länder (LAEA kring landets mittpunkt): en liten ruta får sin
    mittpunkt som ytcentroid. Longitud blir exakt mittpunkten (öst-väst-symmetri
    kring centralmeridianen); latitud ligger en aning söder om mittpunkten,
    eftersom en lon/lat-ruta täcker större yta mot ekvatorn (här ~1e-5 grader
    för en 0,1-gradersruta, växande med rutans storlek)."""
    gdf = gpd.GeoDataFrame(
        {"unit_id": ["SWE.1.1_1"]},
        geometry=[box(14.95, 56.95, 15.05, 57.05)], crs=4326)
    glon, glat = H.geometric_centroids(gdf, "SWE")
    assert glon[0] == pytest.approx(15.0, abs=1e-9)
    assert glat[0] == pytest.approx(57.0, abs=1e-4)
    assert glat[0] < 57.0   # förskjuten mot ekvatorn


# 4. Täckningsregeln --------------------------------------------------------

def test_covered_fraction():
    pop = np.array([60.0, 39.0, 1.0])
    assert H.covered_fraction(pop, np.array([True, True, False])) == pytest.approx(0.99)
    assert H.covered_fraction(pop, np.array([True, True, True])) == pytest.approx(1.0)
    assert H.covered_fraction(pop, np.array([True, False, False])) == pytest.approx(0.60)


# 5. Aggregering per nivå 1-enhet ------------------------------------------

def _adm2_frame():
    """Tre nivå 1-enheter med kända svar:
      AAA.1_1: två nivå 2-enheter, båda med ljus, pop 3 och 1,
               ljuscentroider (0,0) och (90,0) -> hybrid-ljus lon=atan2(1,3);
               geometriska centroider båda (10,0) -> hybrid-geometri (10,0).
               täckning 1.0, används.
      BBB.1_1: pop 60/39/1, enheten med pop 1 saknar ljus -> täckning 0.99, används.
      CCC.1_1: pop 60/38/2, enheten med pop 2 saknar ljus -> täckning 0.98, används ej.
    """
    rows = [
        ("AAA.1.1_1", "AAA.1_1", 3.0, True,  0.0,  0.0, 10.0, 0.0),
        ("AAA.1.2_1", "AAA.1_1", 1.0, True, 90.0,  0.0, 10.0, 0.0),
        ("BBB.1.1_1", "BBB.1_1", 60.0, True, 20.0, 5.0, 20.0, 5.0),
        ("BBB.1.2_1", "BBB.1_1", 39.0, True, 20.0, 5.0, 20.0, 5.0),
        ("BBB.1.3_1", "BBB.1_1",  1.0, False, np.nan, np.nan, 25.0, 9.0),
        ("CCC.1.1_1", "CCC.1_1", 60.0, True, -5.0, 1.0, -5.0, 1.0),
        ("CCC.1.2_1", "CCC.1_1", 38.0, True, -5.0, 1.0, -5.0, 1.0),
        ("CCC.1.3_1", "CCC.1_1",  2.0, False, np.nan, np.nan, -8.0, 2.0),
    ]
    import pandas as pd
    return pd.DataFrame(rows, columns=[
        "unit_id", "adm1", "pop", "has_light", "L_lon", "L_lat", "g_lon", "g_lat"])


def test_hybrid_level1_aggregation():
    out = H.hybrid_level1(_adm2_frame()).set_index("adm1")

    # AAA: känt vinkelmedel för ljus, identiska geometripunkter
    assert out.loc["AAA.1_1", "cover"] == pytest.approx(1.0)
    assert bool(out.loc["AAA.1_1", "used"]) is True
    assert out.loc["AAA.1_1", "hyl_lon"] == pytest.approx(np.degrees(np.arctan2(1.0, 3.0)), abs=1e-9)
    assert out.loc["AAA.1_1", "hyl_lat"] == pytest.approx(0.0, abs=1e-9)
    assert out.loc["AAA.1_1", "hyg_lon"] == pytest.approx(10.0, abs=1e-9)
    assert out.loc["AAA.1_1", "hyg_lat"] == pytest.approx(0.0, abs=1e-9)

    # BBB: precis på gränsen, används; hybrid räknas på enheterna med ljus
    assert out.loc["BBB.1_1", "cover"] == pytest.approx(0.99)
    assert bool(out.loc["BBB.1_1", "used"]) is True
    assert out.loc["BBB.1_1", "hyl_lon"] == pytest.approx(20.0, abs=1e-9)
    assert out.loc["BBB.1_1", "hyl_lat"] == pytest.approx(5.0, abs=1e-9)

    # CCC: under gränsen, används ej
    assert out.loc["CCC.1_1", "cover"] == pytest.approx(0.98)
    assert bool(out.loc["CCC.1_1", "used"]) is False


# 6. Lägesfel: great circle-avstånd ----------------------------------------

def test_great_circle_km_one_degree_latitude():
    """En grad latitud på en sfär med R = 6371 km är R * (pi/180) km."""
    d = H.great_circle_km(0.0, 0.0, 0.0, 1.0)
    assert d == pytest.approx(6371.0 * np.pi / 180.0, rel=1e-9)


def test_great_circle_km_zero_distance():
    assert H.great_circle_km(10.0, 20.0, 10.0, 20.0) == pytest.approx(0.0, abs=1e-9)
