"""Tester av A2-trendskattningen med konstgjorda data där svaret är känt.

Skrivna före analysskriptet (TDD). Måtten (metrics, en, block_ids, bootstrap)
och ändpunkts-/treårsskattningarna återanvänds oförändrade ur
tools/change_detection.py och tools/robustness_checks.py och prövas inte här.
Dessa tester gäller det som är nytt i A2:

  1. Theil–Sen-lutning (median av parvisa lutningar, robust mot en avvikare)
  2. periodindelning kring sensorskarven (DMSP t.o.m. 2012, VIIRS fr.o.m. 2013)
  3. flytt ur lutningarna, inklusive kravet på minst fyra år i varje del
"""

import numpy as np
import pytest

from tools import trend_estimate as T


# 1. Theil–Sen-lutning --------------------------------------------------------

def test_theil_sen_exact_line():
    years = np.arange(2000, 2013)
    assert T.theil_sen_slope(years, 0.5 * years - 3.0) == pytest.approx(0.5)
    assert T.theil_sen_slope(years, -0.2 * years + 10.0) == pytest.approx(-0.2)


def test_theil_sen_robust_to_outlier():
    years = np.arange(2000, 2013)
    y = 0.5 * years
    y[6] += 1000.0   # en grov avvikare mitt i serien
    assert T.theil_sen_slope(years, y) == pytest.approx(0.5)


# 2. Periodindelning kring skarven -------------------------------------------

def test_trend_spec_within_dmsp():
    # 2000-2010 ligger helt inom DMSP: ett fönster y1-2 .. y2+2, högst till 2012
    assert T.trend_spec(2000, 2010) == [(1998, 2012, 10)]


def test_trend_spec_across_seam():
    # 2010-2020 korsar skarven: DMSP y1-2..2012, VIIRS 2013..y2+2 (högst 2025)
    assert T.trend_spec(2010, 2020) == [(2008, 2012, 3), (2013, 2022, 7)]
    assert T.trend_spec(2000, 2020) == [(1998, 2012, 13), (2013, 2022, 7)]
    assert T.trend_spec(2011, 2021) == [(2009, 2012, 2), (2013, 2023, 8)]


# 3. Flytt ur lutningarna -----------------------------------------------------

def test_movement_within_dmsp_known():
    """Linjär serie inom DMSP: flytt = lutning × (y2 − y1)."""
    years = np.arange(1998, 2013)              # 15 år, alla i fönstret
    east = 0.5 * (years - 2005)                # lutning 0,5 km/år
    north = -0.2 * (years - 2005)              # lutning −0,2 km/år
    me, mn = T.trend_movement_series(years, east, north, 2000, 2010)
    assert me == pytest.approx(0.5 * 10)       # 5,0 km
    assert mn == pytest.approx(-0.2 * 10)       # −2,0 km


def test_movement_across_seam_known():
    """Olika lutning i varje del; flytt = s_D×(2013−y1) + s_V×(y2−2013)."""
    years = np.arange(2008, 2023)              # täcker båda delarna
    # öst: 0,3 km/år i DMSP-delen, 0,1 km/år i VIIRS-delen (konstanterna spelar ingen roll)
    east = np.where(years <= 2012, 0.3 * (years - 2010), 0.1 * (years - 2013) + 5.0)
    north = np.where(years <= 2012, -0.4 * (years - 2010), 0.2 * (years - 2013) - 1.0)
    me, mn = T.trend_movement_series(years, east, north, 2010, 2020)
    assert me == pytest.approx(0.3 * 3 + 0.1 * 7)    # 1,6 km
    assert mn == pytest.approx(-0.4 * 3 + 0.2 * 7)   # 0,2 km


def test_movement_needs_four_years_each_part():
    """Färre än fyra år i en del gör att enheten utgår (NaN)."""
    # DMSP-delen har bara tre år (2010, 2011, 2012)
    years = np.array([2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017])
    east = 0.3 * (years - 2013)
    north = 0.0 * years
    me, mn = T.trend_movement_series(years, east, north, 2010, 2020)
    assert np.isnan(me) and np.isnan(mn)


def test_movement_handles_missing_years():
    """NaN-år hoppas över; räcker de kvarvarande till fyra per del räknas flytten."""
    years = np.arange(1998, 2013)
    east = 0.5 * (years - 2005)
    north = np.zeros_like(east, dtype=float)
    east[3] = np.nan                           # ett saknat år, 14 kvar >= 4
    me, mn = T.trend_movement_series(years, east, north, 2000, 2010)
    assert me == pytest.approx(0.5 * 10)
    assert mn == pytest.approx(0.0)
