"""Bevis att Kanada-drivrutinen (tools/canada_confirm.py) räknar som de frysta
skripten: kört på ett befintligt facit (census nivå 2, 2010–2020) ska den
återge de committade talen i docs/change_detection/census_change_detection.csv
exakt (n = 2383, rätt riktning 0,5820, median |e−t| 0,7947). Stämmer det inte
får Kanada-resultaten inte lita på."""

from pathlib import Path

import pytest

from tools.canada_confirm import change_measure

RELEASE = Path("data/release/v2.0.0")
TRUTH = Path("data/validation")


@pytest.mark.skipif(not (TRUTH / "census" / "2010" / "adm2.csv").exists(),
                    reason="facit/utgåva saknas (data/ ingår inte i git)")
def test_driver_matches_frozen_change_detection_census_l2():
    r = change_measure(TRUTH, RELEASE, "census", 2010, 2020, level=2)
    assert r["n"] == 2383
    assert r["metrics"]["dir_right"] == pytest.approx(0.5820, abs=5e-5)
    assert r["metrics"]["err_med"] == pytest.approx(0.7947, abs=5e-5)
