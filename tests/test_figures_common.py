"""Figurernas hjälpfunktioner (tools/figures/common.py) och mellanfilerna."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tools.figures.common import direction_bins, wilson, COL

R = Path(__file__).resolve().parents[1]


def test_direction_bins_known_answer():
    # fyra enheter i [0,1; 0,25): tre rätt, en fel; två i [1, 2): båda rätt
    t = np.array([[0.2, 0], [0, 0.15], [0.12, 0], [0.2, 0], [1.5, 0], [0, -1.2]])
    e = np.array([[0.5, 0.1], [0, 0.4], [0.3, -0.2], [-0.2, 0], [1.0, 0], [0, -3.0]])
    b = direction_bins(t, e).set_index("true_move_km")
    assert list(b.index) == ["[0.1, 0.25)", "[1.0, 2.0)"]
    assert b.loc["[0.1, 0.25)", "n"] == 4 and b.loc["[0.1, 0.25)", "dir_right"] == 0.75
    assert b.loc["[1.0, 2.0)", "dir_right"] == 1.0
    assert b.loc["[1.0, 2.0)", "true_med"] == pytest.approx(1.35)
    assert b.loc["[1.0, 2.0)", "err_med"] == pytest.approx(1.15)


def test_direction_bins_left_closed():
    b = direction_bins(np.array([[0.1, 0.0]]), np.array([[1.0, 0.0]]))
    assert b.true_move_km.tolist() == ["[0.1, 0.25)"]


def test_wilson_known_values():
    lo, hi = wilson(0.5, 100)
    assert lo == pytest.approx(0.4038, abs=1e-4) and hi == pytest.approx(0.5962, abs=1e-4)
    lo, hi = wilson(1.0, 20)
    assert hi == pytest.approx(1.0) and lo == pytest.approx(0.8389, abs=1e-4)


def test_colours_distinct():
    assert len(set(COL.values())) == len(COL)
    assert COL["census"] == "#2a78d6" and COL["inegi"] == "#eb6834" and COL["geostat"] == "#1baf7a"


SH = R / "docs/figures/data/unit_shifts.csv.gz"


@pytest.mark.skipif(not SH.exists(), reason="mellanfil saknas")
def test_bins_reproduce_frozen_change_detection():
    """Klasserna ur mellanfilen återger docs/change_detection (census nivå 2)."""
    d = pd.read_csv(SH)
    fz = pd.read_csv(R / "docs/change_detection/census_change_detection_bins.csv")
    for a, b in [(2000, 2010), (2010, 2020), (2000, 2020)]:
        s = d[(d.dataset == "census") & (d["from"] == a) & (d.to == b)]
        got = direction_bins(s[["true_e", "true_n"]].values, s[["meas_e", "meas_n"]].values).set_index("true_move_km")
        ref = fz[(fz.series == "release") & (fz.level == 2) & (fz["from"] == a) & (fz.to == b)].set_index("true_move_km")
        assert list(got.index) == list(ref.index)
        # mellanfilen är avrundad till 1e-4 km; enstaka enheter kan byta klass eller tecken
        assert (got.n - ref.n).abs().max() <= 2
        assert np.allclose(got.dir_right, ref.dir_right, atol=0.01)


@pytest.mark.skipif(not SH.exists(), reason="mellanfil saknas")
def test_shifts_reproduce_confirmatory_counts():
    d = pd.read_csv(SH)
    for ds, a, b, tag in [("statcan", 2016, 2021, "canada"), ("estat_jp", 2015, 2020, "japan")]:
        s = d[(d.dataset == ds) & (d["from"] == a) & (d.to == b)]
        ref = pd.read_csv(R / f"docs/extended/C_{tag}_main.csv").iloc[0]
        right = ((s.true_e*s.meas_e + s.true_n*s.meas_n) > 0).mean()
        assert len(s) == ref.n and right == pytest.approx(ref.dir_right, abs=5e-4)
