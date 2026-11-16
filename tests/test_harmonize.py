"""Tester av harmoniseringen VIIRS -> DMSP (wpc/harmonize.py)."""
import numpy as np

from wpc import harmonize as H


def test_aggregation_conserves_light_and_centres():
    # 2 DMSP-rader (a..b) -> VIIRS-raderna 2a-1..2b+1 = 5 rader; 7 VIIRS-kolumner
    # motsvarar 4 DMSP-kolumner.
    v = np.zeros((5, 7), "float32")
    v[1, 2] = 4.0                       # jämn rad och kolumn = mitt i DMSP (0, 1)
    a = H.aggregate_to_dmsp(v)
    assert a.shape == (2, 4)
    assert a[0, 1] == 1.0 and a.sum() == 1.0
    v[:] = 0
    v[2, 3] = 4.0                       # udda rad och kolumn = hörn mellan fyra DMSP-celler
    a = H.aggregate_to_dmsp(v)
    assert np.allclose(a[:, 1:3], 0.25) and np.isclose(a.sum(), 1.0)
    assert np.allclose(H.aggregate_to_dmsp(np.ones((5, 7)))[0, 1:3], 1.0)


def test_blur_preserves_sum_inside():
    a = np.zeros((41, 41), "float32"); a[20, 20] = 1
    b = H.blur(a, 2.0)
    assert np.isclose(b.sum(), 1.0, atol=1e-6) and b[20, 20] == b.max()
    assert np.array_equal(H.blur(a, 0), a)


def test_pav_is_monotone_and_weighted():
    y = np.array([1.0, 3.0, 2.0, 4.0]); w = np.array([1, 1, 3, 1])
    out = H.pav(y, w)
    assert np.all(np.diff(out) >= 0)
    assert np.isclose(out[1], 2.25) and np.isclose(out[2], 2.25)


def test_hurdle_fit_gives_zero_background():
    rng = np.random.default_rng(0)
    x = 10 ** rng.uniform(-2.5, 3, 200_000)
    lit = rng.random(x.size) < 1 / (1 + np.exp(-(np.log10(x) - 0.5) * 6))
    y = np.where(lit, np.clip(10 * np.log10(x) + 20, 3, 63), 0.0)
    st = H.BinStats(); st.add(H.bin_index(x), y)
    t = H.fit_table(st)
    assert t["value"][H.ZERO_BIN] == 0
    pred = H.apply_table(np.array([0.0, 1.0, 10.0, 100.0, 1000.0]), t["value"])
    assert pred[0] == 0 and pred[1] == 0            # under tröskeln (p < 1/2): 0
    assert 27 < pred[2] < 33 and pred[3] > pred[2] and pred[4] <= 63.0001
    assert np.all(np.diff(t["value"][1:][t["value"][1:] > 0]) >= -1e-12)
