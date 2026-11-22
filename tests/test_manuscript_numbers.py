"""Tester för tools/manuscript_numbers.py (talrevisionen 2026-11-22).

Syntetiska data med känt svar: avrundning halvt uppåt, intervallformat,
andel närmare, Fishers exakta test, OLS-lutning utan de mest inflytelserika
enheterna, exakt Theil–Sen, överensstämmelse mellan två flyttskattningar och
statuslogiken för ett påstående."""

import numpy as np
import pytest
from scipy import stats

import tools.manuscript_numbers as mn


# --- avrundning och format --------------------------------------------------

@pytest.mark.parametrize("x,dec,want", [
    (0.625, 2, "0.63"),       # exakt halva: uppåt (inte till jämnt)
    (0.2475, 2, "0.25"),
    (0.4750, 2, "0.48"),
    (2.5, 0, "3"),
    (0.0558, 2, "0.06"),
    (-0.0848, 2, "−0.08"),    # minus som i manuset (U+2212)
    (-0.015, 2, "−0.02"),     # halva bort från noll
    (1.45, 1, "1.5"),         # 1.45 är 1.4499… som binärt tal; Decimal(str) ger 1.45
    (0.0, 2, "0.00"),
])
def test_round_half_up(x, dec, want):
    assert mn.fmt(x, dec) == want


def test_fmt_scale_percent():
    assert mn.fmt(0.7188, 0, 100) == "72"
    assert mn.fmt(0.875, 0, 100) == "88"
    assert mn.fmt(0.1875, 0, 100) == "19"


def test_fmt_range():
    assert mn.fmt_range([0.67, 2.57, 1.5], 1) == "0.7–2.6"
    assert mn.fmt_range([0.501, 0.499], 2) == "0.50"          # lika efter avrundning
    assert mn.fmt_range([-0.16, -0.01], 2) == "−0.16 to −0.01"  # negativt: "to"
    assert mn.fmt_range([1234.0, 2441.0], 0) == "1,234–2,441"   # tusentalsavgränsare


# --- andel närmare -------------------------------------------------------------

def test_share_closer_known_and_nan():
    a = np.array([1.0, 2.0, 3.0, np.nan, 5.0])
    b = np.array([2.0, 1.0, 4.0, 1.0, 5.0])
    # giltiga par: (1,2) ja, (2,1) nej, (3,4) ja, (5,5) nej (lika är inte närmare)
    share, n = mn.share_closer(a, b)
    assert n == 4
    assert share == pytest.approx(0.5)


# --- Fishers exakta test ---------------------------------------------------------

def test_fisher_known_value():
    # Fishers teprovning: [[3, 1], [1, 3]] -> tvåsidigt p = 34/70
    assert mn.fisher_two_sided(3, 4, 1, 4) == pytest.approx(34 / 70)


def test_fisher_matches_scipy():
    p = mn.fisher_two_sided(23, 32, 25, 50)
    assert p == pytest.approx(stats.fisher_exact([[23, 9], [25, 25]])[1])


# --- OLS-lutning och inflytande ------------------------------------------------

def _clean(n=200, slope=0.5, seed=1):
    rng = np.random.default_rng(seed)
    t = rng.normal(0, 1, (n, 2))
    e = slope * t + rng.normal(0, 0.1, (n, 2))
    return t, e


def test_pooled_slope_matches_frozen_metric():
    from tools.change_detection import metrics
    t, e = _clean()
    assert mn.pooled_slope(t, e) == pytest.approx(metrics(t, e)["slope"], rel=1e-12)


def test_unit_influence_equals_brute_force():
    t, e = _clean(n=40)
    full = mn.pooled_slope(t, e)
    infl = mn.unit_influence(t, e)
    for i in range(len(t)):
        keep = np.arange(len(t)) != i
        assert infl[i] == pytest.approx(mn.pooled_slope(t[keep], e[keep]) - full, abs=1e-12)


def test_slope_without_top_removes_planted_outliers():
    t, e = _clean(n=300, slope=0.5)
    clean_slope = mn.pooled_slope(t, e)
    # tre enheter med stor verklig flytt och kraftigt överskattad uppmätt flytt
    t_out = np.array([[8.0, 0.0], [0.0, 9.0], [-7.0, -7.0]])
    e_out = 4.0 * t_out
    T, E = np.vstack([t, t_out]), np.vstack([e, e_out])
    assert mn.pooled_slope(T, E) > 1.0
    slope, removed = mn.slope_without_top(T, E, k=3)
    assert sorted(removed.tolist()) == [300, 301, 302]
    assert slope == pytest.approx(clean_slope, rel=1e-12)


# --- Theil–Sen -------------------------------------------------------------------

def test_theil_sen_exact_line():
    x = np.arange(10.0)
    assert mn.theil_sen(x, 2.0 * x + 1.0) == pytest.approx(2.0)


def test_theil_sen_matches_scipy_with_ties_and_outliers():
    rng = np.random.default_rng(7)
    x = np.round(rng.normal(0, 1, 401), 1)       # många lika x (utgår ur paren)
    y = 0.3 * x + rng.normal(0, 0.2, 401)
    y[:15] += 50.0                                # avvikare
    want = stats.theilslopes(y, x)[0]
    assert mn.theil_sen(x, y) == pytest.approx(want, rel=1e-12)
    # chunkad tvåstegsberäkning ska ge samma exakta median
    assert mn.theil_sen(x, y, chunk=37, sample=500) == pytest.approx(want, rel=1e-12)


def test_theil_sen_even_number_of_pairs():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 4.0, 3.0])
    # parvisa lutningar: 1, 2, 1, 3, 1, -1 -> sorterade -1,1,1,1,2,3 -> median 1
    assert mn.theil_sen(x, y) == pytest.approx(1.0)


# --- två skattningar av samma flytt ---------------------------------------------

def test_shift_agreement():
    a = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    b = np.array([[0.05, 0.0], [1.0, 0.2], [0.0, 0.1]])
    # avvikelser 0.05, 0.2, 0.1 -> två av tre inom 0.1 km
    assert mn.shift_agreement(a, b, 0.1) == pytest.approx(2 / 3)


def test_implied_correlation():
    # kvoten median|e2−e1| / median|e| = sqrt(2(1−ρ)) för lika fördelade fel
    assert mn.implied_corr(np.sqrt(2 * (1 - 0.8))) == pytest.approx(0.8)
    rng = np.random.default_rng(3)
    rho, n = 0.75, 200_000
    z1 = rng.normal(size=(n, 2))
    z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.normal(size=(n, 2))
    ratio = np.median(np.hypot(*(z2 - z1).T)) / np.median(np.hypot(*z1.T))
    assert mn.implied_corr(ratio) == pytest.approx(rho, abs=0.01)


# --- status för ett påstående ----------------------------------------------------

def test_status():
    assert mn.status("0.7–2.6", "0.7–2.6") == "stämmer"
    assert mn.status("0.1–0.7", "0.1–0.6") == "avviker"
    assert mn.status("0.39", None) == "saknar källa"
    # bindestreck och minus i texten jämförs som tankstreck/minus
    assert mn.status("0.13-0.54", "0.13–0.54") == "stämmer"
    assert mn.status("-0.08", "−0.08") == "stämmer"


def test_is_tie():
    assert mn.is_tie(0.2850, 2)          # som det står i CSV: exakt halva
    assert not mn.is_tie(0.285024, 2)    # full precision: inte halva
    assert mn.is_tie(0.875, 0, 100)      # 87,5 %
    assert mn.is_tie(-0.015, 2)
    assert not mn.is_tie(0.4749, 2)
