"""Tester av A3-detektionsmodellen med konstgjorda data där svaret är känt.

Skrivna före analysskriptet (TDD). Facitläsning, endpoints-skattning och
Spearmans ρ återanvänds ur befintlig kod och prövas inte här. Dessa tester
gäller det som är nytt i A3:

  1. brusmåttet σ (DMSP inom perioden, kvadratsumma över skarven)
  2. logistisk regression (IRLS) med känt MLE
  3. AUC, Brier och tio-klasskalibrering
  4. minsta upptäckbara förflyttning (|t| där förutsagd sannolikhet är 0,8)
"""

import numpy as np
import pytest

from tools import detection_model as D


# 1. Brusmåttet σ -------------------------------------------------------------

def test_unit_sigma_within_dmsp():
    assert D.unit_sigma(0.84, 0.25, across_seam=False) == pytest.approx(0.84)


def test_unit_sigma_across_seam():
    # kvadratsumman: 3 och 4 ger 5
    assert D.unit_sigma(3.0, 4.0, across_seam=True) == pytest.approx(5.0)


# 2. Logistisk regression -----------------------------------------------------

def test_logistic_intercept_only_is_logit_of_rate():
    """Utan prediktor är MLE interceptet = logit(andelen ettor)."""
    y = np.array([0, 0, 0, 1])                      # andel 0,25
    beta = D.logistic_fit(np.empty((4, 0)), y)
    assert beta[0] == pytest.approx(np.log(0.25 / 0.75))


def test_logistic_symmetric_null_gives_zero_coeffs():
    """Prediktorn saknar information och andelen är 0,5 -> beta = (0, 0)."""
    x = np.array([-1.0, 1.0, -1.0, 1.0])
    y = np.array([0, 0, 1, 1])                       # vid varje x: en 0 och en 1
    beta = D.logistic_fit(x[:, None], y)
    assert beta[0] == pytest.approx(0.0, abs=1e-6)
    assert beta[1] == pytest.approx(0.0, abs=1e-6)


def test_logistic_recovers_known_direction():
    """Monotont samband men icke-separerbart (ändligt MLE): fler ettor vid större
    prediktor ger beta1 > 0, och förutsagd sannolikhet växer med prediktorn."""
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
    y = np.array([0, 0, 1, 0, 1, 1])                 # överlapp vid x = 0 och x = 1
    beta = D.logistic_fit(x[:, None], y)
    assert beta[1] > 0
    xs = np.linspace(-2, 3, 20)
    p = D.predict_prob(beta, xs[:, None])
    assert np.all(np.diff(p) > 0)


# 3. AUC, Brier, kalibrering --------------------------------------------------

def test_auc_known_example():
    scores = np.array([0.1, 0.4, 0.35, 0.8])
    y = np.array([0, 0, 1, 1])
    assert D.auc(scores, y) == pytest.approx(0.75)


def test_auc_perfect_and_reversed():
    s = np.array([0.1, 0.2, 0.9, 0.95]); y = np.array([0, 0, 1, 1])
    assert D.auc(s, y) == pytest.approx(1.0)
    assert D.auc(-s, y) == pytest.approx(0.0)


def test_brier_known():
    p = np.array([0.8, 0.2]); y = np.array([1, 0])
    assert D.brier(p, y) == pytest.approx(((0.8 - 1) ** 2 + (0.2 - 0) ** 2) / 2)


def test_calibration_bins_equal_count():
    p = np.array([0.1, 0.2, 0.8, 0.9])
    y = np.array([0, 0, 1, 1])
    out = D.calibration_bins(p, y, nbins=2)
    assert list(out["n"]) == [2, 2]
    assert out["pred_mean"].tolist() == pytest.approx([0.15, 0.85])
    assert out["obs_prop"].tolist() == pytest.approx([0.0, 1.0])


# 4. Minsta upptäckbara förflyttning -----------------------------------------

def test_min_detectable_move_known():
    """p = 0,8 vid logit(0,8) = beta0 + beta1·log(|t|/σ). Med beta0=0, beta1=1:
    log(|t|/σ) = logit(0,8) -> |t| = σ·0,8/0,2 = 4σ."""
    assert D.min_detectable_move(1.0, 0.0, 1.0, target=0.8) == pytest.approx(4.0)
    assert D.min_detectable_move(2.0, 0.0, 1.0, target=0.8) == pytest.approx(8.0)
