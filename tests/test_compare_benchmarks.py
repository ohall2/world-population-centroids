"""Test av den parvisa läges-jämförelsen i analys B (compare_benchmarks.pairwise)
med konstgjorda data där svaret är känt. Måtten och bootstrapen återanvänds och
prövas inte här."""

import numpy as np
import pytest

from tools.compare_benchmarks import pairwise


def test_pairwise_known():
    err_light = np.array([2.0, 2.0, 2.0, 2.0])
    err_source = np.array([1.0, 3.0, 1.0, 1.0])     # närmare, längre, närmare, närmare
    closer, diff = pairwise(err_light, err_source)
    assert closer == pytest.approx(0.75)             # 3 av 4 enheter
    assert diff == pytest.approx(1.0)                # median av [1, -1, 1, 1]


def test_pairwise_ties_not_counted_as_closer():
    err_light = np.array([1.0, 1.0, 1.0])
    err_source = np.array([1.0, 1.0, 0.5])           # lika räknas inte som närmare
    closer, diff = pairwise(err_light, err_source)
    assert closer == pytest.approx(1 / 3)
    assert diff == pytest.approx(0.0)                # median av [0, 0, 0.5]
