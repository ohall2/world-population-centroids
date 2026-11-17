"""Tester av samregistreringen (wpc/coregister.py)."""
import numpy as np

from wpc import coregister as C


def blobs(n=600, seed=0):
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n].astype(float)
    a = np.zeros((n, n))
    for cy, cx, s, h in zip(rng.uniform(50, n - 50, 60), rng.uniform(50, n - 50, 60),
                            rng.uniform(2, 8, 60), rng.uniform(5, 63, 60)):
        a += h * np.exp(-((y - cy) ** 2 + (x - cx) ** 2) / (2 * s * s))
    return a


def test_recovers_subpixel_shift():
    ref = blobs()
    # img(x) = ref(x - d)  =>  ref(x) = img(x + d)
    img = C.shift_bilinear(ref, -0.4, 0.3)
    dx, dy, peak = C.tile_offset(ref, img)
    # img = ref sampled at (x-0.4, y+0.3) ⇒ ref(x) = img(x+0.4, y-0.3)
    assert abs(dx - 0.4) < 0.08 and abs(dy + 0.3) < 0.08 and peak > 0.9


def test_shift_moves_centre_of_mass_exactly():
    a = blobs(200, 1)
    b = C.shift_bilinear(a, 0.37, -0.21)
    y, x = np.mgrid[0:200, 0:200]
    com = lambda m: (np.sum(m * x) / m.sum(), np.sum(m * y) / m.sum())  # noqa: E731
    (xa, ya), (xb, yb) = com(a), com(b)
    assert np.isclose(b.sum(), a.sum(), rtol=1e-4)
    assert abs((xb - xa) + 0.37) < 1e-3 and abs((yb - ya) - 0.21) < 1e-3


def test_weighted_median():
    assert C.weighted_median(np.array([0.0, 1.0, 10.0]), np.array([1, 5, 1])) == 1.0
