"""Test av download_ee.py utan nätverk: en falsk ee-modul levererar pixlar
ur en känd global matris, så att hela kedjan (listning, rutnät, rutor,
sammansättning, metadata) kan kontrolleras exakt."""
import json, sys, types
from pathlib import Path
import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

W, H = 1000, 413                      # udda storlek: rutorna går inte jämnt ut
SX, SY = 0.0083333333333333, -0.0083333333333333
TX, TY = -180.00416666666666, 75.00416666666666
TRUTH = (np.arange(W * H, dtype=np.int64).reshape(H, W) % 64).astype(np.uint8)
TRUTH[0:3, 0:5] = 0                   # "maskade" pixlar blir 0 via unmask

def make_fake_ee():
    ee = types.ModuleType("ee"); ee.__version__ = "fake"
    calls = {"init": None, "n": 0}
    class Img:
        def __init__(self, aid): self.aid = aid
        def select(self, b): self.band = b; return self
        def unmask(self, v): return self
        def toUint8(self): return self
        def getInfo(self):
            return {"bands": [{"id": "stable_lights", "crs": "EPSG:4326",
                               "crs_transform": [SX, 0, TX, 0, SY, TY],
                               "dimensions": [W, H],
                               "data_type": {"type": "PixelType", "precision": "int", "min": 0, "max": 63}}]}
    class Coll:
        def __init__(self, c): pass
        def aggregate_array(self, p):
            class R:
                def __init__(s, v): s.v = v
                def getInfo(s): return s.v
            if p == "system:index":
                return R(["F101992", "F101994", "F121994", "F182010"])
            return R([694224000000, 757382400000, 757382400000, 1262304000000])
    ee.Image = Img; ee.ImageCollection = Coll
    def Initialize(project=None, url=None): calls["init"] = (project, url)
    ee.Initialize = Initialize
    data = types.ModuleType("ee.data")
    def computePixels(req):
        calls["n"] += 1
        g = req["grid"]; a = g["affineTransform"]
        col0 = round((a["translateX"] - TX) / SX); row0 = round((a["translateY"] - TY) / SY)
        w, h = g["dimensions"]["width"], g["dimensions"]["height"]
        assert req["fileFormat"] == "NUMPY_NDARRAY" and g["crsCode"] == "EPSG:4326"
        assert a["scaleX"] == SX and a["scaleY"] == SY
        out = np.zeros((h, w), dtype=[("stable_lights", "u1")])
        out["stable_lights"] = TRUTH[row0:row0 + h, col0:col0 + w]
        return out
    data.computePixels = computePixels
    data.setDeadline = lambda ms: calls.__setitem__("deadline", ms)
    ee.data = data
    return ee, calls

def test_end_to_end(tmp_path, monkeypatch):
    fake, calls = make_fake_ee()
    monkeypatch.setitem(sys.modules, "ee", fake)
    import download_ee
    download_ee.main(["--project", "test-proj", "--out", str(tmp_path),
                      "--years", "1994", "2010", "--tile", "128", "--workers", "4"])
    assert calls["init"] == ("test-proj", download_ee.HIGHVOLUME_URL)
    assert calls["deadline"] == download_ee.REQUEST_TIMEOUT_MS
    files = sorted(p.name for p in tmp_path.glob("*.tif"))
    assert files == ["DMSP_F101994_stable_lights.tif",
                     "DMSP_F121994_stable_lights.tif",
                     "DMSP_F182010_stable_lights.tif"]
    n_tiles = len(download_ee.make_tiles(W, H, 128))
    assert calls["n"] == 3 * n_tiles
    with rasterio.open(tmp_path / "DMSP_F182010_stable_lights.tif") as r:
        assert (r.width, r.height, r.dtypes[0]) == (W, H, "uint8")
        assert r.crs.to_string() == "EPSG:4326"
        assert np.allclose(list(r.transform)[:6], [SX, 0, TX, 0, SY, TY])
        assert np.array_equal(r.read(1), TRUTH)
    meta = json.loads((tmp_path / "DMSP_F182010_stable_lights.json").read_text())
    assert meta["year"] == 2010 and meta["sha256"] and meta["n_tiles"] == n_tiles
    assert not list(tmp_path.glob("*.part"))
    # Andra körningen hoppar över färdiga filer
    calls["n"] = 0
    download_ee.main(["--project", "test-proj", "--out", str(tmp_path),
                      "--years", "2010", "--tile", "128"])
    assert calls["n"] == 0

def test_tiles_cover_grid_exactly():
    import download_ee
    tiles = download_ee.make_tiles(W, H, 128)
    cover = np.zeros((H, W), int)
    for c, r, w, h in tiles: cover[r:r+h, c:c+w] += 1
    assert (cover == 1).all()

def test_dtype_and_tile_size():
    import download_ee as d
    assert d.numpy_dtype({"precision": "int", "min": 0, "max": 63}) == "uint8"
    assert d.numpy_dtype({"precision": "int", "min": 0, "max": 1000}) == "uint16"
    assert d.numpy_dtype({"precision": "int", "min": -5, "max": 1000}) == "int16"
    assert d.numpy_dtype({"precision": "float"}) == "float32"
    assert d.numpy_dtype({"precision": "double"}) == "float32"
    assert d.tile_size("uint8") == 4096 and d.tile_size("float32") == 2048

def test_parse_years():
    import download_ee as d
    assert d.parse_years(["all"]) is None
    assert d.parse_years(["1992-1994", "2010"]) == [1992, 1993, 1994, 2010]


def test_year_from_id():
    import download_ee as d
    assert d.year_from_id("Harmonized_DN_NTL_2013_simVIIRS") == 2013
    assert d.year_from_id("SVNL_1992") == 1992
    import pytest
    with pytest.raises(ValueError):
        d.year_from_id("noyear")
