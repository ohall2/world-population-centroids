"""Test av årsprodukter (CCNL/VIIRS) i hela pipelinen, med syntetiska data.

Ett globalt 1°-rutnät med ljus i två kända pixlar skrivs som en CCNL-fil med
samma metadata som tools/download_ee.py skriver. Pipelinen ska hitta filen
via metadata, läsa den i filens eget rutnät och ge den väntade centroiden.
"""
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from shapely.geometry import box

from wpc import run
from wpc.lights import Annual, available_years


def write_ccnl(raw: Path, year: int, data: np.ndarray, transform: Affine) -> None:
    raw.mkdir(parents=True, exist_ok=True)
    tif = raw / f"CCNL_{year}_b1.tif"
    with rasterio.open(tif, "w", driver="GTiff", dtype="float32", count=1,
                       width=data.shape[1], height=data.shape[0],
                       crs="EPSG:4326", transform=transform) as dst:
        dst.write(data.astype("float32"), 1)
    tif.with_suffix(".json").write_text(json.dumps({"year": year, "band": "b1"}))


def test_annual_end_to_end(tmp_path):
    t = Affine(1.0, 0, -180.0, 0, -1.0, 90.0)
    data = np.zeros((180, 360), "float32")
    # Pixlar med mitt i (10.5 E, 0.5 N) och (12.5 E, 0.5 N), lika vikt.
    data[89, 190] = 5.0
    data[89, 192] = 5.0
    data[50, 50] = -3.0                       # negativt värde ska bli 0
    raw = tmp_path / "raw" / "ccnl"
    write_ccnl(raw, 2005, data, t)
    write_ccnl(raw, 2006, data * 2, t)
    assert available_years("ccnl", raw, None) == [2005, 2006]

    src = Annual("ccnl", 2005, raw, None)
    assert (src.width, src.height) == (360, 180)
    w = src.read((-140.0, 35.0, -125.0, 45.0))
    assert w.data.min() == 0.0
    src.close()

    gpkg = tmp_path / "b.gpkg"
    gpd.GeoDataFrame({"GID_0": ["AAA"], "COUNTRY": ["Testland"]},
                     geometry=[box(5, -5, 20, 5)], crs="EPSG:4326").to_file(gpkg, layer="ADM_0")
    out = tmp_path / "res"
    run.main(["--source", "ccnl", "--raw-dir", str(raw), "--flares", "none",
              "--boundaries", str(gpkg), "--cache-dir", str(tmp_path / "cache"),
              "--levels", "0", "--out", str(out), "--workers", "1"])
    for year in (2005, 2006):
        df = pd.read_csv(out / "gadm41" / "ccnl" / "adm0" / f"{year}.csv")
        assert df.loc[0, "source"] == "ccnl"
        assert abs(df.loc[0, "aa_lon"] - 11.5) < 1e-6
        # Storcirkeln mellan punkterna buktar mot polen: ca 0.5 + 7.6e-5 grader.
        assert 0.5 < df.loc[0, "aa_lat"] < 0.5 + 1e-3
        assert df.loc[0, "aa_flag"] == 0
    meta = json.loads((out / "gadm41" / "ccnl" / "adm0" / "2005.json").read_text())
    assert meta["light_source"]["source"] == "ccnl"
