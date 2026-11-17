"""Datamängdens filformat (tools/build_release.py)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from build_release import GPKG_LAYER, sha256, write_csv, write_gpkg  # noqa: E402
from compare_series import read_results  # noqa: E402


def _frame():
    return pd.DataFrame({
        "unit_id": ["NA", "SWE", "SJM"], "COUNTRY": ["Namibia", "Sweden", "Svalbard"],
        "year": [2020, 2020, 2020], "sensor": ["VIIRS"] * 3,
        "bm_lon": [17.1, 15.0, np.nan], "bm_lat": [-22.5, 59.0, np.nan],
        "aa_lon": [17.2, 15.1, np.nan], "aa_lat": [-22.4, 59.1, np.nan],
        "bm_flag": [0, 0, -1], "max_lat": [65.0] * 3})


def test_csv_gz_roundtrip_and_reproducible(tmp_path):
    d = _frame()
    f = tmp_path / "x.csv.gz"
    write_csv(d, f)
    h = sha256(f)
    write_csv(d, f)
    assert sha256(f) == h  # bitidentisk (gzip utan tidsstämpel)
    back = read_results(f)
    assert back.unit_id.tolist() == ["NA", "SWE", "SJM"]  # "NA" är Namibia, inte saknat
    assert back.bm_lon.isna().tolist() == [False, False, True]


def test_gpkg_points_and_empty_geometry(tmp_path):
    gpd = pytest.importorskip("geopandas")
    d = _frame()
    g = tmp_path / "x.gpkg"
    write_gpkg(d, g)
    write_gpkg(d, g)  # skrivs över, inte dubbleras
    G = gpd.read_file(g, layer=GPKG_LAYER)
    assert len(G) == 3 and G.crs.to_epsg() == 4326
    assert G.geometry.iloc[:2].x.tolist() == pytest.approx([17.1, 15.0])
    assert G.geometry.iloc[2] is None or G.geometry.iloc[2].is_empty
    assert G.unit_id.tolist()[0] == "NA"
