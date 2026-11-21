"""Test av Japan-läsaren (--dataset estat_jp) med känt svar. Prövar JIS-meshkodens
mittpunkt (1 km och 500 m) och att läsaren hittar befolkningskolumnen via den
japanska etiketten, härleder rutans mittpunkt ur meshkoden och låter hemlighållna
rutor (icke-numerisk befolkning) utgå."""

import numpy as np
import pytest

from tools.validation_centroids import _jis_mesh_centroid, read_estat_jp


def test_mesh_centroid_1km_known():
    # 1 km-mesh 53394547 (centrala Tokyo): mittpunkt (139,71875, 35,7041667)
    lon, lat = _jis_mesh_centroid("53394547")
    assert lon == pytest.approx(139.71875, abs=1e-6)
    assert lat == pytest.approx(35.70416667, abs=1e-6)


def test_mesh_centroid_500m_quadrants():
    # halvmesh: 1 = SV-kvadranten, 4 = NÖ-kvadranten av 1 km-rutan ovan
    lon1, lat1 = _jis_mesh_centroid("533945471")
    lon4, lat4 = _jis_mesh_centroid("533945474")
    assert (lon1, lat1) == pytest.approx((139.715625, 35.70208333), abs=1e-6)
    assert (lon4, lat4) == pytest.approx((139.721875, 35.70625), abs=1e-6)


def _write_mesh_csv(folder):
    folder.mkdir(parents=True, exist_ok=True)
    # Två rubrikrader (koder, japanska etiketter), sedan data; cp932 som e-Stat.
    lines = [
        "KEY_CODE,HTKSYORI,GASSAN,T001",
        ",秘匿処理,合算,人口（総数）",
        "53394547,,,100",
        "53394557,,,50",
        "53394567,1,,*",            # hemlighållen ruta: befolkning "*"
    ]
    (folder / "mesh.csv").write_text("\n".join(lines) + "\n", encoding="cp932")


def test_read_estat_jp(tmp_path):
    _write_mesh_csv(tmp_path / "2020")
    df, info = read_estat_jp(tmp_path, 2020)

    assert list(df.columns) == ["lon", "lat", "pop"]
    assert len(df) == 3                                   # en rad per mesh
    assert info["pop_column"] == "T001" and "人口" in info["pop_label"]
    assert info["mesh_len"] == 8
    assert set(info["secrecy_columns"]) == {"HTKSYORI", "GASSAN"}
    assert info["meshes_nonnumeric_pop"] == 1             # den hemlighållna rutan
    assert info["pop_total"] == pytest.approx(150.0)      # 100 + 50 (hemlig = NaN)
    # första rutans mittpunkt härleds ur meshkoden
    assert (df.lon.iloc[0], df.lat.iloc[0]) == pytest.approx((139.71875, 35.70416667), abs=1e-6)
    assert df["pop"].iloc[0] == pytest.approx(100.0)
    assert np.isnan(df["pop"].iloc[2])                    # hemlig ruta -> NaN (utgår senare)
