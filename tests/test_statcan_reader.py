"""Test av StatCan GAF-läsaren (--dataset statcan) med en syntetisk fil där
svaret är känt. Läsaren lägger ihop dissemination block-rader till dissemination
areas: en punkt per DA (DA-representativ punkt) med DA-befolkning = summan av
blockens. Kolumnnamnen matchas med regex (tvåspråkiga GAF-rubriker)."""

import pandas as pd
import pytest

from tools.validation_centroids import read_statcan


def _write_gaf(folder, rows):
    folder.mkdir(parents=True, exist_ok=True)
    cols = ["DBUID_IDIDU", "DAUID_ADIDU", "DBPOP2021_IDPOP2021", "DARPLAT_ADLAT", "DARPLONG_ADLONG"]
    pd.DataFrame(rows, columns=cols).to_csv(folder / "gaf.csv", index=False)


def test_read_statcan_aggregates_blocks_to_da(tmp_path):
    # DA "A": två block (10 + 5) vid punkten (45,0, -75,0); DA "B": ett block (20);
    # DA "C": ett block med 0 invånare (ska bli kvar med pop 0 i läsaren).
    _write_gaf(tmp_path / "2021", [
        ["35000001", "35000A", 10, 45.0, -75.0],
        ["35000002", "35000A", 5, 45.0, -75.0],
        ["35000003", "35000B", 20, 46.0, -76.0],
        ["35000004", "35000C", 0, 47.0, -77.0],
    ])
    df, info = read_statcan(tmp_path, 2021)

    assert list(df.columns) == ["lon", "lat", "pop"]
    assert len(df) == 3                                  # en rad per DA
    assert info["blocks"] == 4 and info["dissemination_areas"] == 3
    a, b, c = df.iloc[0], df.iloc[1], df.iloc[2]         # ordning som i filen
    assert (a.lon, a.lat, a["pop"]) == pytest.approx((-75.0, 45.0, 15.0))
    assert (b.lon, b.lat, b["pop"]) == pytest.approx((-76.0, 46.0, 20.0))
    assert c["pop"] == pytest.approx(0.0)
    assert info["pop_total"] == pytest.approx(35.0)


def test_read_statcan_missing_column_errors(tmp_path):
    folder = tmp_path / "2016"
    folder.mkdir(parents=True)
    pd.DataFrame({"DAUID": ["1"], "DARPLAT": [45.0], "DARPLONG": [-75.0]}).to_csv(
        folder / "gaf.csv", index=False)      # saknar befolkningskolumn
    with pytest.raises(ValueError, match="ingen kolumn matchar"):
        read_statcan(tmp_path, 2016)
