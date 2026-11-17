"""Test av skarvningen (tools/splice_series.py)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import splice_series as S  # noqa: E402


def frame(year, pts):
    return pd.DataFrame([{"unit_id": u, "year": year, "aa_lon": lo, "aa_lat": la,
                          "bm_lon": lo, "bm_lat": la} for u, (lo, la) in pts.items()])


def test_move_inverts_east_north():
    lon, lat = np.array([18.0, 179.9]), np.array([59.0, -10.0])
    de, dn = S.east_north_km(lon, lat, lon + 0.3, lat - 0.2)
    lon2, lat2 = S.move(lon, lat, de, dn)
    assert np.allclose(lat2, lat - 0.2) and np.allclose(lon2, (lon + 0.3 + 180) % 360 - 180, atol=1e-9)


def test_splice_is_continuous_and_keeps_dmsp(tmp_path):
    dmsp = {2012: frame(2012, {"A": (10.0, 50.0), "B": (20.0, 0.0)}),
            2013: frame(2013, {"A": (10.1, 50.0), "B": (20.0, 0.1)})}
    viirs = {2013: frame(2013, {"A": (10.5, 50.2), "B": (np.nan, np.nan)}),
             2014: frame(2014, {"A": (10.55, 50.2), "B": (21.0, 0.0)})}
    out, delta = S.splice_level(dmsp, viirs, 2013)
    assert set(out) == {2012, 2013, 2014}
    assert float(out[2013].set_index("unit_id").loc["A", "bm_lon"]) == 10.1
    # VIIRS(2013) + Δ = DMSP(2013) exakt
    d = delta.loc["A"]
    lo, la = S.move(10.5, 50.2, d.d_east_km_bm, d.d_north_km_bm)
    assert np.isclose(lo, 10.1) and np.isclose(la, 50.0)
    a = out[2014].set_index("unit_id").loc["A"]
    # VIIRS-hoppet 2013 -> 2014 (0.05° österut) ska bevaras, flyttat till DMSP 2013.
    de, dn = S.east_north_km(10.1, 50.0, a.bm_lon, a.bm_lat)
    de_v, dn_v = S.east_north_km(10.5, 50.2, 10.55, 50.2)
    assert abs(de - de_v) < 0.02 and abs(dn - dn_v) < 0.02
    b = out[2014].set_index("unit_id").loc["B"]
    assert b.splice_ok == 0 and b.bm_lon == 21.0            # kan inte skarvas: oförändrad
    assert a.splice_ok == 1


def test_read_results_keeps_NA_id(tmp_path):
    f = tmp_path / "x.csv"
    f.write_text("source,unit_id,COUNTRY,bm_lon,bm_lat\nv,NA,NA,1.5,\nv,SWE,Sweden,,2\n")
    d = S.read_results(f)
    assert list(d.unit_id) == ["NA", "SWE"] and d.COUNTRY[0] == "NA"
    assert d.bm_lon.dtype.kind == "f" and np.isnan(d.bm_lat[0]) and np.isnan(d.bm_lon[1])
