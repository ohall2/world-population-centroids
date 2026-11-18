#!/usr/bin/env python3
"""
Kontroll inför artikeln: varför ligger ljusets centroid fel på landsnivå?

Ljusets centroid (viktat medel av 3D-enhetsvektorer, Aboufadel–Austin) och
befolkningens skiljer sig med vektorn D = m_L − m_P. Med ett gemensamt
rutnät (celler c) och ljus per invånare r_c = L_c / P_c i befolkade celler
delas D exakt i fyra delar (Price-identiteten, Price 1972, ekv. A7):

  D = A + s·(Bb + Bw + C)

  A   ljus i celler utan befolkning: (1 − s)(m_Ltom − m_P),
      s = andelen av ljuset som ligger i befolkade celler
  Bb  skillnader i ljus per invånare MELLAN regioner (GADM nivå 1):
      Cov_P(u, r̄_g) / E_P(r)
  Bw  skillnader i ljus per invånare INOM regioner:
      E_P[Cov_P(u, r | g)] / E_P(r)
  C   läget inom cellerna: ljusets och befolkningens tyngdpunkt i samma
      cell skiljer sig, Σ_c (L_c / L_bef)(u_c^L − u_c^P)

u är enhetsvektorer; m är viktade medelvektorer. Delarna redovisas i km
öster och norr i tangentplanet vid befolkningens centroid (R·T/|m_P|).
Summan jämförs med det exakta storcirkelavståndet mellan centroiderna.

Indata: VIIRS VNL v2 (average_masked), fackelmask, GADM 4.1 nivå 0 och 1
för landet, och facit som punkter (lon, lat, pop).

Användning (se docs/HARMONIZATION.md §14):
  python tools/check_identity.py --country MEX --year 2020 --facit inegi ...
"""
from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

R = 6371.0


def unit_vec(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def to_lonlat(m):
    return np.degrees(np.arctan2(m[1], m[0])), np.degrees(np.arcsin(m[2] / np.linalg.norm(m)))


def gc_km(a, b):
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    return R * np.arccos(np.clip(a @ b, -1, 1))


def read_inegi(path: Path):
    """Samma läsare som facit i tools/validation_centroids.py."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from validation_centroids import read_inegi as _read
    year = int(re.search(r"(20\d\d)", path.parent.name or path.name)[1]) if re.search(r"(20\d\d)", path.parent.name) else 2020
    d, _ = _read(path.parent.parent, year)
    return d.dropna()


def read_geostat_csv(path: Path):
    from pyproj import Transformer
    d = pd.read_csv(path, dtype={"GRD_ID": str, "CNTR_ID": str})
    d = d[(d["T"] > 0) & ~d.GRD_ID.str.endswith("_unallocated")]
    m = d.GRD_ID.str.extract(r"N(\d+)E(\d+)$")
    n, e = m[0].astype(float) + 500.0, m[1].astype(float) + 500.0
    lon, lat = Transformer.from_crs(3035, 4326, always_xy=True).transform(e.values, n.values)
    return pd.DataFrame({"lon": lon, "lat": lat, "pop": d["T"].astype(float).values})


def run(country, pts, viirs, units, flares, grid, max_lat=None):
    import geopandas as gpd
    import rasterio
    from rasterio.features import rasterize
    from rasterio.windows import from_bounds

    a0 = gpd.read_file(units, layer="adm0", engine="pyogrio")
    a0 = a0[a0.GID_0 == country]
    a1 = gpd.read_file(units, layer="adm1", engine="pyogrio")
    a1 = a1[a1.GID_0 == country].reset_index(drop=True)
    geom = a0.geometry.iloc[0]
    minx, miny, maxx, maxy = geom.bounds
    if country in ("NOR",):  # Svalbard och Jan Mayen ingår inte i GEOSTAT
        maxx, maxy = min(maxx, 32.0), min(maxy, 72.0)
    src = rasterio.open(viirs)
    win = from_bounds(minx, miny, maxx, maxy, src.transform).round_offsets().round_lengths()
    tr = src.window_transform(win)
    shp = (int(win.height), int(win.width))
    fl = gpd.read_file(flares, engine="pyogrio", bbox=(minx, miny, maxx, maxy))
    flshapes = [(g, 1) for g in fl.geometry]
    # Läs i band om 1024 rader för att hålla nere minnet; spara bara tända pixlar.
    Lw, Llon, Llat = [], [], []
    lon_all = tr.c + tr.a * (np.arange(shp[1]) + 0.5)
    for r0 in range(0, shp[0], 1024):
        h = min(1024, shp[0] - r0)
        sw = rasterio.windows.Window(win.col_off, win.row_off + r0, shp[1], h)
        st = src.window_transform(sw)
        light = src.read(1, window=sw).astype("float64")
        light[~np.isfinite(light) | (light < 0)] = 0
        cm = rasterize([(geom, 1)], out_shape=(h, shp[1]), transform=st, fill=0, dtype="uint8").astype(bool)
        light[~cm] = 0
        if flshapes:
            light[rasterize(flshapes, out_shape=(h, shp[1]), transform=st, fill=0, dtype="uint8") == 1] = 0
        lat_r = st.f + st.e * (np.arange(h) + 0.5)
        if max_lat is not None:
            light[lat_r > max_lat] = 0
        w = light * np.cos(np.radians(lat_r))[:, None]
        rr, cc = np.nonzero(w > 0)
        Lw.append(w[rr, cc]); Llon.append(lon_all[cc]); Llat.append(lat_r[rr])
    Lw, Llon, Llat = np.concatenate(Lw), np.concatenate(Llon), np.concatenate(Llat)
    cmask_geom = geom

    # befolkningspunkter i landet och deras region (punkt i polygon)
    pg = gpd.GeoDataFrame(pts.reset_index(drop=True), geometry=gpd.points_from_xy(pts.lon, pts.lat), crs=4326)
    j = gpd.sjoin(pg, a1[["geometry"]].assign(reg=a1.index + 1), how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    inside = j.reg.notna()               # punkt i någon av landets regioner
    P = pts.reset_index(drop=True)[inside.values].copy()
    P["reg"] = j.reg[inside.values].fillna(0).astype(int).values

    # gemensamt rutnät
    x0, y0 = tr.c, tr.f
    def cell(lo, la):
        return np.floor((lo - x0) / grid).astype(np.int64) * 100000 + np.floor((y0 - la) / grid).astype(np.int64)
    cl, cp = cell(Llon, Llat), cell(P.lon.values, P.lat.values)
    uL, uP = unit_vec(Llon, Llat), unit_vec(P.lon.values, P.lat.values)
    pw = P["pop"].values

    # summor per cell
    cells = np.union1d(cl, cp)
    iL, iP = np.searchsorted(cells, cl), np.searchsorted(cells, cp)
    n = len(cells)
    Lc = np.bincount(iL, Lw, n)
    Pc = np.bincount(iP, pw, n)
    VL = np.stack([np.bincount(iL, Lw * uL[:, k], n) for k in range(3)], 1)
    VP = np.stack([np.bincount(iP, pw * uP[:, k], n) for k in range(3)], 1)
    # region per cell: befolkningens största region i cellen
    regc = np.zeros(n, int)
    dfr = pd.DataFrame({"i": iP, "reg": P.reg.values, "w": pw}).groupby(["i", "reg"]).w.sum().reset_index()
    dfr = dfr.sort_values("w").drop_duplicates("i", keep="last")
    regc[dfr.i.values] = dfr.reg.values

    popd = Pc > 0
    Ltot, Ptot = Lc.sum(), Pc.sum()
    Lpop = Lc[popd].sum()
    s = Lpop / Ltot
    mL = VL.sum(0) / Ltot
    mP = VP.sum(0) / Ptot
    mLpop = VL[popd].sum(0) / Lpop
    mLtom = VL[~popd].sum(0) / (Ltot - Lpop) if Ltot > Lpop else mP
    A = (1 - s) * (mLtom - mP)
    # celltyngdpunkter
    uPc = VP[popd] / Pc[popd, None]
    uLc = np.where(Lc[popd, None] > 0, VL[popd] / np.maximum(Lc[popd, None], 1e-300), uPc)
    rc = Lc[popd] / Pc[popd]
    q = Pc[popd] / Ptot                     # befolkningsandelar
    Er = (q * rc).sum()
    # regionmedel av r (befolkningsviktat) = L_g / P_g
    g = regc[popd]
    Lg = pd.Series(Lc[popd]).groupby(g).transform("sum").values
    Pg = pd.Series(Pc[popd]).groupby(g).transform("sum").values
    rg = Lg / Pg
    Eu = (q[:, None] * uPc).sum(0)
    cov = lambda a: (q[:, None] * (uPc - Eu) * (a - Er)[:, None]).sum(0)
    Bb = cov(rg) / Er
    Bw = cov(rc - rg + Er) / Er            # = cov(rc) - cov(rg)
    C = ((Lc[popd] / Lpop)[:, None] * (uLc - uPc)).sum(0)
    B = s * (Bb + Bw + C)
    D = mL - mP
    assert np.allclose(A + B, D, atol=1e-12), (A + B, D)

    # tangentplan vid befolkningens centroid
    lo0, la0 = to_lonlat(mP)
    e_hat = np.array([-np.sin(np.radians(lo0)), np.cos(np.radians(lo0)), 0.0])
    n_hat = np.array([-np.sin(np.radians(la0)) * np.cos(np.radians(lo0)),
                      -np.sin(np.radians(la0)) * np.sin(np.radians(lo0)), np.cos(np.radians(la0))])
    scale = R / np.linalg.norm(mP)
    def en(T):
        return float(T @ e_hat * scale), float(T @ n_hat * scale)
    out = {"country": country, "grid_deg": grid, "max_lat": max_lat,
           "pop_points": int(len(P)), "pop": float(Ptot),
           "light_share_in_populated_cells": float(s),
           "error_gc_km": float(gc_km(mL, mP)),
           "light_centroid": to_lonlat(mL), "pop_centroid": (lo0, la0)}
    for k, T in (("D_total", D), ("A_light_unpopulated", A), ("Bb_between_regions", s * Bb),
                 ("Bw_within_regions", s * Bw), ("C_within_cells", s * C)):
        e, nn = en(T)
        out[k] = {"east_km": round(e, 2), "north_km": round(nn, 2), "length_km": round(float(np.hypot(e, nn)), 2)}
    # ljus per invånare per region, för tolkning
    rr = (pd.DataFrame({"reg": g, "L": Lc[popd], "P": Pc[popd]}).groupby("reg").sum())
    rr["r_rel"] = (rr.L / rr.P) / (Lpop / Ptot)
    rr["name"] = [a1.NAME_1.iloc[i - 1] if i > 0 else "(ingen)" for i in rr.index]
    out["light_per_capita_by_region_rel"] = rr.sort_values("r_rel")[["name", "r_rel", "P"]].round(3).to_dict("records")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--country", required=True)
    p.add_argument("--facit", choices=["inegi", "geostat", "points"], required=True)
    p.add_argument("--facit-file", type=Path, required=True)
    p.add_argument("--viirs", type=Path, required=True)
    p.add_argument("--units", type=Path, required=True, help="GeoPackage med lagren adm0 och adm1")
    p.add_argument("--flares", type=Path, required=True)
    p.add_argument("--grid", type=float, default=1 / 30, help="Gemensamt rutnät i grader (standard 2')")
    p.add_argument("--max-lat", type=float, default=None)
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args(argv)
    pts = (read_inegi(a.facit_file) if a.facit == "inegi" else read_geostat_csv(a.facit_file)
           if a.facit == "geostat" else pd.read_csv(a.facit_file, usecols=["lon", "lat", "pop"]))
    res = run(a.country, pts, a.viirs, a.units, a.flares, a.grid, a.max_lat)
    txt = json.dumps(res, indent=1, ensure_ascii=False)
    print(txt)
    if a.out:
        a.out.write_text(txt)


if __name__ == "__main__":
    main()
