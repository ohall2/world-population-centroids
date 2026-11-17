#!/usr/bin/env python3
"""Fristående kontroll: gäller skarven DMSP->VIIRS bara centroider?

Oberoende av pipelinen (egen kod, egen rastrering med pixelmittpunkt i polygon).
För varje enhet och år beräknas ur momentsummor:
  T      total ljusvikt (ljus x pixelns yta)
  cent   Aboufadel–Austin-centroid (viktat medel av 3D-enhetsvektorer)
  SD     standardavstånd på sfären, R*sqrt(2(1-|m|)) (km)
  Neff   effektivt antal tända pixlar, exp(entropi av vikterna)
  nlit   antal tända pixlar
Årsförändringar: centroidflytt (km), |Δ log| för T, SD, Neff, nlit.

Användning: python tools/check_seam_measures.py <serie> <utfil.csv> <år:fil> ...
"""
import os, sys, time
import numpy as np, geopandas as gpd, rasterio
from rasterio.features import rasterize
from rasterio.windows import Window

R = 6371.0
# Indata: GADM 4.1 nivå 1, förenklad med tolerans 0,002° (shapely simplify,
# preserve_topology=False), och fackelmasken med 2 km buffert.
UNITS = os.environ.get("UNITS", "data/derived/check/adm1_simpl.gpkg")
THRESH = float(os.environ.get("THRESH", "0"))  # Li m.fl. rekommenderar DN > 7
FLARES = os.environ.get("FLARES", "data/processed/flares_viirs_r2km.gpkg")


def moments(path, units, flares, stripe=1024):
    src = rasterio.open(path)
    n = len(units) + 1
    acc = {k: np.zeros(n) for k in ("T", "x", "y", "z", "wlogw", "nlit")}
    shapes = list(zip(units.geometry, units["lab"]))
    fl = [(g, 1) for g in flares.geometry]
    t0 = time.time()
    for r0 in range(0, src.height, stripe):
        h = min(stripe, src.height - r0)
        win = Window(0, r0, src.width, h)
        tr = src.window_transform(win)
        lab = rasterize(shapes, out_shape=(h, src.width), transform=tr, fill=0, dtype="int32")
        if not lab.any():
            continue
        data = src.read(1, window=win).astype("float64")
        data[data < 0] = 0
        data[data <= THRESH] = 0
        if src.nodata is not None:
            data[src.read(1, window=win) == src.nodata] = 0
        fm = rasterize(fl, out_shape=(h, src.width), transform=tr, fill=0, dtype="uint8")
        data[fm == 1] = 0
        rows = np.arange(h)
        lat = tr.f + tr.e * (rows + 0.5)
        lon = tr.c + tr.a * (np.arange(src.width) + 0.5)
        area = (np.cos(np.radians(lat)))[:, None] * np.ones((1, src.width))
        w = data * area
        sel = (lab > 0) & (w > 0)
        if not sel.any():
            continue
        L = lab[sel]; W = w[sel]
        la = np.radians(np.broadcast_to(lat[:, None], lab.shape)[sel])
        lo = np.radians(np.broadcast_to(lon[None, :], lab.shape)[sel])
        cl = np.cos(la)
        for k, v in (("T", W), ("x", W * cl * np.cos(lo)), ("y", W * cl * np.sin(lo)),
                     ("z", W * np.sin(la)), ("wlogw", W * np.log(W)), ("nlit", np.ones_like(W))):
            acc[k] += np.bincount(L, weights=v, minlength=n)
        print(f"\r  {path.split('/')[-1]} rad {r0 + h}/{src.height}  {time.time() - t0:5.0f} s", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    return acc


def main():
    series, out = sys.argv[1], sys.argv[2]
    pairs = [a.split(":", 1) for a in sys.argv[3:]]
    units = gpd.read_file(UNITS, engine="pyogrio")
    units = units[~units.geometry.is_empty & units.geometry.notna()].reset_index(drop=True)
    units["lab"] = np.arange(1, len(units) + 1)
    flares = gpd.read_file(FLARES, engine="pyogrio")
    import pandas as pd
    rows = []
    for yr, path in pairs:
        a = moments(path, units, flares)
        T = a["T"][1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            mx, my, mz = a["x"][1:] / T, a["y"][1:] / T, a["z"][1:] / T
            mn = np.sqrt(mx**2 + my**2 + mz**2)
            H = np.log(T) - a["wlogw"][1:] / T
            rows.append(pd.DataFrame({
                "series": series, "year": int(yr), "GID_1": units.GID_1, "T": T,
                "lon": np.degrees(np.arctan2(my, mx)), "lat": np.degrees(np.arcsin(mz / mn)),
                "SD_km": R * np.sqrt(np.clip(2 * (1 - mn), 0, None)),
                "Neff": np.exp(H), "nlit": a["nlit"][1:]}))
    pd.concat(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
