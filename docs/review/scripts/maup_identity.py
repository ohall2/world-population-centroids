#!/usr/bin/env python3
"""Granskning: hur beror uppdelningen Bb/Bw i check_identity.py på regionindelningen?

Regionerna (GADM nivå 1) slås ihop till k grupper med k-means på
regionernas tyngdpunkter (k = 1 ... antal regioner), och identiteten
räknas om. Bb + Bw är oberoende av indelningen; fördelningen mellan dem är det inte.
"""
import sys, json, tempfile, os
import numpy as np, geopandas as gpd, pandas as pd
sys.path.insert(0, "/home/claude/wpc/tools")
import check_identity as ci

country, facit, ffile, viirs, units, flares, grid = sys.argv[1:8]
grid = float(grid)
max_lat = float(sys.argv[8]) if len(sys.argv) > 8 else None
pts = (ci.read_inegi(ffile and __import__("pathlib").Path(ffile)) if facit == "inegi"
       else ci.read_geostat_csv(ffile) if facit == "geostat"
       else pd.read_csv(ffile, usecols=["lon", "lat", "pop"]))
a0 = gpd.read_file(units, layer="adm0"); a0 = a0[a0.GID_0 == country]
a1 = gpd.read_file(units, layer="adm1"); a1 = a1[a1.GID_0 == country].reset_index(drop=True)
xy = np.c_[a1.geometry.representative_point().x, a1.geometry.representative_point().y]
from scipy.cluster.vq import kmeans2
res = []
for k in sorted({1, 2, 3, 5, 8, len(a1) // 2, len(a1)}):
    if k < len(a1):
        np.random.seed(0)
        _, lab = kmeans2(xy, k, minit="++", seed=0)
        g = a1.assign(grp=lab).dissolve("grp").reset_index()
        g["GID_0"] = country; g["NAME_1"] = g.grp.astype(str); g["GID_1"] = g.NAME_1
    else:
        g = a1
    tmp = tempfile.mktemp(suffix=".gpkg")
    a0.to_file(tmp, layer="adm0", driver="GPKG")
    g[["GID_0", "GID_1", "NAME_1", "geometry"]].to_file(tmp, layer="adm1", driver="GPKG")
    o = ci.run(country, pts, viirs, tmp, flares, grid, max_lat)
    os.remove(tmp)
    res.append({"k": int(g.shape[0]), "D": o["D_total"]["length_km"],
                "Bb": o["Bb_between_regions"]["length_km"], "Bw": o["Bw_within_regions"]["length_km"],
                "A": o["A_light_unpopulated"]["length_km"]})
    print(res[-1], flush=True)
json.dump(res, open(f"/home/claude/review/out/maup_{country}_{max_lat}.json", "w"), indent=1)
