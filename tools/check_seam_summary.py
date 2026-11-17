#!/usr/bin/env python3
"""Sammanfatta seam_measures-utfilerna: skarv 2013->2014 mot vanliga årshopp."""
import sys
import numpy as np, pandas as pd

PAIRS = [(2010, 2011), (2011, 2012), (2012, 2013), (2013, 2014), (2014, 2015), (2015, 2016), (2016, 2017)]
SEAM = (2013, 2014)


def hav(lon1, lat1, lon2, lat2):
    la1, lo1, la2, lo2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(a))


def summarize(df, min_nlit=10):
    w = df.pivot(index="GID_1", columns="year")
    # enheter med minst min_nlit tända pixlar alla år
    ok = (w["nlit"] >= min_nlit).all(axis=1)
    w = w[ok]
    out = {}
    for a, b in PAIRS:
        d = {
            "centroid_km": hav(w["lon"][a], w["lat"][a], w["lon"][b], w["lat"][b]),
            "SD_log": np.abs(np.log(w["SD_km"][b] / w["SD_km"][a])),
            "Neff_log": np.abs(np.log(w["Neff"][b] / w["Neff"][a])),
            "nlit_log": np.abs(np.log(w["nlit"][b] / w["nlit"][a])),
            "T_log": np.abs(np.log(w["T"][b] / w["T"][a])),
            # signerad förändring: systematisk riktning vid skarven?
            "SD_logs": np.log(w["SD_km"][b] / w["SD_km"][a]),
            "Neff_logs": np.log(w["Neff"][b] / w["Neff"][a]),
            "T_logs": np.log(w["T"][b] / w["T"][a]),
        }
        out[(a, b)] = {k: float(np.nanmedian(v)) for k, v in d.items()}
    t = pd.DataFrame(out).T
    t.index = [f"{a}-{b}" for a, b in t.index]
    normal = t.drop(index=f"{SEAM[0]}-{SEAM[1]}")
    ratio = t.loc[f"{SEAM[0]}-{SEAM[1]}"] / normal.median()
    return t, ratio, int(ok.sum())


if __name__ == "__main__":
    for f in sys.argv[1:]:
        df = pd.read_csv(f)
        t, ratio, n = summarize(df)
        print(f"\n=== {f}  (enheter: {n})")
        print(t.round(4).to_string())
        print("skarv / median av övriga hopp (absoluta mått):")
        print(ratio[["centroid_km", "SD_log", "Neff_log", "nlit_log", "T_log"]].round(2).to_string())
