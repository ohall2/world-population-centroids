"""Fig. 6: årshopp i centroiden per enhet (GADM nivå 1), 2010–2017, fem serier.
Gemensamt urval: enheter med minst 10 tända pixlar alla år i alla fyra jämförelseserier
(docs/check_seam/*.csv.gz, Aboufadel–Austin) och med centroid i v2.0.0 (aa_lon/aa_lat) alla år.
Obs: jämförelseserierna är beräknade med tools/check_seam_measures.py (förenklade polygoner),
v2.0.0 med pipelinen; skillnaden i centroid är liten (HARMONIZATION §12: 4,87 mot 4,75 km för Li).
v2.0.0 läses från utgåvans adm1-fil; sökväg som argument (standard data/release/v2.0.0)."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R = Path(__file__).resolve().parents[2]; S = R / "docs/check_seam"; O = R / "docs/figures"
REL = Path(sys.argv[1]) if len(sys.argv) > 1 else R / "data/release/v2.0.0"
Y = list(range(2010, 2018))
def hav(lo1, la1, lo2, la2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*6371*np.arcsin(np.sqrt(a))
series = {"li": "Li et al. 2020", "li_dn7": "Li et al. 2020, DN > 7", "chiov": "Chiovelli et al. 2026", "pcnl": "PCNL (Li, S. et al. 2023)"}
W = {}
for k in series:
    w = pd.read_csv(S / f"{k}.csv.gz").pivot(index="GID_1", columns="year")
    W[k] = w[(w["nlit"] >= 10).all(axis=1)]
r = pd.read_csv(REL / "wpc_v2_gadm41_adm1.csv.gz", keep_default_na=False, usecols=["unit_id", "year", "aa_lon", "aa_lat"])
r = r[r.year.isin(Y)]
for c in ("aa_lon", "aa_lat"): r[c] = pd.to_numeric(r[c], errors="coerce")
v = r.pivot(index="unit_id", columns="year").dropna()
v = v.rename(columns={"aa_lon": "lon", "aa_lat": "lat"}, level=0)
W["v2"] = v; series["v2"] = "This series (v2.0.0)"
ids = set.intersection(*[set(w.index) for w in W.values()])
ids = sorted(ids); print("common units:", len(ids))
rows = []
for k, w in W.items():
    w = w.loc[ids]
    for a, b in zip(Y[:-1], Y[1:]):
        rows.append({"series": k, "pair": f"{a}–{str(b)[2:]}", "a": a,
                     "median_km": float(np.median(hav(w["lon"][a], w["lat"][a], w["lon"][b], w["lat"][b])))})
t = pd.DataFrame(rows); t.to_csv(O / "data/fig6_seam_jumps.csv", index=False)
piv = t.pivot(index="series", columns="a", values="median_km")
ratio = piv[2013] / piv[[2010, 2011, 2012]].median(axis=1)
print(piv.round(2).to_string()); print("seam / median DMSP-era jump:", ratio.round(2).to_dict())
col = {"li": "#e87ba4", "li_dn7": "#eb6834", "chiov": "#6250d6", "pcnl": "#1baf7a", "v2": "#2a78d6"}
mk = {"li": "v", "li_dn7": "^", "chiov": "s", "pcnl": "D", "v2": "o"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.axvspan(2.5, 3.5, color="#f0efec", zorder=0)
ax.text(3, 0.12, "sensor change\nDMSP → VIIRS", ha="center", va="bottom", fontsize=7.5, color="#52514e")
x = np.arange(len(Y) - 1)
for k in ["li", "li_dn7", "chiov", "pcnl", "v2"]:
    y = piv.loc[k, Y[:-1]].values
    ax.plot(x, y, color=col[k], marker=mk[k], ms=5, lw=2.4 if k == "v2" else 1.6,
            label=f"{series[k]}  (seam ×{ratio[k]:.1f})", zorder=4 if k == "v2" else 3,
            markeredgecolor="white", markeredgewidth=0.8)
ax.set_xticks(x); ax.set_xticklabels([f"{a}–{str(b)[2:]}" for a, b in zip(Y[:-1], Y[1:])])
ax.set_ylabel("Median year-to-year shift of the\nlight-weighted centroid (km)"); ax.set_ylim(0, 5.3)
ax.set_xlabel("Year pair (GADM level 1 units, n = %d)" % len(ids))
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=7.3, loc="upper right", bbox_to_anchor=(1.0, 0.82))
fig.text(0.01, 0.005, "Seam ratio: 2013–14 shift divided by the median of the 2010–11, 2011–12 and 2012–13 shifts.",
         fontsize=6.8, color="#52514e")
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("pdf", "png"): fig.savefig(O / f"fig6_seam.{ext}", dpi=300)
