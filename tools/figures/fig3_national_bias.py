"""Fig. 3: varför landscentroiderna ligger fel.
(a) Uppdelning av avvikelsen D = A + Bb + Bw + C (längder i km), rutnät 10′, utan maskning.
(b) Mexiko: delarnas längd beroende på rutnätet (1′–30′).
(c) Ljus per invånare per region i förhållande till landets medel (10′).
Läser docs/check_identity/id_*.json (tools/check_identity.py)."""
import json, numpy as np, matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R = Path(__file__).resolve().parents[2]; I = R / "docs/check_identity"; O = R / "docs/figures"
J = lambda k: json.load(open(I / f"id_{k}.json"))
C = [("Mexico", "MEX_0.1666667"), ("Brazil", "BRA_0.1666667"), ("Sweden", "SWE_none"),
     ("Norway", "NOR_none"), ("Finland", "FIN_none")]
parts = [("D_total", "Total error D", "#52514e"), ("Bb_between_regions", "Between regions (Bb)", "#2a78d6"),
         ("Bw_within_regions", "Within regions (Bw)", "#eb6834"), ("A_light_unpopulated", "Light in unpopulated cells (A)", "#1baf7a")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], width_ratios=[1.25, 1], hspace=0.65, wspace=0.35)
ax = fig.add_subplot(gs[0, 0])
h = 0.19
for i, (lab, k) in enumerate(C):
    d = J(k); y = len(C) - 1 - i
    for j, (key, pl, col) in enumerate(parts):
        v = d[key]["length_km"]
        ax.barh(y + (1.5 - j)*h, v, height=h*0.9, color=col, label=pl if i == 0 else None)
        if key in ("D_total",): ax.text(v + 2, y + 1.5*h, f"{v:.0f}", va="center", fontsize=7)
ax.set_yticks(range(len(C))); ax.set_yticklabels([c[0] for c in C][::-1])
ax.set_xlabel("Length of component (km)"); ax.set_xlim(0, 170)
ax.set_title("a  National error by component (10′)", loc="left", fontsize=9, fontweight="bold")
ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2)
ax2 = fig.add_subplot(gs[0, 1])
gr = ["0.0166667", "0.0333333", "0.0833333", "0.1666667", "0.5"]; gl = [1, 2, 5, 10, 30]
for key, pl, col in parts[1:]:
    ax2.plot(gl, [J(f"MEX_{g}")[key]["length_km"] for g in gr], marker="o", ms=4, lw=2, color=col, label=pl)
ax2.axhline(J("MEX_0.1666667")["D_total"]["length_km"], color="#52514e", lw=1, ls=(0, (4, 3)))
ax2.text(30, 134, "D = 131 km", ha="right", fontsize=7, color="#52514e")
ax2.set_xscale("log"); ax2.set_xticks(gl); ax2.set_xticklabels([f"{g}′" for g in gl])
ax2.set_xlabel("Grid cell size"); ax2.set_ylabel("km"); ax2.set_ylim(0, 145)
ax2.set_title("b  Mexico: split by grid size", loc="left", fontsize=9, fontweight="bold")
ax3 = fig.add_subplot(gs[1, :])
lab_pick = {"Mexico": {"Distrito Federal": 1, "Sonora": 1},
            "Brazil": {"Amazonas": 1, "Mato Grosso": 1},
            "Sweden": {"Stockholm": 1, "Norrbotten": 1}, "Norway": {"Oslo": 1, "Troms": 1, "Finnmark": 1},
            "Finland": {"Southern Finland": 1, "Lapland": 1}}
for i, (lab, k) in enumerate(C):
    r = J(k)["light_per_capita_by_region_rel"]; y = len(C) - 1 - i
    v = np.array([x["r_rel"] for x in r]); p = np.array([x["P"] for x in r])
    ax3.scatter(v, np.full(len(v), y), s=8 + 120*p/p.max(), color="#2a78d6", alpha=0.55, edgecolor="white", lw=0.6, zorder=3)
    for x in r:
        if x["name"] in lab_pick[lab]:
            nm = "Mexico City" if x["name"] == "Distrito Federal" and lab == "Mexico" else ("Brasília DF" if x["name"] == "Distrito Federal" else x["name"])
            o = lab_pick[lab][x["name"]]
            ax3.annotate(nm, (x["r_rel"], y), xytext=(0, 9*o), textcoords="offset points", ha="center", va="bottom" if o > 0 else "top", fontsize=6.8)
ax3.axvline(1, color="#52514e", lw=0.8)
ax3.set_xscale("log"); ax3.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax3.set_xlim(0.22, 5); ax3.set_xticks([0.25, 0.5, 1, 2, 4]); ax3.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax3.set_yticks(range(len(C))); ax3.set_yticklabels([c[0] for c in C][::-1]); ax3.set_ylim(-0.6, len(C) - 0.2)
ax3.set_xlabel("Light per resident relative to the national mean (log scale); dot area ∝ population")
ax3.set_title("c  Light per resident by region (GADM level 1, 10′ grid)", loc="left", fontsize=9, fontweight="bold")
for a in (ax, ax2, ax3):
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.grid(axis="x", color="#e4e3df", lw=0.6); a.set_axisbelow(True)
fig.text(0.01, 0.005, "No masking of light north of 65° N. Lengths of vector components; D = A + Bb + Bw + C holds as vectors. "
         "C (within cells) ≤ 0.4 km, not shown.", fontsize=6.6, color="#52514e")
for ext in ("pdf", "png"): fig.savefig(O / f"fig3_national_bias.{ext}", dpi=300, bbox_inches="tight")
d = J("BRA_0.1666667")["light_per_capita_by_region_rel"]; print([(x["name"], x["r_rel"]) for x in d if x["name"] in ("Rio de Janeiro", "São Paulo")])
