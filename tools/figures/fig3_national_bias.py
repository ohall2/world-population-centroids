"""Fig. 3: varför landscentroiderna ligger fel.
(a) Uppdelning av avvikelsen D = A + Bb + Bw + C (längder i km), rutnät 10′, utan maskning.
(b) Mexiko: delarnas längd beroende på rutnätet (1′–30′).
(Ljus per invånare per region visas i Fig. 1.)
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
fig = plt.figure(figsize=(7.2, 3.9))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.35, bottom=0.3, top=0.92, left=0.12, right=0.98)
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
fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 0.06), ncol=4)
ax2 = fig.add_subplot(gs[0, 1])
gr = ["0.0166667", "0.0333333", "0.0833333", "0.1666667", "0.5"]; gl = [1, 2, 5, 10, 30]
for key, pl, col in parts[1:]:
    ax2.plot(gl, [J(f"MEX_{g}")[key]["length_km"] for g in gr], marker="o", ms=4, lw=2, color=col, label=pl)
ax2.axhline(J("MEX_0.1666667")["D_total"]["length_km"], color="#52514e", lw=1, ls=(0, (4, 3)))
ax2.text(30, 134, "D = 131 km", ha="right", fontsize=7, color="#52514e")
ax2.set_xscale("log"); ax2.set_xticks(gl); ax2.set_xticklabels([f"{g}′" for g in gl])
ax2.set_xlabel("Grid cell size"); ax2.set_ylabel("km"); ax2.set_ylim(0, 145)
ax2.set_title("b  Mexico: split by grid size", loc="left", fontsize=9, fontweight="bold")
for a in (ax, ax2):
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.grid(axis="x", color="#e4e3df", lw=0.6); a.set_axisbelow(True)
fig.text(0.01, 0.005, "No masking of light north of 65° N. Lengths of vector components; D = A + Bb + Bw + C holds as vectors. "
         "C (within cells) ≤ 0.4 km, not shown.", fontsize=6.6, color="#52514e")
for ext in ("pdf", "png"): fig.savefig(O / f"fig3_national_bias.{ext}", dpi=300)
print("ok")
