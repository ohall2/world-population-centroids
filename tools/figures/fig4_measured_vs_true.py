"""Fig. 4 (värden utanför ±4 km läggs på kanten och räknas i texten): mätt mot sann förflyttning per kommun (nivå 2).
x = den sanna förflyttningens längd; y = ljusets förflyttning projicerad på den sanna
riktningen (negativ = fel riktning). Linjen y = x är perfekt mätning.
Läser docs/figures/data/unit_shifts.csv.gz (tools/figures/unit_shifts.py)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R = Path(__file__).resolve().parents[2]; O = R / "docs/figures"
d = pd.read_csv(O / "data/unit_shifts.csv.gz")
d["true"] = np.hypot(d.true_e, d.true_n)
d["proj"] = (d.meas_e*d.true_e + d.meas_n*d.true_n) / d.true.where(d.true > 0)
cases = [("census", 2000, 2020, "USA, 2000–2020", "#2a78d6"),
         ("inegi", 2010, 2020, "Mexico, 2010–2020", "#eb6834"),
         ("geostat", 2011, 2021, "Europe, 2011–2021", "#1baf7a")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), sharey=True)
L = 4.0
for ax, (ds, a, b, lab, c) in zip(axs, cases):
    s = d[(d.dataset == ds) & (d["from"] == a) & (d.to == b)].dropna(subset=["proj"])
    ax.axhspan(-100, 0, color="#f0efec", zorder=0)
    ax.scatter(s.true.clip(upper=L), s.proj.clip(-L, L), s=3, color=c, alpha=0.35, lw=0, rasterized=True, zorder=2)
    ax.plot([0, L], [0, L], color="#0b0b0b", lw=1, zorder=3)
    ax.axhline(0, color="#52514e", lw=0.8, zorder=1)
    ax.set_xlim(0, L); ax.set_ylim(-L, L)
    out = ((s.true > L) | (s.proj.abs() > L)).sum()
    right = (s.proj > 0).mean()
    ax.set_title(lab, loc="left", fontsize=8.5, fontweight="bold")
    ax.text(0.97, 0.04, f"n = {len(s):,}; {out} outside\nright direction {right:.0%}".replace(",", " "),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.text(3.3, 3.75, "y = x", fontsize=7.5, ha="right")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0].set_ylabel("Light-measured shift along the\ntrue direction (km)")
fig.supxlabel("True shift of the census population centroid (km); counties and municipalities (level 2)", fontsize=8.5, y=0.02)
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("pdf", "png"): fig.savefig(O / f"fig4_measured_vs_true.{ext}", dpi=300)
print("ok")
