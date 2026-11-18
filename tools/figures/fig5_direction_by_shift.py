"""Fig. 5: andel kommuner (nivå 2) där ljusets förflyttning har rätt riktning,
efter storleken på den sanna förflyttningen. Läser docs/change_detection/*_bins.csv.
Klasser med n < 20 visas inte. Felstaplar: Wilson 95 % för en andel."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R = Path(__file__).resolve().parents[2]
D = R / "docs/change_detection"; O = R / "docs/figures"; O.mkdir(exist_ok=True)

def wilson(p, n, z=1.96):
    c = (p + z*z/(2*n)) / (1 + z*z/n)
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n)
    return c-h, c+h

cases = [("census", 2000, 2020, "USA 2000–2020", "#2a78d6", "o"),
         ("inegi", 2010, 2020, "Mexico 2010–2020", "#eb6834", "s"),
         ("geostat", 2011, 2021, "Europe 2011–2021", "#1baf7a", "^")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig, ax = plt.subplots(figsize=(6.0, 3.8))
ax.axhspan(0, 0.5, color="#f0efec", zorder=0)
ax.axhline(0.5, color="#52514e", lw=0.8, zorder=1)
ax.axhline(0.8, color="#52514e", lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.text(24, 0.49, "chance", ha="right", va="top", fontsize=8, color="#52514e")
ax.text(24, 0.805, "80 %", ha="right", va="bottom", fontsize=8, color="#52514e")
for ds, f, t, lab, col, mk in cases:
    b = pd.read_csv(D / f"{ds}_change_detection_bins.csv")
    b = b[(b.series == "release") & (b.level == 2) & (b["from"] == f) & (b["to"] == t) & (b.n >= 20)]
    lo, hi = wilson(b.dir_right.values, b.n.values)
    ax.errorbar(b.true_med, b.dir_right, yerr=[b.dir_right-lo, hi-b.dir_right], color=col,
                marker=mk, ms=6, lw=2, capsize=0, elinewidth=1, label=lab, zorder=3,
                markeredgecolor="white", markeredgewidth=1)
ax.set_xscale("log"); ax.set_xlim(0.04, 25); ax.set_ylim(0.3, 1.0)
ax.set_xticks([0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10]); ax.set_xticklabels(["0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10"])
ax.set_xlabel("True shift of the census centroid (km, median of class; log scale)")
ax.set_ylabel("Share of municipalities with\nmeasured shift in the right direction")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left", fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(O / f"fig5_direction_by_shift.{ext}", dpi=300)
print("ok")
