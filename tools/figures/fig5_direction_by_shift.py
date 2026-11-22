"""Fig. 5: andel kommuner (nivå 2) där ljusets förflyttning har rätt riktning,
efter storleken på den sanna förflyttningen. USA, Mexiko och Europa ur
docs/change_detection/*_bins.csv; Kanada och Japan, som saknar klassfiler, klassas här
ur docs/figures/data/unit_shifts.csv.gz med samma klassgränser (common.direction_bins).
Klasser med n < 20 visas inte. Felstaplar: Wilson 95 % för en andel."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from common import COL, MK, INK2, SHADE, wilson, direction_bins, style
R = Path(__file__).resolve().parents[2]
D = R / "docs/change_detection"; O = R / "docs/figures"; O.mkdir(exist_ok=True)

cases = [("census", 2000, 2020, "USA 2000–2020"), ("inegi", 2010, 2020, "Mexico 2010–2020"),
         ("geostat", 2011, 2021, "Europe 2011–2021"), ("statcan", 2016, 2021, "Canada 2016–2021"),
         ("estat_jp", 2015, 2020, "Japan 2015–2020")]
SH = pd.read_csv(O / "data/unit_shifts.csv.gz")

def bins(ds, f, t):
    if (D / f"{ds}_change_detection_bins.csv").exists():
        b = pd.read_csv(D / f"{ds}_change_detection_bins.csv")
        return b[(b.series == "release") & (b.level == 2) & (b["from"] == f) & (b["to"] == t)]
    s = SH[(SH.dataset == ds) & (SH["from"] == f) & (SH.to == t)]
    return direction_bins(s[["true_e", "true_n"]].values, s[["meas_e", "meas_n"]].values)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig, ax = plt.subplots(figsize=(6.0, 3.8))
ax.axhspan(0, 0.5, color=SHADE, zorder=0)
ax.axhline(0.5, color=INK2, lw=0.8, zorder=1)
ax.axhline(0.8, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.text(0.021, 0.49, "chance", ha="left", va="top", fontsize=8, color=INK2)
ax.text(24, 0.805, "80 %", ha="right", va="bottom", fontsize=8, color=INK2)
for ds, f, t, lab in cases:
    b = bins(ds, f, t); b = b[b.n >= 20]
    lo, hi = wilson(b.dir_right.values, b.n.values)
    ax.errorbar(b.true_med, b.dir_right, yerr=[b.dir_right-lo, hi-b.dir_right], color=COL[ds],
                marker=MK[ds], ms=6, lw=2, capsize=0, elinewidth=1, label=lab, zorder=3,
                markeredgecolor="white", markeredgewidth=1)
    print(lab, b[["true_move_km", "n", "dir_right"]].round(3).values.tolist())
ax.set_xscale("log"); ax.set_xlim(0.02, 25); ax.set_ylim(0.25, 1.0)
ax.set_xticks([0.025, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10]); ax.set_xticklabels(["0.025", "0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10"])
ax.set_xlabel("True shift of the census centroid (km, median of class; log scale)")
ax.set_ylabel("Share of municipalities with\nmeasured shift in the right direction")
style(ax)
ax.legend(frameon=False, loc="lower right", fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(O / f"fig5_direction_by_shift.{ext}", dpi=300)
print("ok")
