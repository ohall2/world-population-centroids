"""Fig. 4 (värden utanför axlarna läggs på kanten och räknas i panelen): mätt mot sann
förflyttning per kommun (nivå 2).
x = den sanna förflyttningens längd; y = ljusets förflyttning projicerad på den sanna
riktningen (negativ = fel riktning). Linjen y = x är perfekt mätning.
Övre raden (USA, Mexiko, Europa) och Kanada på ±4 km; Japan på ±1 km, eftersom dess
förflyttningar är mycket mindre (median 0,06 km mot 0,13–0,54 km).
Läser docs/figures/data/unit_shifts.csv.gz (tools/figures/unit_shifts.py)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from common import COL, INK, INK2, SHADE, style
R = Path(__file__).resolve().parents[2]; O = R / "docs/figures"
d = pd.read_csv(O / "data/unit_shifts.csv.gz")
d["true"] = np.hypot(d.true_e, d.true_n)
d["proj"] = (d.meas_e*d.true_e + d.meas_n*d.true_n) / d.true.where(d.true > 0)
cases = [("census", 2000, 2020, "a  USA, 2000–2020", 4.0),
         ("inegi", 2010, 2020, "b  Mexico, 2010–2020", 4.0),
         ("geostat", 2011, 2021, "c  Europe, 2011–2021", 4.0),
         ("statcan", 2016, 2021, "d  Canada, 2016–2021", 4.0),
         ("estat_jp", 2015, 2020, "e  Japan, 2015–2020 (note scale)", 1.0)]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, axs = plt.subplots(2, 3, figsize=(7.2, 6.0))
for ax, (ds, a, b, lab, L) in zip(axs.flat, cases):
    s = d[(d.dataset == ds) & (d["from"] == a) & (d.to == b)].dropna(subset=["proj"])
    ax.axhspan(-100, 0, color=SHADE, zorder=0)
    ax.scatter(s.true.clip(upper=L), s.proj.clip(-L, L), s=4 if len(s) < 500 else 3, color=COL[ds],
               alpha=0.6 if len(s) < 500 else 0.35, lw=0, rasterized=True, zorder=2)
    ax.plot([0, L], [0, L], color=INK, lw=1, zorder=3)
    ax.axhline(0, color=INK2, lw=0.8, zorder=1)
    ax.set_xlim(0, L); ax.set_ylim(-L, L)
    out = ((s.true > L) | (s.proj.abs() > L)).sum()
    right = (s.proj > 0).mean()
    ax.set_title(lab, loc="left", fontsize=8.5, fontweight="bold")
    ax.text(0.97, 0.04, f"n = {len(s):,}; {out} outside\nright direction {right:.0%}".replace(",", " "),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.text(0.825*L, 0.94*L, "y = x", fontsize=7.5, ha="right")
    style(ax, grid=None)
    if L == 1.0:
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticklabels(["−1", "−0.5", "0", "0.5", "1"])
for ax in axs[:, 0]:
    ax.set_ylabel("Light-measured shift along the\ntrue direction (km)")
for ax in (axs[0, 1], axs[0, 2]):
    ax.set_yticklabels([])
ax = axs[1, 2]; ax.set_axis_off()
ax.text(0.0, 0.5, "Shaded: the measured shift\npoints the wrong way.\n\nPanels a–d share the\n±4 km scale; panel e is\ndrawn at ±1 km.\n\nValues beyond the axes\nare drawn at the edge\nand counted.",
        transform=ax.transAxes, va="center", fontsize=7.5, color=INK2)
fig.supxlabel("True shift of the census population centroid (km); counties and municipalities (level 2)", fontsize=8.5, y=0.015)
fig.tight_layout(rect=(0, 0.025, 1, 1))
for ext in ("pdf", "png"): fig.savefig(O / f"fig4_measured_vs_true.{ext}", dpi=300)
print("ok")
