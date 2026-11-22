"""Fig. 1: lägesfel per enhet (ljusets mot befolkningens centroid), per nivå och facit.
Läser docs/figures/data/unit_errors.csv.gz (från tools/figures/unit_errors.py).
Senaste facitåret per system. USA nivå 1 delas i delstater och Puerto Ricos municipios.
Färger och markörer per facit ur tools/figures/common.py."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from common import COL, LATEST
R = Path(__file__).resolve().parents[2]; O = R / "docs/figures"
d = pd.read_csv(O / "data/unit_errors.csv.gz")
fmtkm = lambda v: f"{v:.2f}" if v < 10 else f"{v:.0f}"
d = d[d.year == d.dataset.map(LATEST)]
def grp(r):
    if r.dataset == "census":
        if r.level == 1: return "USA states" if r.gid0 == "USA" else "Puerto Rico municipios"
        return ("USA" if r.gid0 == "USA" else "Puerto Rico") if r.level == 0 else "USA counties"
    return {"geostat": "Europe (GEOSTAT)", "inegi": "Mexico", "ibge": "Brazil",
            "statcan": "Canada", "estat_jp": "Japan"}[r.dataset]
d["grp"] = d.apply(grp, axis=1)
col = {"USA": COL["census"], "Puerto": COL["census"], "Europe": COL["geostat"], "Mexico": COL["inegi"],
       "Brazil": COL["ibge"], "Canada": COL["statcan"], "Japan": COL["estat_jp"]}
rest = ["Europe (GEOSTAT)", "Mexico", "Brazil", "Canada", "Japan"]
order = {0: ["USA", "Puerto Rico"] + rest, 1: ["USA states", "Puerto Rico municipios"] + rest,
         2: ["USA counties"] + rest}
title = {0: "Level 0 (countries)", 1: "Level 1 (states, regions)", 2: "Level 2 (counties, municipalities)"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, axs = plt.subplots(3, 1, figsize=(6.4, 7.4), sharex=True,
                        gridspec_kw={"height_ratios": [7, 7, 6]})
rng = np.random.default_rng(1)
for ax, lv in zip(axs, (0, 1, 2)):
    labs = order[lv]
    for i, g in enumerate(labs):
        e = d[(d.level == lv) & (d.grp == g)].err_km.values
        c = col[g.split()[0]]; y = len(labs) - 1 - i
        if len(e) >= 20:
            q10, q25, q50, q75, q90 = np.percentile(e, [10, 25, 50, 75, 90])
            ax.plot([q10, q90], [y, y], color=c, lw=1.2, solid_capstyle="butt")
            ax.add_patch(plt.Rectangle((q25, y-0.22), q75-q25, 0.44, color=c, alpha=0.35, lw=0))
            ax.plot([q50, q50], [y-0.3, y+0.3], color=c, lw=2.4)
            ax.text(q90*1.15, y, f"{fmtkm(q50)} km  (n = {len(e):,})".replace(",", " "), va="center", fontsize=7.5)
        else:
            ax.scatter(e, y + rng.uniform(-0.15, 0.15, len(e)), s=12, color=c, edgecolor="white", lw=0.5, zorder=3)
            ax.text(e.max()*1.25, y, (f"{fmtkm(e[0])} km" if len(e) == 1 else f"{fmtkm(np.median(e))} km  (n = {len(e)})"), va="center", fontsize=7.5)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs[::-1]); ax.set_ylim(-0.6, len(labs)-0.4)
    ax.set_title(title[lv], loc="left", fontsize=9, fontweight="bold")
    ax.set_xscale("log"); ax.grid(axis="x", color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axs[-1].set_xlim(0.03, 2000)
axs[-1].set_xticks([0.1, 1, 10, 100, 1000]); axs[-1].set_xticklabels(["0.1", "1", "10", "100", "1000"])
axs[-1].set_xlabel("Distance between light-weighted and census population centroid (km, log scale)")
fig.text(0.01, 0.004, "Bars: median; boxes: 25–75th percentile; lines: 10–90th percentile. Groups with n < 20 shown as points.\n"
         "Latest census year: USA 2020, Europe 2021, Mexico 2020, Brazil 2022, Canada 2021, Japan 2020.", fontsize=6.8, color="#52514e")
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("pdf", "png"): fig.savefig(O / f"fig1_location_error.{ext}", dpi=300)
print(d.groupby(["level", "grp"]).err_km.agg(["size", "median"]).round(2).to_string())
