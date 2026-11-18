"""Fig. 1: ljus per invånare per region (GADM nivå 1, rutnät 10′) och förskjutningen av
landets ljuscentroid från befolkningens centroid. Mexiko, Brasilien, Sverige, Norge, Finland.
Läser docs/check_identity/id_*.json och GADM nivå 1 (förenklad) som argument
(standard data/derived/check/adm1_simpl.gpkg)."""
import sys, json, numpy as np, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from pathlib import Path
R = Path(__file__).resolve().parents[2]; I = R / "docs/check_identity"; O = R / "docs/figures"
G = gpd.read_file(sys.argv[1] if len(sys.argv) > 1 else R / "data/derived/check/adm1_simpl.gpkg")
J = lambda k: json.load(open(I / f"id_{k}.json"))
cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#a9c8ef", "#f0efec", "#f6b89d", "#eb6834"])
norm = LogNorm(0.25, 4)
EXAG = 3
panels = [("a  Mexico", ["MEX_0.1666667"]), ("b  Brazil", ["BRA_0.1666667"]),
          ("c  Sweden, Norway, Finland", ["SWE_none", "NOR_none", "FIN_none"])]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, axs = plt.subplots(1, 3, figsize=(7.4, 4.2), gridspec_kw={"width_ratios": [1.25, 1.05, 0.75]})
for i, (ax, (title, keys)) in enumerate(zip(axs, panels)):
    fig.text([0.01, 0.40, 0.73][i], 0.95, title, fontsize=9, fontweight="bold", va="top")
    lat0 = []
    for k in keys:
        d = J(k); c = k[:3]
        rr = {x["name"]: x["r_rel"] for x in d["light_per_capita_by_region_rel"]}
        g = G[G.GID_0 == c].copy(); g["r"] = g.NAME_1.map(rr)
        g.plot(ax=ax, column="r", cmap=cmap, norm=norm, edgecolor="white", linewidth=0.3)
        (px, py), (lx, ly) = d["pop_centroid"], d["light_centroid"]
        lx, ly = px + EXAG*(lx - px), py + EXAG*(ly - py)   # pilen förlängd EXAG gånger
        ax.annotate("", xy=(lx, ly), xytext=(px, py),
                    arrowprops=dict(arrowstyle="-|>", color="#0b0b0b", lw=1.4, mutation_scale=9, shrinkA=0, shrinkB=0), zorder=5)
        ax.scatter([px], [py], s=22, color="white", edgecolor="#0b0b0b", lw=1.1, zorder=6)
        off = {"BRA": (8, -4), "SWE": (6, -2), "FIN": (-4, 8)}.get(c, (5, 3))
        ax.annotate(f"{d['D_total']['length_km']:.0f} km", (px, py) if c in ("BRA", "SWE") else (lx, ly), xytext=off, textcoords="offset points",
                    fontsize=7.5, fontweight="bold")
        lat0.append(py)
    ax.set_aspect(1 / np.cos(np.radians(np.mean(lat0))))
    ax.set_axis_off()
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cax = fig.add_axes([0.25, 0.17, 0.5, 0.025])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0.25, 0.5, 1, 2, 4])
cb.ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"]); cb.ax.minorticks_off()
cb.set_label("Light per resident relative to the national mean (GADM level 1, 10′ grid)", fontsize=8)
fig.text(0.01, 0.005, "Arrows run from the census population centroid (white dot) towards the light-weighted centroid, drawn 3 times their true length;\n"
         "labels give the true distance. VIIRS 2020 (Mexico), 2022 (Brazil), 2021 (Sweden, Norway, Finland); no masking north of 65° N.", fontsize=6.6, color="#52514e")
fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.25, wspace=0.05)
for ext in ("pdf", "png"): fig.savefig(O / f"fig1_map.{ext}", dpi=300)
print("ok")
