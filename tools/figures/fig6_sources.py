"""Fig. 6: jämförelse av källor på nivå 2, senaste folkräkningsår per facit.
(a) medianfel (km) mot facit per källa och facit, logaritmisk axel;
(b) andel enheter där källan ligger närmare facit än ljuset, med 95 %-blockbootstrap.
Läser docs/extended/<facit>_benchmark_{position,pairs}.csv (USA, Europa, Mexiko,
Brasilien) och docs/extended/C_{canada,japan}_{position,benchmark_pairs}.csv.
Positionsfilerna avser enheterna som är gemensamma för alla källor (B, selection.csv)."""
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from common import COL, MK, NAME, ORDER, LATEST, INK2, SHADE, style
R = Path(__file__).resolve().parents[2]; E = R / "docs/extended"; O = R / "docs/figures"
C = {"statcan": "canada", "estat_jp": "japan"}

def load(kind):
    out = []
    for ds in ORDER:
        if ds in C:
            d = pd.read_csv(E / f"C_{C[ds]}_{'position' if kind == 'position' else 'benchmark_pairs'}.csv")
        else:
            d = pd.read_csv(E / f"{ds}_benchmark_{kind}.csv")
        d = d[(d.level == 2) & (d.year == LATEST[ds])].assign(dataset=ds)
        assert len(d), (ds, kind)
        out.append(d)
    return pd.concat(out)

pos, pairs = load("position"), load("pairs")
SRC = {"light": "Nighttime light", "ghsbuilts": "Built-up surface\n(GHS-BUILT-S)",
       "ghspop": "GHS-POP*", "worldpop": "WorldPop*", "geometric": "Geometric centre"}
CENSUS = ("ghspop", "worldpop")   # innehåller folkräkningstal
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 4.2), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.6})
off = {ds: (i - 2.5) * 0.12 for i, ds in enumerate(ORDER)}

def rows(ax, srcs):
    for i, s in enumerate(srcs):
        y = len(srcs) - 1 - i
        if s in CENSUS: ax.axhspan(y - 0.45, y + 0.45, color=SHADE, zorder=0, lw=0)
    ax.set_yticks(range(len(srcs))); ax.set_ylim(-0.6, len(srcs) - 0.4)

s1 = list(SRC); rows(a1, s1); a1.set_yticklabels([SRC[s] for s in s1][::-1])
for i, s in enumerate(s1):
    for ds in ORDER:
        r = pos[(pos.source == s) & (pos.dataset == ds)].iloc[0]
        y = len(s1) - 1 - i - off[ds]
        a1.plot([r.err_med_lo, r.err_med_hi], [y, y], color=COL[ds], lw=1, zorder=2)
        a1.plot(r.err_med, y, marker=MK[ds], color=COL[ds], ms=5.5, mec="white", mew=0.7, ls="", zorder=3)
a1.set_xscale("log"); a1.set_xlim(0.15, 40)
a1.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20]); a1.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
a1.set_xlabel("Median distance to census centroid (km, log scale)")
a1.set_title("a  Location error by source", loc="left", fontsize=9, fontweight="bold")
style(a1, grid="x")

s2 = [s for s in SRC if s != "light"]; rows(a2, s2); a2.set_yticklabels([SRC[s] for s in s2][::-1])
a2.axvline(0.5, color=INK2, lw=0.8, zorder=1)
for i, s in enumerate(s2):
    for ds in ORDER:
        r = pairs[(pairs.source == s) & (pairs.dataset == ds)].iloc[0]
        y = len(s2) - 1 - i - off[ds]
        a2.plot([r.closer_lo, r.closer_hi], [y, y], color=COL[ds], lw=1, zorder=2)
        a2.plot(r.closer_share, y, marker=MK[ds], color=COL[ds], ms=5.5, mec="white", mew=0.7, ls="", zorder=3)
a2.set_xlim(0, 1); a2.set_xticks([0, 0.25, 0.5, 0.75, 1]); a2.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
a2.set_xlabel("Share of units where the source is\ncloser to the census centroid than light")
a2.set_title("b  Paired with light", loc="left", fontsize=9, fontweight="bold")
a2.text(0.51, -0.55, "50 %", fontsize=7.5, color=INK2, va="bottom")
style(a2, grid="x")
h = [Line2D([], [], marker=MK[ds], color=COL[ds], ls="", ms=6, mec="white", mew=0.7,
            label=f"{NAME[ds]} {LATEST[ds]}") for ds in ORDER]
fig.legend(handles=h, frameon=False, ncol=6, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, 0.075),
           handletextpad=0.2, columnspacing=1.0)
fig.text(0.01, 0.005, "* Contain census counts (shaded rows); WorldPop also uses nighttime lights. Second-level units, latest census year;\n"
         "units common to all sources. Lines: 95 % block-bootstrap intervals (blocks: first-level units).",
         fontsize=6.8, color=INK2)
fig.subplots_adjust(left=0.17, right=0.98, top=0.93, bottom=0.27)
for ext in ("pdf", "png"): fig.savefig(O / f"fig6_sources.{ext}", dpi=300)
print(pos.pivot(index="dataset", columns="source", values="err_med").round(2).to_string())
print(pairs.pivot(index="dataset", columns="source", values="closer_share").round(2).to_string())
