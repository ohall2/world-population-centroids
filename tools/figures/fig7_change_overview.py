"""Fig. 7: förändring i översikt, alla facit och perioder, ordnade efter sensorperiod
(DMSP, över skarven 2012→2013, bara VIIRS).
(a) andel rätt riktning på nivå 2 med 95 %-blockbootstrap; linjer vid 0,50 (slump) och
0,65 (tolkningsregelns första rad, planen §7).
(b) nivå 1: ljusets centroid, hybrid-ljus och hybrid-geom (A1). Hybridskattningarnas
förändring är inte beräknad för Kanada och Japan (tools/canada_confirm.py redovisar
hybriden bara för läge), så där visas bara ljuset.
Läser docs/change_detection/<facit>_change_detection.csv (serie release),
docs/extended/<facit>_hybrid_change.csv och docs/extended/C_{canada,japan}_change.csv."""
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from common import COL, MK, NAME, INK, INK2, SHADE, style
R = Path(__file__).resolve().parents[2]; E = R / "docs/extended"; D = R / "docs/change_detection"; O = R / "docs/figures"
CASES = [("DMSP", "census", 2000, 2010),
         ("Across the join", "census", 2010, 2020), ("Across the join", "census", 2000, 2020),
         ("Across the join", "inegi", 2010, 2020), ("Across the join", "geostat", 2011, 2021),
         ("VIIRS", "statcan", 2016, 2021), ("VIIRS", "estat_jp", 2015, 2020)]
C = {"statcan": "canada", "estat_jp": "japan"}
q = lambda d, a, b, **k: d[(d["from"] == a) & (d.to == b) & pd.concat([d[c] == v for c, v in k.items()], axis=1).all(axis=1)].iloc[0]

def lvl2(ds, a, b):
    if ds in C: return q(pd.read_csv(E / f"C_{C[ds]}_change.csv"), a, b, level=2)
    return q(pd.read_csv(D / f"{ds}_change_detection.csv"), a, b, series="release", level=2)

def lvl1(ds, a, b):
    if ds in C: return {"light": q(pd.read_csv(E / f"C_{C[ds]}_change.csv"), a, b, level=1)}
    h = pd.read_csv(E / f"{ds}_hybrid_change.csv")
    return {s: q(h, a, b, series=s, level=1) for s in ("light", "hybrid_light", "hybrid_geom")}

EST = {"light": ("Light-weighted centroid", "o", "full"),
       "hybrid_light": ("Hybrid, level-2 light centroids", "s", "none"),
       "hybrid_geom": ("Hybrid, level-2 geometric centres", "^", "none")}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.8, 6.6), sharex=True, gridspec_kw={"height_ratios": [1, 1.1]})
x = range(len(CASES))
rows = []
for ax in (a1, a2):
    ax.axhspan(0, 0.5, color=SHADE, zorder=0)
    ax.axhline(0.5, color=INK2, lw=0.8, zorder=1)
    ax.axhline(0.65, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for xb in (0.5, 4.5): ax.axvline(xb, color=INK2, lw=0.6, ls=":", zorder=1)
    style(ax)
for i, (per, ds, a, b) in enumerate(CASES):
    r = lvl2(ds, a, b)
    a1.errorbar(i, r.dir_right, yerr=[[r.dir_right - r.dir_right_lo], [r.dir_right_hi - r.dir_right]], color=COL[ds],
                marker=MK[ds], ms=7, elinewidth=1.4, capsize=0, mec="white", mew=0.8, zorder=3)
    a1.text(i + 0.12, r.dir_right, f"{r.dir_right:.2f}", va="center", fontsize=7, color=INK)
    rows.append({"level": 2, "dataset": ds, "from": a, "to": b, "series": "light", "n": int(r.n),
                 "dir_right": r.dir_right, "lo": r.dir_right_lo, "hi": r.dir_right_hi})
    l1 = lvl1(ds, a, b)
    for j, (s, rr) in enumerate(l1.items()):
        lab, mk, fs = EST[s]; xo = i + (j - 1) * 0.22 if len(l1) > 1 else i
        a2.errorbar(xo, rr.dir_right, yerr=[[rr.dir_right - rr.dir_right_lo], [rr.dir_right_hi - rr.dir_right]],
                    color=COL[ds], marker=mk, ms=6, fillstyle=fs, mfc=COL[ds] if fs == "full" else "white",
                    elinewidth=1.2, capsize=0, mew=1.2, zorder=3)
        rows.append({"level": 1, "dataset": ds, "from": a, "to": b, "series": s, "n": int(rr.n),
                     "dir_right": rr.dir_right, "lo": rr.dir_right_lo, "hi": rr.dir_right_hi})
a1.set_ylim(0.4, 0.75); a2.set_ylim(0.0, 1.04)
a1.text(-0.45, 0.652, "0.65", fontsize=7, color=INK2, va="bottom", ha="left")
a1.text(-0.45, 0.497, "chance", fontsize=7, color=INK2, va="top", ha="left")
for xc, lab in ((0, "DMSP"), (2.5, "Across the sensor join"), (5.5, "VIIRS only")):
    a1.text(xc, 0.745, lab, ha="center", va="top", fontsize=8, color=INK2, fontweight="bold")
a1.set_ylabel("Share in the right direction")
a1.set_title("a  Level 2 (counties, municipalities)", loc="left", fontsize=9, fontweight="bold")
a2.set_ylabel("Share in the right direction")
a2.set_title("b  Level 1 (states, regions): light and hybrid estimates", loc="left", fontsize=9, fontweight="bold")
a2.set_xticks(list(x)); a2.set_xticklabels([f"{NAME[ds]}\n{a}–{b}" for _, ds, a, b in CASES], fontsize=8)
h = [Line2D([], [], color=INK2, marker=mk, ls="", ms=6, fillstyle=fs, mfc=INK2 if fs == "full" else "white", mew=1.2,
            label=lab) for lab, mk, fs in EST.values()]
fig.legend(handles=h, frameon=False, fontsize=7.5, loc="lower center", ncol=3, handletextpad=0.3, columnspacing=1.2,
           bbox_to_anchor=(0.5, 0.045))
fig.text(0.01, 0.005, "Lines: 95 % block-bootstrap intervals. Dashed: 0.65, the first row of the pre-specified rule for Canada and Japan.\n"
         "Hybrid estimates of change were not computed for Canada and Japan; only light is shown there. Colours as in Figures 1, 4–6.",
         fontsize=6.6, color=INK2)
fig.tight_layout(rect=(0, 0.075, 1, 1))
for ext in ("pdf", "png"): fig.savefig(O / f"fig7_change_overview.{ext}", dpi=300)
print(pd.DataFrame(rows).round(3).to_string(index=False))
