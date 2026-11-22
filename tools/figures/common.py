"""Gemensamt för figurerna: färger, markörer och namn per facit, Wilson-intervall
och storleksklasser för den sanna förflyttningen.

Färger: USA, Mexiko och Europa som tidigare (referenspalettens plats 1–3).
Brasilien, Kanada och Japan valda så att varje par skiljer sig med ΔE ≥ 8,8
(OKLab × 100) under simulerad protanopi, deuteranopi och tritanopi (Machado m.fl.
2009, full grad) och ΔE ≥ 15 för normalt färgseende. Markörerna skiljer sig
dessutom per facit, så identiteten bärs inte av färgen ensam."""
import numpy as np, pandas as pd

COL = {"census": "#2a78d6", "inegi": "#eb6834", "geostat": "#1baf7a",
       "ibge": "#882255", "statcan": "#ddaa33", "estat_jp": "#3b1f6e"}
MK = {"census": "o", "inegi": "s", "geostat": "^", "ibge": "D", "statcan": "v", "estat_jp": "h"}
NAME = {"census": "USA", "inegi": "Mexico", "geostat": "Europe", "ibge": "Brazil",
        "statcan": "Canada", "estat_jp": "Japan"}
ORDER = ["census", "geostat", "inegi", "ibge", "statcan", "estat_jp"]
LATEST = {"census": 2020, "geostat": 2021, "inegi": 2020, "ibge": 2022, "statcan": 2021, "estat_jp": 2020}
INK, INK2, GRID, SHADE = "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"

# Klassgränser som i tools/change_detection.py (frysta), vänsterslutna.
EDGES = [0, 0.1, 0.25, 0.5, 1, 2, 5, 10, np.inf]


def wilson(p, n, z=1.96):
    p, n = np.asarray(p, float), np.asarray(n, float)
    c = (p + z*z/(2*n)) / (1 + z*z/n)
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n)
    return c-h, c+h


def direction_bins(t, e):
    """Andel rätt riktning (cos > 0) per storleksklass av |t|, som binningen i
    tools/change_detection.py. t, e: (n, 2) km öst/nord. Returnerar DataFrame med
    true_move_km, n, dir_right, true_med, err_med."""
    t, e = np.asarray(t, float), np.asarray(e, float)
    tn = np.hypot(*t.T); err = np.hypot(*(e - t).T)
    cos = (t * e).sum(1) / np.maximum(tn * np.hypot(*e.T), 1e-12)
    cl = pd.cut(tn, EDGES, right=False)
    rows = [{"true_move_km": str(c), "n": len(g), "dir_right": g.ok.mean(),
             "true_med": g.tn.median(), "err_med": g.err.median()}
            for c, g in pd.DataFrame({"cl": cl, "ok": cos > 0, "err": err, "tn": tn}).groupby("cl", observed=True)]
    return pd.DataFrame(rows)


def style(ax, grid="y"):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if grid: ax.grid(axis=grid, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
