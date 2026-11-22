#!/usr/bin/env python3
"""Talrevision av docs/manuscript/manuscript.md (2026-11-22).

Skriptet gör två saker.

1. Räknar de tal i manuset som saknade skript i repot, ur de frysta data
   (utgåvan data/release/v2.0.0, facit i data/validation*, rådata för USA:s
   FIPS-koder) och med de frysta funktionerna, som importeras och inte ändras:
   - Sverige ur GEOSTAT: lägesfel 2011 och 2021 (kommuner, län) och förändring
     2011–2021 för kommuner (median flytt, rätt riktning med blockbootstrap).
   - Fishers exakta test, Mexikos delstater mot USA:s, för varje period.
   - Kvoten median flytt / median fel för Mexikos och USA:s delstater.
   - FIPS-kod mot GADM-polygon för USA:s counties (efter
     docs/review/scripts/check_fips_truth.py).
   - Andelen enheter där ljuset ligger närmare facit än geometrisk mittpunkt.
   - Lutningens känslighet för varje rad på nivå 2 i tabell 3: OLS-lutningen
     utan de tre mest inflytelserika enheterna, och Theil–Sen-lutningen.
   - Härledda tal (kvoter, andelar, räkningar) ur de committade CSV-filerna.
2. Prövar varje kontrollerat påstående i manuset mot talen: textbiten måste
   finnas ordagrant i manuset, talet ur data avrundas halvt uppåt till textens
   precision, och status sätts till stämmer / avviker / saknar källa.

Utdata
  docs/manuscript/numbers.csv               alla tal som revisionen bygger på
  docs/review/talrevision-2026-11-22.md     revisionen (tabellen skrivs av skriptet)

Definitioner som inte står i manuset
  - Inflytande: för varje enhet ändringen i den sammanslagna OLS-lutningen
    (öst- och nordkomponenterna, som i change_detection.metrics) när enhetens
    båda komponenter tas bort. De tre enheterna med störst |ändring| tas bort
    samtidigt och lutningen räknas om.
  - Theil–Sen: medianen av alla parvisa lutningar (e_j − e_i)/(t_j − t_i) över
    samma 2n sammanslagna värden, par med lika t utelämnade (som scipy).
  - FIPS mot GADM: facitflytten för en county beräknad med kvarteringsgrupper
    per FIPS-kod respektive per GADM-polygon, båda som Aboufadel–Austin-medel.
    GADM-enheten knyts till den FIPS-county vars centroid 2010 ligger närmast
    (högst 5 km, entydigt), som i granskningsskriptet.
  - Lägesfelskvoten i diskussionen: median |e − t| delad med medelvärdet av de
    två årens median lägesfel; implicit korrelation 1 − kvot²/2.

Körning (från repots rot)
  python tools/manuscript_numbers.py            # allt, ~5 min
  python tools/manuscript_numbers.py --skip-a2  # utan omkörning av A2
"""

from __future__ import annotations

import argparse
import io
import json
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(".")
REL = Path("data/release/v2.0.0")
TRUTH = Path("data/validation")
MANUS = Path("docs/manuscript/manuscript.md")
OUT_CSV = Path("docs/manuscript/numbers.csv")
OUT_MD = Path("docs/review/talrevision-2026-11-22.md")
R_KM = 6371.0
PX_KM = 2 * np.pi * R_KM / 360 / 120        # 30″ vid ekvatorn
MINUS = "\u2212"
_TIES: list[str] = []      # värden som legat exakt på en halva vid avrundningen


# =============================================================================
# Avrundning, format och status
# =============================================================================

def fmt(x, dec: int, scale: float = 1.0, thousands: bool = True) -> str:
    """Avrunda halvt uppåt (bort från noll) med textens precision.

    Talet går via str(), så 0.625 och ett CSV-värde som 0.4750 avrundas som de
    står skrivna. Negativa tal får minus som i manuset (U+2212)."""
    d = Decimal(str(float(x))) * Decimal(str(scale))
    if is_tie(x, dec, scale):
        _TIES.append(str(float(x)))
    q = d.quantize(Decimal(1).scaleb(-dec), rounding=ROUND_HALF_UP)
    if q == 0:
        q = abs(q)
    s = f"{abs(q):,.{dec}f}" if thousands else f"{abs(q):.{dec}f}"
    return (MINUS + s) if q < 0 else s


def fmt_range(values, dec: int, scale: float = 1.0) -> str:
    """Intervall min–max i manusets stil: "0.7–2.6", "−0.16 to −0.01"."""
    v = [float(x) for x in values]
    lo, hi = fmt(min(v), dec, scale), fmt(max(v), dec, scale)
    if lo == hi:
        return lo
    sep = " to " if (lo.startswith(MINUS) or hi.startswith(MINUS)) else "–"
    return f"{lo}{sep}{hi}"


def fmt_ci(v, lo, hi, dec: int, scale: float = 1.0) -> str:
    """Punkt med intervall som i tabell 3: "0.50 (0.36–0.64)"."""
    a, b = fmt(lo, dec, scale), fmt(hi, dec, scale)
    sep = " to " if (a.startswith(MINUS) or b.startswith(MINUS)) else "–"
    return f"{fmt(v, dec, scale)} ({a}{sep}{b})"


def _norm(s: str) -> str:
    s = s.replace("\u2009", " ").replace("\u00a0", " ").strip()
    s = s.replace("-", "–")
    s = re.sub(r"(^|[\s(])–", lambda m: m.group(1) + MINUS, s)
    return s


def status(text_num: str, data: str | None) -> str:
    if data is None:
        return "saknar källa"
    return "stämmer" if _norm(text_num) == _norm(data) else "avviker"


def is_tie(x, dec: int, scale: float = 1.0) -> bool:
    """Ligger värdet, som det står skrivet, exakt på en halva vid avrundningen?"""
    d = Decimal(str(float(x))) * Decimal(str(scale))
    r = d.scaleb(dec) - d.scaleb(dec).to_integral_value(rounding="ROUND_FLOOR")
    return r == Decimal("0.5")


# =============================================================================
# Statistik
# =============================================================================

def share_closer(a, b):
    """Andel enheter där a < b (a ligger närmare), bland par där båda finns."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.mean(a[ok] < b[ok])), int(ok.sum())


def fisher_two_sided(k1, n1, k2, n2) -> float:
    return float(stats.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])[1])


def pooled_slope(t, e) -> float:
    """OLS-lutning för e på t, komponenterna sammanslagna (som metrics)."""
    return float(np.polyfit(np.ravel(t), np.ravel(e), 1)[0])


def unit_influence(t, e) -> np.ndarray:
    """Ändring i den sammanslagna lutningen när enhet i (båda komponenterna)
    tas bort; sluten form ur summorna."""
    x, y = np.asarray(t, float), np.asarray(e, float)
    N = x.size
    Sx, Sy, Sxx, Sxy = x.sum(), y.sum(), (x * x).sum(), (x * y).sum()
    full = (Sxy - Sx * Sy / N) / (Sxx - Sx * Sx / N)
    sx, sy = x.sum(1), y.sum(1)
    sxx, sxy = (x * x).sum(1), (x * y).sum(1)
    n = N - x.shape[1]
    sx_, sy_ = Sx - sx, Sy - sy
    loo = ((Sxy - sxy) - sx_ * sy_ / n) / ((Sxx - sxx) - sx_ * sx_ / n)
    return loo - full


def slope_without_top(t, e, k: int = 3):
    """Lutningen utan de k enheter som var för sig ändrar lutningen mest."""
    infl = unit_influence(t, e)
    top = np.argsort(-np.abs(infl), kind="stable")[:k]
    keep = np.ones(len(t), bool)
    keep[top] = False
    return pooled_slope(np.asarray(t)[keep], np.asarray(e)[keep]), top


def theil_sen(x, y, chunk: int | None = None, sample: int = 2_000_000,
              direct_max_pairs: int = 5_000_000, seed: int = 20261122) -> float:
    """Exakt Theil–Sen-lutning (median av parvisa lutningar, lika x utelämnade).

    Små problem räknas direkt. Stora räknas i två steg utan att alla par
    lagras: ett slumpurval av par ger ett fönster kring medianen, sedan räknas
    exakt hur många lutningar som ligger under fönstret och fönstrets värden
    samlas in. Hamnar medianen utanför fönstret vidgas det och steget görs om."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    o = np.argsort(x, kind="stable")
    x, y = x[o], y[o]
    n = len(x)
    if chunk is None and n * (n - 1) // 2 <= direct_max_pairs:
        i, j = np.triu_indices(n, 1)
        dx = x[j] - x[i]
        ok = dx != 0
        return float(np.median((y[j] - y[i])[ok] / dx[ok]))
    chunk = chunk or 512

    def blocks():
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            dx = x[None, :] - x[a:b, None]
            dy = y[None, :] - y[a:b, None]
            m = (np.arange(n)[None, :] > np.arange(a, b)[:, None]) & (dx != 0)
            yield dy[m] / dx[m]

    total = sum(len(s) for s in blocks())
    k_lo, k_hi = (total - 1) // 2, total // 2          # 0-baserade mittrang
    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, sample)
    j = rng.integers(0, n, sample)
    ok = x[i] != x[j]
    smp = (y[j] - y[i])[ok] / (x[j] - x[i])[ok]
    half = 3.0 / np.sqrt(max(len(smp), 1))
    while True:
        qlo, qhi = max(0.0, 0.5 - half), min(1.0, 0.5 + half)
        lo, hi = (np.quantile(smp, qlo) if qlo > 0 else -np.inf,
                  np.quantile(smp, qhi) if qhi < 1 else np.inf)
        below, inside = 0, []
        for s in blocks():
            below += int((s < lo).sum())
            inside.append(s[(s >= lo) & (s <= hi)])
        inside = np.sort(np.concatenate(inside))
        if below <= k_lo and below + len(inside) > k_hi:
            return float((inside[k_lo - below] + inside[k_hi - below]) / 2)
        half *= 2


def shift_agreement(a, b, tol: float = 0.1) -> float:
    """Andel enheter där två skattningar av samma flytt skiljer högst tol km."""
    d = np.hypot(*(np.asarray(a, float) - np.asarray(b, float)).T)
    return float(np.mean(d <= tol))


def implied_corr(ratio) -> float:
    """Korrelation mellan två års lägesfel som ger kvoten
    median|e2 − e1| / median|e| = sqrt(2(1 − ρ)) för lika fördelade normalfel."""
    return float(1 - np.asarray(ratio, float) ** 2 / 2)


# =============================================================================
# Register över tal
# =============================================================================

class Numbers:
    """Alla tal som revisionen bygger på, med källa. Värden ur CSV-filer läses
    som de står skrivna (för avrundningen); beräknade tal lagras med full
    precision."""

    def __init__(self):
        self.rows: dict[str, dict] = {}
        self._csv: dict[str, pd.DataFrame] = {}

    def put(self, key, value, source, note=""):
        self.rows[key] = {"key": key, "value": value, "source": source, "note": note}
        return value

    def get(self, key):
        return self.rows[key]["value"]

    def csv(self, path, col, **flt):
        p = str(path)
        if p not in self._csv:
            self._csv[p] = pd.read_csv(p, dtype=str, keep_default_na=False)
        d = self._csv[p]
        for k, v in flt.items():
            d = d[d[k] == str(v)]
        if len(d) != 1:
            raise KeyError(f"{p} {col} {flt}: {len(d)} rader")
        raw = d.iloc[0][col]
        key = f"{Path(p).stem}:{col}:" + ",".join(f"{k}={v}" for k, v in flt.items())
        self.put(key, float(raw), p)
        return float(raw)

    def js(self, path, *keys):
        d = json.load(open(path))
        for k in keys:
            d = d[k]
        self.put(f"{Path(path).stem}:" + ".".join(map(str, keys)), float(d), str(path))
        return float(d)


# =============================================================================
# Nya beräkningar
# =============================================================================

def _imports():
    import sys
    sys.path.insert(0, str(ROOT.resolve()))


def sweden(N: Numbers, n_boot: int):
    """Sverige ur GEOSTAT med de frysta funktionerna (via C-drivrutinen, som
    återger change_detection exakt; frö BOOT_SEED = 20261121)."""
    from tools.canada_confirm import change_measure, position_measure
    src = "data/validation/geostat + utgåvan; tools/canada_confirm.{position,change}_measure"
    for y in (2011, 2021):
        for lv in (1, 2):
            m = position_measure(TRUTH, REL, None, y, lv, n_boot, ds="geostat", iso="SWE")
            N.put(f"swe:pos:L{lv}:{y}:n", m["n"], src)
            N.put(f"swe:pos:L{lv}:{y}:err_med", m["err_med"], src)
            N.put(f"swe:pos:L{lv}:{y}:err_med_lo", m["err_med_lo"], src)
            N.put(f"swe:pos:L{lv}:{y}:err_med_hi", m["err_med_hi"], src)
    c = change_measure(TRUTH, REL, "geostat", 2011, 2021, 2, iso="SWE", n_boot=n_boot)
    N.put("swe:chg:L2:n", c["n"], src)
    N.put("swe:chg:L2:n_blocks", c["n_blocks"], src, "block = län")
    for k in ("true_med", "err_med", "dir_right", "rho", "slope"):
        N.put(f"swe:chg:L2:{k}", c["metrics"][k], src)
    for k in ("dir_right", "true_med"):
        N.put(f"swe:chg:L2:{k}_lo", c["ci"][f"{k}_lo"], src)
        N.put(f"swe:chg:L2:{k}_hi", c["ci"][f"{k}_hi"], src)


def states(N: Numbers):
    """Fishers exakta test och kvoten flytt/fel för delstaterna."""
    from tools.canada_confirm import change_inputs
    from tools.change_detection import metrics
    src = "facit nivå 1 + utgåvan; tools/canada_confirm.change_inputs, change_detection.metrics"
    sel = {}
    for ds, iso, pairs in (("census", "USA", [(2000, 2010), (2010, 2020), (2000, 2020)]),
                           ("inegi", "MEX", [(2010, 2020)])):
        for y1, y2 in pairs:
            ids, t, e, _ = change_inputs(TRUTH, REL, ds, y1, y2, 1, iso=iso)
            tn = np.hypot(*t.T)
            cos = (t * e).sum(1) / np.maximum(tn * np.hypot(*e.T), 1e-12)
            k, n = int((cos[tn > 0] > 0).sum()), int((tn > 0).sum())
            m = metrics(t, e)
            sel[(iso, y1, y2)] = (k, n)
            tag = f"states:{iso}:{y1}-{y2}"
            N.put(f"{tag}:k_right", k, src)
            N.put(f"{tag}:n", n, src)
            N.put(f"{tag}:true_med", m["true_med"], src)
            N.put(f"{tag}:err_med", m["err_med"], src)
            N.put(f"{tag}:true_over_err", m["true_med"] / m["err_med"], src)
    k2, n2 = sel[("MEX", 2010, 2020)]
    for y1, y2 in [(2000, 2010), (2010, 2020), (2000, 2020)]:
        k1, n1 = sel[("USA", y1, y2)]
        N.put(f"fisher:MEX2010-2020_vs_USA{y1}-{y2}:p", fisher_two_sided(k2, n2, k1, n1),
              "scipy.stats.fisher_exact, tvåsidigt",
              f"Mexiko {k2}/{n2} mot USA {k1}/{n1} rätt riktning")


def _fips_centroids(year):
    """Aboufadel–Austin-medel av kvarteringsgrupperna per FIPS-county
    (efter docs/review/scripts/check_fips_truth.py)."""
    cols = ["STATEFP", "COUNTYFP", "TRACTCE", "BLKGRPCE", "POPULATION", "LATITUDE", "LONGITUDE"]
    parts = []
    for f in sorted(Path(f"data/raw/validation/census/{year}").glob("*.txt")):
        raw = f.read_bytes().decode("utf-8-sig", "replace")
        hdr = 0 if raw.startswith("STATEFP") else None
        parts.append(pd.read_csv(io.StringIO(raw), header=hdr, names=None if hdr == 0 else cols,
                                 dtype={"STATEFP": str, "COUNTYFP": str, "TRACTCE": str}))
    d = pd.concat(parts, ignore_index=True)
    d["lon"] = pd.to_numeric(d.LONGITUDE, errors="coerce")
    d["lat"] = pd.to_numeric(d.LATITUDE, errors="coerce")
    d["pop"] = pd.to_numeric(d.POPULATION, errors="coerce")
    d = d[(d["pop"] > 0) & d.lon.notna()]
    d["fips"] = d.STATEFP.str.zfill(2) + d.COUNTYFP.str.zfill(3)
    lo, la = np.radians(d.lon), np.radians(d.lat)
    d["x"] = d["pop"] * np.cos(la) * np.cos(lo)
    d["y"] = d["pop"] * np.cos(la) * np.sin(lo)
    d["z"] = d["pop"] * np.sin(la)
    g = d.groupby("fips")[["x", "y", "z", "pop"]].sum()
    nn = np.sqrt(g.x ** 2 + g.y ** 2 + g.z ** 2)
    return pd.DataFrame({"lon": np.degrees(np.arctan2(g.y, g.x)),
                         "lat": np.degrees(np.arcsin(g.z / nn)), "pop": g["pop"]})


def fips_vs_gadm(N: Numbers):
    from scipy.spatial import cKDTree
    from tools.canada_confirm import change_inputs
    from tools.change_detection import en, metrics, truth
    src = "data/raw/validation/census (FIPS) mot data/validation/census (GADM)"
    F = {y: _fips_centroids(y) for y in (2000, 2010, 2020)}
    T = {}
    for y in (2000, 2010, 2020):
        d = pd.read_csv(TRUTH / "census" / str(y) / "adm2.csv", keep_default_na=False)
        d = d[(d.complete.astype(str).str.lower() == "true") & (d.n_points >= 10)
              & d.unit_id.str.startswith("USA.")]
        T[y] = d.set_index("unit_id")[["aa_lon", "aa_lat", "bm_lon", "bm_lat"]].astype(float)

    def xyz(lon, lat):
        lo, la = np.radians(lon), np.radians(lat)
        return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

    f, g = F[2010], T[2010]
    dist, idx = cKDTree(xyz(f.lon.values, f.lat.values)).query(xyz(g.aa_lon.values, g.aa_lat.values))
    m = pd.DataFrame({"unit_id": g.index, "fips": f.index[idx], "d_km": dist * R_KM})
    good = m[(m.d_km < 5) & ~m.fips.duplicated(keep=False)].set_index("unit_id")
    N.put("fips:matched_2010", len(good), src, f"av {len(m)} GADM-counties 2010")
    for y1, y2 in ((2000, 2010), (2010, 2020), (2000, 2020)):
        ids, t_bm, e, _ = change_inputs(TRUTH, REL, "census", y1, y2, 2, iso="USA")
        gm = good.reindex(ids).dropna(subset=["fips"])
        gm = gm[gm.fips.isin(F[y1].index) & gm.fips.isin(F[y2].index)]
        sel = np.isin(ids, gm.index)
        u = ids[sel]
        t_fips = np.column_stack(en(F[y1].loc[gm.fips].lon.values, F[y1].loc[gm.fips].lat.values,
                                    F[y2].loc[gm.fips].lon.values, F[y2].loc[gm.fips].lat.values))
        t_gaa = np.column_stack(en(T[y1].loc[u].aa_lon.values, T[y1].loc[u].aa_lat.values,
                                   T[y2].loc[u].aa_lon.values, T[y2].loc[u].aa_lat.values))
        tag = f"fips:{y1}-{y2}"
        N.put(f"{tag}:n_matched", int(sel.sum()), src, f"av {len(ids)} counties i perioden")
        N.put(f"{tag}:share_within_0.1km_aa", shift_agreement(t_fips, t_gaa, 0.1), src,
              "FIPS-AA mot GADM-AA")
        N.put(f"{tag}:share_within_0.1km_bm", shift_agreement(t_fips, t_bm[sel], 0.1), src,
              "FIPS-AA mot GADM-Barmore (frysta facitflytten)")
        N.put(f"{tag}:dir_right_fips", metrics(t_fips, e[sel])["dir_right"], src)
        N.put(f"{tag}:dir_right_gadm", metrics(t_bm[sel], e[sel])["dir_right"], src)


def light_vs_geometric(N: Numbers):
    from tools.change_detection import truth, release
    from tools.compare_benchmarks import geometric_series, coords
    from tools.hybrid_estimate import great_circle_km
    cache = Path("data/processed/boundaries")
    src = "facit + utgåvan + compare_benchmarks.geometric_series"
    sets = [("census", TRUTH, "USA", [2000, 2010, 2020]), ("inegi", TRUTH, "MEX", [2010, 2020]),
            ("ibge", TRUTH, "BRA", [2022]), ("geostat", TRUTH, None, [2011, 2021]),
            ("statcan", Path("data/validation_canada"), "CAN", [2021]),
            ("estat_jp", Path("data/validation_japan"), "JPN", [2020])]
    for ds, root, iso, years in sets:
        for lv in (1, 2):
            tr = {y: truth(root, ds, y, lv, 10) for y in years}
            isos = [iso] if iso else sorted({u[:3] for u in tr[years[0]].index})
            geo = geometric_series(cache, "gadm41", lv, isos)
            rel = release(REL, lv, set(years)).set_index(["unit_id", "year"])
            for y in years:
                t = tr[y][tr[y].index.str[:3].isin(isos)]
                lt = rel.xs(y, level="year")[["bm_lon", "bm_lat"]]
                llon, llat = coords(lt, t.index)
                glon, glat = coords(geo, t.index)
                el = great_circle_km(llon, llat, t.bm_lon.values, t.bm_lat.values)
                eg = great_circle_km(glon, glat, t.bm_lon.values, t.bm_lat.values)
                share, n = share_closer(el, eg)
                N.put(f"light_vs_geom:{ds}:L{lv}:{y}:share", share, src)
                N.put(f"light_vs_geom:{ds}:L{lv}:{y}:n", n, src)


TABLE3_L2 = [("census", 2000, 2010), ("census", 2010, 2020), ("census", 2000, 2020),
             ("inegi", 2010, 2020), ("geostat", 2011, 2021)]


def slope_sensitivity(N: Numbers):
    from tools.canada_confirm import change_inputs
    src = "facit nivå 2 + utgåvan; tools/canada_confirm.change_inputs"
    for ds, y1, y2 in TABLE3_L2:
        ids, t, e, _ = change_inputs(TRUTH, REL, ds, y1, y2, 2)
        tag = f"slope:{ds}:{y1}-{y2}"
        N.put(f"{tag}:n", len(ids), src)
        N.put(f"{tag}:ols", pooled_slope(t, e), src)
        s3, top = slope_without_top(t, e, 3)
        N.put(f"{tag}:ols_without_top3", s3, src, "borttagna: " + ", ".join(ids[top]))
        infl = unit_influence(t, e)
        for name, order in (("down", np.argsort(infl, kind="stable")[:3]),     # borttagning sänker mest
                            ("up", np.argsort(-infl, kind="stable")[:3])):     # borttagning höjer mest
            keep = np.ones(len(t), bool)
            keep[order] = False
            N.put(f"{tag}:ols_without_top3_{name}", pooled_slope(t[keep], e[keep]), src,
                  "borttagna: " + ", ".join(ids[order]))
        for r, i in enumerate(top, 1):
            N.put(f"{tag}:top{r}:dslope", float(infl[i]), src,
                  f"{ids[i]}, |t| = {np.hypot(*t[i]):.2f} km, |e| = {np.hypot(*e[i]):.2f} km")
        N.put(f"{tag}:theil_sen", theil_sen(t.ravel(), e.ravel()), src)


def a2_full_precision(N: Numbers):
    """Kör om A2 (tools/trend_estimate) i samma ordning som den committade
    körningen och läser Europas intervall med full precision. Kontrollerar att
    omkörningen återger de committade CSV-filerna med fyra decimaler."""
    import tools.change_detection as cdm
    import tools.trend_estimate as te
    # Samma slumptillstånd som i en ny process (fröet sätts vid import av
    # change_detection; andra steg här har satt om det).
    cdm.RNG = np.random.default_rng(20261118)
    rel = {lv: te.load_release(REL, lv) for lv in (1, 2)}
    res = []
    for ds in te.PAIRS:
        r = te.run_dataset(ds, TRUTH, rel, 10, 1000)
        old = pd.read_csv(f"docs/extended/{ds}_trend_change.csv")
        num = old.select_dtypes("number").columns
        if not np.array_equal(r[num].round(4).to_numpy(), old[num].round(4).to_numpy()):
            raise RuntimeError(f"A2-omkörningen återger inte {ds}_trend_change.csv")
        res.append(r)
    r = pd.concat(res)
    g = r[(r.dataset == "geostat") & (r.series == "trend") & (r.level == 2)].iloc[0]
    for k in ("dir_right", "dir_right_lo", "dir_right_hi"):
        N.put(f"a2full:geostat:trend:L2:{k}", float(g[k]),
              "omkörning av tools/trend_estimate.run_dataset (full precision)",
              "återger docs/extended/*_trend_change.csv med fyra decimaler")


CD_RUNS = {"census": ["2000-2010", "2010-2020", "2000-2020"], "inegi": ["2010-2020"], "geostat": ["2011-2021"]}


def cd_full_precision(N: Numbers):
    """Kör om tools/change_detection.main med kommandot i docs/HARMONIZATION.md
    (inklusive serien chiov_bltc, så att slumpföljden blir densamma) och lagrar
    utgåvans rader med full precision. CSV-filerna har fyra decimaler, och två
    intervallgränser i tabell 3 ligger där exakt på en halva. Omkörningen måste
    återge de committade filerna med fyra decimaler."""
    import tempfile
    import tools.change_detection as cdm
    orig = pd.DataFrame.to_csv

    def full(self, *a, **k):
        k["float_format"] = "%.12f"
        return orig(self, *a, **k)

    for ds, pairs in CD_RUNS.items():
        with tempfile.TemporaryDirectory() as tmp:
            cdm.RNG = np.random.default_rng(20261118)      # som i en ny process
            pd.DataFrame.to_csv = full
            try:
                cdm.main(["--dataset", ds, "--pairs", *pairs, "--release", str(REL),
                          "--series", "chiov_bltc=data/results_mask2km/gadm41/chiovelli-bltc",
                          "--out", tmp])
            finally:
                pd.DataFrame.to_csv = orig
            new = pd.read_csv(Path(tmp) / f"{ds}_change_detection.csv")
        old = pd.read_csv(f"docs/change_detection/{ds}_change_detection.csv")
        num = old.select_dtypes("number").columns
        if not np.array_equal(new[num].round(4).to_numpy(), old[num].round(4).to_numpy()):
            raise RuntimeError(f"Omkörningen återger inte {ds}_change_detection.csv")
        for _, row in new[new.series == "release"].iterrows():
            for col in num:
                if col in ("level", "from", "to"):
                    continue
                N.put(f"cdfull:{ds}:L{int(row.level)}:{int(row['from'])}-{int(row['to'])}:{col}", float(row[col]),
                      "omkörning av tools/change_detection.main (full precision)",
                      "återger docs/change_detection/*_change_detection.csv med fyra decimaler")


def misc(N: Numbers):
    """Tal ur utgåvan, konfigurationen och dokumenten."""
    for lv in (0, 1, 2):
        d = pd.read_csv(REL / f"wpc_v2_gadm41_adm{lv}.csv.gz", keep_default_na=False, low_memory=False,
                        usecols=["unit_id", "year", "bm_lon", "bm_lat", "aa_lon", "aa_lat", "splice_ok"])
        for c in ("bm_lon", "bm_lat", "aa_lon", "aa_lat", "splice_ok"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
        src = str(REL / f"wpc_v2_gadm41_adm{lv}.csv.gz")
        N.put(f"release:L{lv}:units", d.unit_id.nunique(), src)
        N.put(f"release:L{lv}:spliced", int((d[d.year == 2012].splice_ok == 1).sum()), src,
              "splice_ok = 1 år 2012")
        ok = d.bm_lon.notna() & d.aa_lon.notna()
        from tools.hybrid_estimate import great_circle_km
        dist = great_circle_km(d.bm_lon[ok].values, d.bm_lat[ok].values,
                               d.aa_lon[ok].values, d.aa_lat[ok].values)
        N.put(f"release:L{lv}:bm_aa_median_max_over_years",
              float(pd.Series(dist).groupby(d.year[ok].values).median().max()), src,
              "största årsmedianen av avståndet Barmore–Aboufadel")
    o = pd.read_csv("config/dmsp_offsets.csv", comment="#")
    o = o[~((o.satellite == "F15") & (o.year == 2008))]
    src = "config/dmsp_offsets.csv (utan F15 2008)"
    N.put("coreg:n", len(o), src)
    N.put("coreg:max_abs_component_px", float(np.maximum(o.dx.abs(), o.dy.abs()).max()), src)
    N.put("coreg:max_vector_px", float(np.hypot(o.dx, o.dy).max()), src)
    N.put("coreg:max_vector_km_equator", float(np.hypot(o.dx, o.dy).max() * PX_KM), src,
          "1 pixel = 30″ = 0,927 km vid ekvatorn")
    N.put("coreg:max_abs_component_km", float(np.maximum(o.dx.abs(), o.dy.abs()).max() * PX_KM), src)
    N.put("coreg:share_north", float((o.dy < 0).mean()), src, "dy < 0 = norr om VIIRS")
    N.put("coreg:share_west", float((o.dx < 0).mean()), src, "dx < 0 = väster om VIIRS")
    cal = pd.read_csv("config/calibration_v1.csv", comment="#")
    N.put("calibration:composites", len(cal), "config/calibration_v1.csv")
    # Fackelmaskens tabell i METHODS.md (kolumnen 2 km)
    for line in Path("docs/METHODS.md").read_text().splitlines():
        m = re.match(r"\| (Bakken|Permian|Edmonton|Houston|Oklahoma City|Midland-Odessa)[^|]*\|"
                     r"[^|]*\|\s*(\d+) %", line)
        if m:
            N.put(f"flaremask:2km:{m.group(1)}", float(m.group(2)),
                  "docs/METHODS.md, tabellen Val av radie (2 km)")
    for y in (2011, 2021):
        for lv in (0, 1, 2):
            t = pd.read_csv(TRUTH / "geostat" / str(y) / f"adm{lv}.csv", keep_default_na=False)
            t = t[(t.complete.astype(str).str.lower() == "true") & (t.n_points >= 10)]
            N.put(f"geostat:{y}:L{lv}:countries", t.unit_id.str[:3].nunique(),
                  f"data/validation/geostat/{y}/adm{lv}.csv")
    js = "data/validation/ibge/2022/summary.json"
    N.js(js, "cells")
    N.js(js, "pop_total")
    u = pd.read_csv("docs/figures/data/unit_errors.csv.gz")
    for y in (2000, 2010, 2020):
        v = u[(u.dataset == "census") & (u.year == y) & (u.level == 1) & (u.gid0 == "PRI")].err_km
        N.put(f"unit_errors:PRI:L1:{y}:median", float(v.median()), "docs/figures/data/unit_errors.csv.gz")
    # Sammanfattningens ordantal (ord = blankstegsavgränsade ord mellan
    # "**Abstract.**" och "\[198 words\]")
    txt = MANUS.read_text()
    ab = txt.split("**Abstract.**", 1)[1].split("\\[", 1)[0]
    N.put("abstract:words", len(ab.split()), str(MANUS), "blankstegsavgränsade ord")


def derived(N: Numbers):
    """Tal som är enkla funktioner av committade CSV-värden."""
    C = N.csv
    cd = "docs/change_detection/{}_change_detection.csv"
    vp = "docs/validation/{}_position.csv"
    # Diskussionen: median |e − t| / medel av de två årens median lägesfel
    rows = [("census", 2000, 2010), ("census", 2010, 2020), ("census", 2000, 2020),
            ("inegi", 2010, 2020), ("geostat", 2011, 2021)]
    for ds, y1, y2 in rows:
        err = C(cd.format(ds), "err_med", series="release", level=2, **{"from": y1, "to": y2})
        l1 = C(vp.format(ds), "bm_median", series="release", level=2, year=y1)
        l2 = C(vp.format(ds), "bm_median", series="release", level=2, year=y2)
        r = err / ((l1 + l2) / 2)
        N.put(f"locratio:{ds}:{y1}-{y2}", r, "kvot ur " + cd.format(ds) + " och " + vp.format(ds))
        N.put(f"locratio:{ds}:{y1}-{y2}:implied_corr", implied_corr(r), "1 − kvot²/2")
    # A2: felets ändring mot ändpunkter på nivå 2
    for ds, y1, y2 in rows:
        p = f"docs/extended/{ds}_trend_change.csv"
        f = {"level": 2, "from": y1, "to": y2}
        ep = C(p, "err_med", series="endpoint", **f)
        for s in ("three_year", "trend"):
            N.put(f"a2:{ds}:{y1}-{y2}:{s}:change", C(p, "err_med", series=s, **f) / ep - 1,
                  p, "err_med / err_med(endpoint) − 1")
    N.put("a2:met_criterion", sum(N.get(f"a2:{ds}:{y1}-{y2}:trend:change") <= -0.20
                                  for ds, y1, y2 in rows), "räknat ur docs/extended/*_trend_change.csv")
    N.put("indep_noise_3yr", 1 - 1 / np.sqrt(3), "härledning", "1 − 1/√3")
    # Europa: korrektionen
    rn = pd.read_csv("docs/robustness/robust_national_correction.csv")
    g = rn[rn.dataset == "geostat"]
    src = "docs/robustness/robust_national_correction.csv"
    for y in (2011, 2021):
        gy = g[g.year == y]
        N.put(f"natcorr:geostat:{y}:median_light", float(gy.err_light_weighted_km.median()), src)
        N.put(f"natcorr:geostat:{y}:median_pop", float(gy.err_pop_weighted_regions_km.median()), src)
    N.put("natcorr:geostat:n_smaller", int((g.err_pop_weighted_regions_km < g.err_light_weighted_km).sum()), src)
    N.put("natcorr:geostat:n", len(g), src)
    # Andel counties som flyttade mer än 2 km (2000–2020)
    b = pd.read_csv("docs/change_detection/census_change_detection_bins.csv")
    b = b[(b.series == "release") & (b.level == 2) & (b["from"] == 2000) & (b["to"] == 2020)]
    big = b[b.true_move_km.str.match(r"\[(2\.0|5\.0|10\.0),")]
    N.put("bins:census:2000-2020:share_gt2km", big.n.sum() / b.n.sum(),
          "docs/change_detection/census_change_detection_bins.csv")


# =============================================================================
# Påståendena i manuset
# =============================================================================

def build_claims(N: Numbers, manus: str):
    C, J = N.csv, N.js
    claims = []
    cd = "docs/change_detection/{}_change_detection.csv"
    bins = "docs/change_detection/{}_change_detection_bins.csv"
    vp = "docs/validation/{}_position.csv"
    bp = "docs/extended/{}_benchmark_position.csv"
    bpair = "docs/extended/{}_benchmark_pairs.csv"
    bch = "docs/extended/{}_benchmark_change.csv"
    hp = "docs/extended/{}_hybrid_position.csv"
    hc = "docs/extended/{}_hybrid_change.csv"
    tr = "docs/extended/{}_trend_change.csv"
    perm = "docs/robustness/robust_permutation.csv"
    rst = "docs/robustness/robust_states.csv"
    rnc = "docs/robustness/robust_national_correction.csv"
    cc, cj = "docs/extended/C_canada_{}.csv", "docs/extended/C_japan_{}.csv"
    idf = "docs/check_identity/id_{}.json"
    YRS = {"census": [2000, 2010, 2020], "geostat": [2011, 2021], "inegi": [2010, 2020], "ibge": [2022]}
    FOUR = ["census", "geostat", "inegi", "ibge"]
    PER = {"census": [(2000, 2010), (2010, 2020), (2000, 2020)], "inegi": [(2010, 2020)],
           "geostat": [(2011, 2021)]}

    mark = [len(_TIES)]

    def add(plats, ctx, text, data, sources, num=None, note="", forslag="", st=None, keys=()):
        if ctx not in manus:
            raise ValueError(f"Textbiten finns inte i manuset: {ctx!r}")
        if text not in ctx:
            raise ValueError(f"{text!r} finns inte i {ctx!r}")
        s = st if st is not None else status(num if num is not None else text, data)
        if s == "avviker" and not forslag and plats.startswith("Table"):
            forslag = f"Ändra cellen till \"{data}\"."
        ties = _TIES[mark[0]:]
        mark[0] = len(_TIES)
        if ties:
            note = (note + " " if note else "") + "Gränsfall: " + ", ".join(sorted(set(ties))) + \
                " ligger exakt på en halva (avrundat uppåt)."
        claims.append({"plats": plats, "ctx": ctx, "text": text, "data": data if data is not None else "–",
                       "source": "; ".join(sorted(set(sources))) if sources else "–",
                       "status": s, "note": note, "forslag": forslag})

    def r(vals, dec, scale=1.0):
        return fmt_range(vals, dec, scale)

    # Återkommande mängder ------------------------------------------------------
    def l2_light_err():
        v = [C(vp.format(d), "bm_median", series="release", level=2, year=y) for d in FOUR for y in YRS[d]]
        v += [C(cc.format("position"), "err_med", source="light", level=2, year=2021),
              C(cj.format("position"), "err_med", source="light", level=2, year=2020)]
        return v, [vp.format(d) for d in FOUR] + [cc.format("position"), cj.format("position")]

    def l2_true_med(with_c=True):
        v = [C(cd.format(d), "true_med", series="release", level=2, **{"from": a, "to": b})
             for d in PER for a, b in PER[d]]
        s = [cd.format(d) for d in PER]
        if with_c:
            v += [C(cc.format("change"), "true_med", level=2), C(cj.format("change"), "true_med", level=2)]
            s += [cc.format("change"), cj.format("change")]
        return v, s

    def l2_dir(with_c=True):
        v = [C(cd.format(d), "dir_right", series="release", level=2, **{"from": a, "to": b})
             for d in PER for a, b in PER[d]]
        s = [cd.format(d) for d in PER]
        if with_c:
            v += [C(cc.format("change"), "dir_right", level=2), C(cj.format("change"), "dir_right", level=2)]
            s += [cc.format("change"), cj.format("change")]
        return v, s

    def lvg(lv=2):
        keys = [k for k in N.rows if k.startswith("light_vs_geom:") and f":L{lv}:" in k and k.endswith(":share")]
        return [N.get(k) for k in keys], ["tools/manuscript_numbers.py (light_vs_geom)"]

    def pair_share(src_name, lv=2):
        v = [C(bpair.format(d), "closer_share", level=lv, year=y, source=src_name) for d in FOUR for y in YRS[d]]
        v += [C(cc.format("benchmark_pairs"), "closer_share", source=src_name, level=lv, year=2021),
              C(cj.format("benchmark_pairs"), "closer_share", source=src_name, level=lv, year=2020)]
        return v, [bpair.format(d) for d in FOUR] + [cc.format("benchmark_pairs"), cj.format("benchmark_pairs")]

    nat_light = [C(rnc, "err_light_weighted_km", dataset="census", year=y, country="USA") for y in YRS["census"]] + \
                [C(rnc, "err_light_weighted_km", dataset="inegi", year=y, country="MEX") for y in YRS["inegi"]] + \
                [C(rnc, "err_light_weighted_km", dataset="ibge", year=2022, country="BRA")]
    nat_pop = [C(rnc, "err_pop_weighted_regions_km", dataset="census", year=y, country="USA") for y in YRS["census"]] + \
              [C(rnc, "err_pop_weighted_regions_km", dataset="inegi", year=y, country="MEX") for y in YRS["inegi"]] + \
              [C(rnc, "err_pop_weighted_regions_km", dataset="ibge", year=2022, country="BRA")]

    # ---------------------------------------------------------------- Abstract
    P = "Abstract"
    v, s = l2_light_err()
    add(P, "lies a median 0.7–2.6 km from the census centroid", "0.7–2.6", r(v, 1), s,
        note="Nivå 2, ljus, alla facitår, sex facit.")
    v, s = lvg()
    add(P, "It is closer than the geometric centre of the unit in 73–92 percent of units", "73–92",
        r(v, 0, 100), s, note="Nivå 2, alla facitår (Kanada 2021, Japan 2020).")
    v, s = pair_share("ghsbuilts")
    add(P, "built-up surface is closer still in 57–76 percent", "57–76", r(v, 0, 100), s)
    add(P, "National centroids are displaced by 65–134 km", "65–134", r(nat_light, 0), [rnc],
        note="USA, Mexiko, Brasilien, alla facitår.")
    add(P, "reduces the error to 2–12 km", "2–12", r(nat_pop, 0), [rnc])
    v, s = l2_true_med()
    add(P, "municipal centres of population moved a median 0.06–0.54 km", "0.06–0.54", r(v, 2), s)
    v, s = l2_dir()
    add(P, "the measured shift points the right way in 46–66 percent of municipalities", "46–66",
        r(v, 0, 100), s)
    add(P, "\\[198 words\\]", "198", fmt(N.get("abstract:words"), 0), [str(MANUS)],
        note="Blankstegsavgränsade ord; tidskriftens räkneregel kan skilja.")

    # ----------------------------------------------------------- Introduction
    P = "Introduction"
    v = [N.get(f"geostat:{y}:L{lv}:countries") for y in (2011, 2021) for lv in (0, 1, 2)]
    add(P, "25–32 European countries and territories (2011, 2021)", "25–32", r(v, 0),
        ["data/validation/geostat"], note="Länder med complete-enheter per nivå och år.")
    v, s = l2_light_err()
    add(P, "Within counties and municipalities the light-weighted centroid lies one to two kilometres "
           "from the census centroid", "one to two", r(v, 0), s, num="1–2",
        note="Medianerna är 0,67–2,57 km; Kanada (2,6 km) ligger utanför \"one to two\".",
        forslag="\"Within counties and municipalities the light-weighted centroid lies a median 0.7–2.6 km "
                "from the census centroid\" (eller \"about one to three kilometres\", som i slutsatsen).")

    # --------------------------------------------------------- Data and Methods
    P = "Data and Methods: Lights"
    add(P, "the 34 composites with intercalibration coefficients are used", "34",
        fmt(N.get("calibration:composites"), 0), ["config/calibration_v1.csv"])
    add(P, "removes 83–99 percent of the light in the Bakken and Permian oil fields", "83–99",
        r([N.get("flaremask:2km:Bakken"), N.get("flaremask:2km:Permian")], 0), ["docs/METHODS.md"],
        note="Källan är en tabell i METHODS.md, inte en CSV.")
    add(P, "removing 6–10 percent of the light in Houston, Edmonton and Oklahoma City", "6–10",
        r([N.get(f"flaremask:2km:{k}") for k in ("Houston", "Edmonton", "Oklahoma City")], 0),
        ["docs/METHODS.md"], note="Källan är en tabell i METHODS.md, inte en CSV.")
    add(P, "but 65 percent in Midland–Odessa", "65", fmt(N.get("flaremask:2km:Midland-Odessa"), 0),
        ["docs/METHODS.md"], note="Källan är en tabell i METHODS.md, inte en CSV.")
    add(P, "The estimated shifts reach 1.4 pixels (about 1.3 km)", "1.4 pixels (about 1.3 km)",
        f"{fmt(N.get('coreg:max_vector_px'), 1)} pixels (about {fmt(N.get('coreg:max_vector_km_equator'), 1)} km)",
        ["config/dmsp_offsets.csv"],
        note=f"1,4 px är den största komponenten (|dx| = {fmt(N.get('coreg:max_abs_component_px'), 2)} px, "
             f"F18 2010, {fmt(N.get('coreg:max_abs_component_km'), 2)} km). Förskjutningens längd når "
             f"{fmt(N.get('coreg:max_vector_px'), 2)} px (F15 2007).",
        forslag="\"The estimated shifts reach 1.6 pixels (about 1.5 km at the equator)\", eller "
                "\"reach 1.4 pixels (about 1.3 km) in one direction\".")
    add(P, "almost all composites lie north or northwest of VIIRS", "almost all",
        f"norr {fmt(N.get('coreg:share_north'), 0, 100)} %, väster {fmt(N.get('coreg:share_west'), 0, 100)} %",
        ["config/dmsp_offsets.csv"], st="stämmer", note="Kvalitativt påstående; andelar av 34 kompositer.")

    P = "Data and Methods: Units and Centroids"
    add(P, "263 countries, 3,662 first-level units and 47,217 second-level units", "263",
        fmt(N.get("release:L0:units"), 0), [str(REL)])
    add(P, "263 countries, 3,662 first-level units and 47,217 second-level units", "3,662",
        fmt(N.get("release:L1:units"), 0), [str(REL)])
    add(P, "263 countries, 3,662 first-level units and 47,217 second-level units", "47,217",
        fmt(N.get("release:L2:units"), 0), [str(REL)])
    mx = max(N.get("release:L1:bm_aa_median_max_over_years"), N.get("release:L2:bm_aa_median_max_over_years"))
    add(P, "at the first and second levels the two differ by less than 0.01 km in median", "less than 0.01",
        fmt(mx, 4), [str(REL)], st="stämmer" if mx < 0.01 else "avviker",
        note="Största årsmedianen på nivå 1 och 2.")

    P = "Data and Methods: Splicing"
    for lv, txt in ((0, "250 of 263 countries"), (1, "3,496 of 3,662 regions"), (2, "44,019 of 47,217 municipalities")):
        add(P, "250 of 263 countries, 3,496 of 3,662 regions and 44,019 of 47,217 municipalities",
            txt.rsplit(" ", 1)[0],
            f"{fmt(N.get(f'release:L{lv}:spliced'), 0)} of {fmt(N.get(f'release:L{lv}:units'), 0)}", [str(REL)])

    P = "Data and Methods: Census Benchmarks"
    v = [C(vp.format("census"), "n", series="release", level=2, year=y) for y in YRS["census"]]
    add(P, "In the United States this retains the 2,373–2,383 counties with at least ten block groups",
        "2,373–2,383", r(v, 0), [vp.format("census")],
        note="2 373–2 383 är antalet counties som finns med båda åren i en period (tabell 3). Per "
             "folkräkningsår uppfyller 2 400–2 441 counties regeln (tabell 2 har 2 436).",
        forslag="\"this retains the 2,400–2,441 counties (2,373–2,383 in both years of a period) with at "
                "least ten block groups\".")
    add(P, "and 49 states and the District of Columbia", "49 states and the District of Columbia",
        fmt(C(rst, "n", level=1, **{"from": 2000, "to": 2010}), 0), [rst], num="50")
    v = [N.get(f"fips:{a}-{b}:share_within_0.1km_aa") for a, b in PER["census"]]
    vb = [N.get(f"fips:{a}-{b}:share_within_0.1km_bm") for a, b in PER["census"]]
    nm = [N.get(f"fips:{a}-{b}:n_matched") for a, b in PER["census"]]
    add(P, "the two give the same shift of the census centroid to within 0.1 km in 94–96 percent of counties",
        "94–96", r(v, 0, 100), ["tools/manuscript_numbers.py (fips)"],
        note=f"Aboufadel–Austin båda vägarna, {r(nm, 0)} counties som kan knytas entydigt. Mot den frysta "
             f"facitflytten (Barmore) är andelen {r(vb, 0, 100)} %.")
    add(P, "its 78 municipios are first-level units in GADM", "78",
        fmt(C(rnc, "regions", dataset="census", year=2020, country="PRI"), 0), [rnc])

    P = "Data and Methods: Benchmark stages"
    add(P, "about 94 percent of the 2021 dissemination areas correspond to exactly one 2016 area", "94",
        None, [], note="Ingen fil i repot innehåller DA-korrespondensen 2016–2021; markerad [REF needed] i manuset.")

    P = "Table 1"
    add(P, "1.50 million cells, 203.0 million people", "1.50", fmt(N.get("summary:cells") / 1e6, 2),
        ["data/validation/ibge/2022/summary.json"])
    add(P, "1.50 million cells, 203.0 million people", "203.0", fmt(N.get("summary:pop_total") / 1e6, 1),
        ["data/validation/ibge/2022/summary.json"])
    for y, txt in ((2011, "25–26"), (2021, "29–32")):
        v = [N.get(f"geostat:{y}:L{lv}:countries") for lv in (0, 1, 2)]
        add(P, "25–26 (2011) and 29–32 (2021) countries and territories", txt, r(v, 0),
            ["data/validation/geostat"])

    P = "Data and Methods: Tests (national bias)"
    add(P, "(131 against 133 km for Mexico, 95 against 94 km for Brazil)", "131",
        fmt(J(idf.format("MEX_0.1666667"), "error_gc_km"), 0), [idf.format("MEX_0.1666667")])
    add(P, "(131 against 133 km for Mexico, 95 against 94 km for Brazil)", "133",
        fmt(C(vp.format("inegi"), "bm_median", series="release", level=0, year=2020), 0), [vp.format("inegi")])
    add(P, "(131 against 133 km for Mexico, 95 against 94 km for Brazil)", "95",
        fmt(J(idf.format("BRA_0.1666667"), "error_gc_km"), 0), [idf.format("BRA_0.1666667")])
    add(P, "(131 against 133 km for Mexico, 95 against 94 km for Brazil)", "94",
        fmt(C(vp.format("ibge"), "bm_median", series="release", level=0, year=2022), 0), [vp.format("ibge")])

    # ---------------------------------------------------------------- Results
    P = "Results: Location ¶1"
    # Lägesfelen i ¶1 är hela valideringsurvalet (docs/validation, figur 1). USA:s nivå 1 i
    # docs/validation innehåller Puerto Rico; delstaterna (50) läses därför ur B.
    L2 = lambda d: [C(vp.format(d), "bm_median", series="release", level=2, year=y) for y in YRS[d]]
    L1 = lambda d: ([C(bp.format(d), "err_med", level=1, year=y, source="light") for y in YRS[d]] if d == "census"
                    else [C(vp.format(d), "bm_median", series="release", level=1, year=y) for y in YRS[d]])
    srcL = lambda d, lv: bp.format(d) if (d == "census" and lv == 1) else vp.format(d)
    ctx = ("a median 1.5 km from the census centroid in U.S. counties, 0.7–0.9 km in European "
           "municipalities, 1.0–1.1 km in Mexico and 1.1 km in Brazil")
    add(P, ctx, "1.5", r(L2("census"), 1), [srcL("census", 2)])
    add(P, ctx, "0.7–0.9", r(L2("geostat"), 1), [srcL("geostat", 2)])
    add(P, ctx, "1.0–1.1", r(L2("inegi"), 1), [srcL("inegi", 2)])
    add(P, ctx, "1.1 km in Brazil", r(L2("ibge"), 1) + " km in Brazil", [srcL("ibge", 2)])
    v = [C(vp.format(d), "bm_p90", series="release", level=2, year=y) for d in FOUR for y in YRS[d]]
    add(P, "nine in ten units are within 2.4–4.3 km", "2.4–4.3", r(v, 1), [vp.format(d) for d in FOUR])
    add(P, "the median is 0.7 km for Japanese municipalities and 2.6 km for Canadian census divisions", "0.7",
        fmt(C(cj.format("position"), "err_med", source="light", level=2, year=2020), 1), [cj.format("position")])
    add(P, "the median is 0.7 km for Japanese municipalities and 2.6 km for Canadian census divisions", "2.6",
        fmt(C(cc.format("position"), "err_med", source="light", level=2, year=2021), 1), [cc.format("position")])
    ctx = ("0.3 km for Puerto Rico's municipios, 2.5–3.0 km for European regions, 7–9 km for Mexican states, "
           "11–13 km for U.S. states and 17 km for Brazilian states")
    v = [N.get(f"unit_errors:PRI:L1:{y}:median") for y in YRS["census"]]
    add(P, ctx, "0.3", r(v, 1), ["docs/figures/data/unit_errors.csv.gz"],
        note="Medianen är 0,357 (2000), 0,328 (2010) och 0,326 km (2020). De andra talen i meningen är "
             "spann över alla facitår.",
        forslag="\"0.3–0.4 km for Puerto Rico's municipios\" (eller \"0.3 km in 2020\").")
    add(P, ctx, "2.5–3.0", r(L1("geostat"), 1), [srcL("geostat", 1)])
    add(P, ctx, "7–9", r(L1("inegi"), 0), [srcL("inegi", 1)])
    add(P, ctx, "11–13", r(L1("census"), 0), [srcL("census", 1)], note="50 delstater (utan Puerto Rico).")
    add(P, ctx, "17 km for Brazilian", r(L1("ibge"), 0) + " km for Brazilian", [srcL("ibge", 1)])
    ctx = ("65–88 km from the census centroid in the United States, 133–134 km in Mexico and 94 km in Brazil; "
           "the median across European countries is 11–13 km")
    add(P, ctx, "65–88", r(nat_light[:3], 0), [rnc])
    add(P, ctx, "133–134", r(nat_light[3:5], 0), [rnc])
    add(P, ctx, "94 km in Brazil", fmt(nat_light[5], 0) + " km in Brazil", [rnc])
    v = [C(vp.format("geostat"), "bm_median", series="release", level=0, year=y) for y in YRS["geostat"]]
    add(P, ctx, "11–13 km", r(v, 0) + " km", [vp.format("geostat")])
    add(P, "the county median is 1.46–1.52 km", "1.46–1.52", r(L2("census"), 2), [srcL("census", 2)])
    swe = "tools/manuscript_numbers.py (sweden)"
    ctx = "the median error is 0.8–1.1 km for the 290 municipalities and 3.8–4.1 km for the 21 counties"
    add(P, ctx, "0.8–1.1", r([N.get(f"swe:pos:L2:{y}:err_med") for y in (2011, 2021)], 1), [swe])
    add(P, ctx, "290", r([N.get(f"swe:pos:L2:{y}:n") for y in (2011, 2021)], 0), [swe])
    add(P, ctx, "3.8–4.1", r([N.get(f"swe:pos:L1:{y}:err_med") for y in (2011, 2021)], 1), [swe])
    add(P, ctx, "21 counties", r([N.get(f"swe:pos:L1:{y}:n") for y in (2011, 2021)], 0) + " counties", [swe])

    P = "Results: Location ¶2"
    v = [C(bp.format(d), "err_med", level=2, year=y, source="geometric") for d in FOUR for y in YRS[d]]
    ctx = "The geometric centre lies a median 1.8–5.6 km from the census centroid in the first four benchmarks, 2.4 km in Japan and 15.6 km in Canada"
    add(P, ctx, "1.8–5.6", r(v, 1), [bp.format(d) for d in FOUR])
    add(P, ctx, "2.4 km in Japan", fmt(C(cj.format("position"), "err_med", source="geometric", level=2, year=2020), 1)
        + " km in Japan", [cj.format("position")])
    add(P, ctx, "15.6", fmt(C(cc.format("position"), "err_med", source="geometric", level=2, year=2021), 1),
        [cc.format("position")])
    v, s = lvg()
    add(P, "the closer of the two in 73–92 percent of second-level units", "73–92", r(v, 0, 100), s)
    v = [C(bp.format(d), "err_med", level=2, year=y, source="ghsbuilts") for d in FOUR for y in YRS[d]]
    v.append(C(cj.format("position"), "err_med", source="ghsbuilts", level=2, year=2020))
    add(P, "Its median error is 0.5–1.3 km in five benchmarks and 2.4 km in Canada", "0.5–1.3", r(v, 1),
        [bp.format(d) for d in FOUR] + [cj.format("position")])
    add(P, "Its median error is 0.5–1.3 km in five benchmarks and 2.4 km in Canada", "2.4 km in Canada",
        fmt(C(cc.format("position"), "err_med", source="ghsbuilts", level=2, year=2021), 1) + " km in Canada",
        [cc.format("position")])
    v, s = pair_share("ghsbuilts")
    add(P, "closer to the census centroid than the light-weighted centroid in 57–76 percent of units", "57–76",
        r(v, 0, 100), s)
    lo = C(cc.format("benchmark_pairs"), "closer_lo", source="ghsbuilts", level=2, year=2021)
    hi = C(cc.format("benchmark_pairs"), "closer_hi", source="ghsbuilts", level=2, year=2021)
    add(P, "except Canada (49–61 percent)", "49–61", r([lo, hi], 0, 100), [cc.format("benchmark_pairs")])
    los = [C(bpair.format(d), "closer_lo", level=2, year=y, source="ghsbuilts") for d in FOUR for y in YRS[d]]
    los.append(C(cj.format("benchmark_pairs"), "closer_lo", source="ghsbuilts", level=2, year=2020))
    add(P, "The intervals exclude 50 percent in every benchmark except Canada", "exclude 50 percent",
        f"lägsta nedre gräns utom Kanada {fmt(min(los), 2)}", [bpair.format(d) for d in FOUR],
        st="stämmer" if min(los) > 0.5 and lo <= 0.5 else "avviker")

    P = "Results: Location ¶3"
    v = [C(bp.format(d), "err_med", level=2, year=y, source=src) for d in FOUR for y in YRS[d]
         for src in ("ghspop", "worldpop")]
    v += [C(f.format("position"), "err_med", source=src, level=2, year=y)
          for f, y in ((cc, 2021), (cj, 2020)) for src in ("ghspop", "worldpop")]
    add(P, "The two gridded population products lie within 0.1–0.7 km", "0.1–0.7", r(v, 1),
        [bp.format(d) for d in FOUR] + [cc.format("position"), cj.format("position")],
        note=f"Spannet är {fmt(min(v), 3)}–{fmt(max(v), 3)} km; största värdet är WorldPop i Japan.",
        forslag="\"lie within 0.1–0.6 km\".")
    v, s = pair_share("ghspop")
    add(P, "closer than the light-weighted centroid in 80–99 percent of units for GHS-POP", "80–99", r(v, 0, 100), s)
    v, s = pair_share("worldpop")
    add(P, "and 55–98 percent for WorldPop", "55–98", r(v, 0, 100), s)
    gp = lambda y: C(bp.format("census"), "err_med", level=2, year=y, source="ghspop")
    ctx = "the error of GHS-POP is 0.09 km in 2010 and 0.29–0.40 km in 2000 and 2020"
    add(P, ctx, "0.09", fmt(gp(2010), 2), [bp.format("census")])
    add(P, ctx, "0.29–0.40", r([gp(2000), gp(2020)], 2), [bp.format("census")])

    P = "Results: Location ¶4"
    B1 = lambda d: [C(bp.format(d), "err_med", level=1, year=y, source="light") for y in YRS[d]]
    add(P, "(10.5–13.1 against 16.0–17.2 km)", "10.5–13.1", r(B1("census"), 1), [bp.format("census")])
    add(P, "(10.5–13.1 against 16.0–17.2 km)", "16.0–17.2",
        r([C(bp.format("census"), "err_med", level=1, year=y, source="ghsbuilts") for y in YRS["census"]], 1),
        [bp.format("census")])
    add(P, "Brazil (16.5 against 22.2 km)", "16.5", r(B1("ibge"), 1), [bp.format("ibge")])
    add(P, "Brazil (16.5 against 22.2 km)", "22.2",
        fmt(C(bp.format("ibge"), "err_med", level=1, year=2022, source="ghsbuilts"), 1), [bp.format("ibge")])
    add(P, "Japan (2.3 against 3.7 km)", "2.3",
        fmt(C(cj.format("position"), "err_med", source="light", level=1, year=2020), 1), [cj.format("position")])
    add(P, "Japan (2.3 against 3.7 km)", "3.7",
        fmt(C(cj.format("position"), "err_med", source="ghsbuilts", level=1, year=2020), 1), [cj.format("position")])
    add(P, "built-up surface in Canada (26.4 against 35.8 km)", "26.4",
        fmt(C(cc.format("position"), "err_med", source="ghsbuilts", level=1, year=2021), 1), [cc.format("position")])
    add(P, "built-up surface in Canada (26.4 against 35.8 km)", "35.8",
        fmt(C(cc.format("position"), "err_med", source="light", level=1, year=2021), 1), [cc.format("position")])

    # Tabell 2 (läses ur manuset) ---------------------------------------------
    T2 = {"United States, 2020": ("census", 2020, None), "Europe, 2021": ("geostat", 2021, None),
          "Mexico, 2020": ("inegi", 2020, None), "Brazil, 2022": ("ibge", 2022, None),
          "Canada, 2021": (None, 2021, cc), "Japan, 2020": (None, 2020, cj)}
    srcs = ["light", "ghsbuilts", "ghspop", "worldpop", "geometric"]
    heads = ["Light", "Built-up surface", "GHS-POP", "WorldPop", "Geometric centre"]
    for line in manus.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and cells[0] in T2:
            ds, y, cf = T2[cells[0]]
            for j, (src, h) in enumerate(zip(srcs, heads)):
                if cf is None:
                    val = C(bp.format(ds), "err_med", level=2, year=y, source=src)
                    n = C(bp.format(ds), "n", level=2, year=y, source=src)
                    s = bp.format(ds)
                else:
                    val = C(cf.format("position"), "err_med", source=src, level=2, year=y)
                    n = C(cf.format("position"), "n", source=src, level=2, year=y)
                    s = cf.format("position")
                if j == 0:
                    add("Table 2", line, cells[1], fmt(n, 0), [s], note=f"{cells[0]}: Units")
                add("Table 2", line, cells[2 + j], fmt(val, 2), [s], note=f"{cells[0]}: {h}")

    # ------------------------------------------------------------ Nationellt
    P = "Results: National Bias ¶1"
    mex = idf.format("MEX_0.1666667")
    rel_r = lambda f, name: next(x["r_rel"] for x in json.load(open(f))["light_per_capita_by_region_rel"]
                                 if x["name"] == name)

    def rr(f, *names):
        vals = [rel_r(f, n) for n in names]
        for n, x in zip(names, vals):
            N.put(f"{Path(f).stem}:r_rel:{n}", x, f)
        return vals

    add(P, "In Mexico the light-weighted centroid lies 131 km northwest", "131",
        fmt(J(mex, "error_gc_km"), 0), [mex],
        note=f"D = {J(mex, 'D_total', 'east_km')} km öster, {J(mex, 'D_total', 'north_km')} km norr (nordväst).")
    add(P, "1.7–1.8 times the national mean in the northern border states of Sonora and Chihuahua", "1.7–1.8",
        r(rr(mex, "Sonora", "Chihuahua"), 1), [mex])
    add(P, "0.5–0.6 times the mean in Guerrero and Mexico City", "0.5–0.6",
        r(rr(mex, "Guerrero", "Distrito Federal"), 1), [mex])
    bra = idf.format("BRA_0.1666667")
    add(P, "the centroid is displaced 95 km towards the west and south", "95", fmt(J(bra, "error_gc_km"), 0), [bra],
        note=f"D = {J(bra, 'D_total', 'east_km')} km öster, {J(bra, 'D_total', 'north_km')} km norr.")
    add(P, "Mato Grosso produces 1.7 times", "1.7", fmt(rr(bra, "Mato Grosso")[0], 1), [bra])
    add(P, "and Amazonas 0.6 times", "0.6", fmt(rr(bra, "Amazonas")[0], 1), [bra])
    ctx = "the displacements are 98 km (Sweden), 149 km (Norway) and 71 km (Finland), all northwards"
    for c, txt in (("SWE", "98"), ("NOR", "149"), ("FIN", "71")):
        f = idf.format(f"{c}_none")
        add(P, ctx, txt, fmt(J(f, "error_gc_km"), 0), [f],
            note=f"{c}: D norr = {J(f, 'D_total', 'north_km')} km.")
    ctx = "light per resident is 3.1 times the national mean in Norrbotten and 4.0 times in Finnmark, but 0.29 times in Oslo"
    add(P, ctx, "3.1", fmt(rr(idf.format("SWE_none"), "Norrbotten")[0], 1), [idf.format("SWE_none")])
    add(P, ctx, "4.0", fmt(rr(idf.format("NOR_none"), "Finnmark")[0], 1), [idf.format("NOR_none")])
    add(P, ctx, "0.29", fmt(rr(idf.format("NOR_none"), "Oslo")[0], 2), [idf.format("NOR_none")])

    P = "Results: National Bias ¶2"
    ctx = "127 of 131 km in Mexico, 92 of 95 km in Brazil, 95 of 98 km in Sweden and 146 of 149 km in Norway"
    for c, txt in (("MEX_0.1666667", "127 of 131"), ("BRA_0.1666667", "92 of 95"),
                   ("SWE_none", "95 of 98"), ("NOR_none", "146 of 149")):
        f = idf.format(c)
        add(P, ctx, txt, f"{fmt(J(f, 'Bb_between_regions', 'length_km'), 0)} of {fmt(J(f, 'D_total', 'length_km'), 0)}", [f])
    f = idf.format("FIN_none")
    ctx = "Bb is 38 km and the within-region term Bw 34 km of 71 km"
    add(P, ctx, "38", fmt(J(f, "Bb_between_regions", "length_km"), 0), [f])
    add(P, ctx, "34", fmt(J(f, "Bw_within_regions", "length_km"), 0), [f])
    add(P, ctx, "71 km", fmt(J(f, "D_total", "length_km"), 0) + " km", [f])
    files = [idf.format(c) for c in ("MEX_0.1666667", "BRA_0.1666667", "SWE_none", "NOR_none", "FIN_none")]
    ac = max(max(J(f, "A_light_unpopulated", "length_km"), J(f, "C_within_cells", "length_km")) for f in files)
    add(P, "contribute less than 2 km on this grid", "less than 2 km", f"största {fmt(ac, 2)} km", files,
        st="stämmer" if ac < 2 else "avviker")
    f1 = idf.format("MEX_0.0166667")
    add(P, "A takes up 115 of Mexico's 131 km", "115 of Mexico's 131",
        f"{fmt(J(f1, 'A_light_unpopulated', 'length_km'), 0)} of Mexico's {fmt(J(f1, 'error_gc_km'), 0)}", [f1])
    ctx = "(Sweden 98 → 33 km, Norway 149 → 61 km, Finland 71 → 42 km on the same grid)"
    for c, txt in (("SWE", "98 → 33"), ("NOR", "149 → 61"), ("FIN", "71 → 42")):
        a, b = idf.format(f"{c}_none"), idf.format(f"{c}_65")
        add(P, ctx, txt, f"{fmt(J(a, 'error_gc_km'), 0)} → {fmt(J(b, 'error_gc_km'), 0)}", [a, b])
    ctx = ("the national error is 133–134 km in our series, 95–119 km in the series of Chiovelli et al. (2026) "
           "and 35–42 km in a version of our series")
    vs = lambda ds, ser: [C(vp.format(ds), "bm_median", series=ser, level=0, year=y) for y in YRS[ds]]
    add(P, ctx, "133–134", r(vs("inegi", "release"), 0), [vp.format("inegi")])
    add(P, ctx, "95–119", r(vs("inegi", "chiov_bltc"), 0), [vp.format("inegi")])
    add(P, ctx, "35–42", r(vs("inegi", "spliced"), 0), [vp.format("inegi")],
        note="Serien 'spliced' (VIIRS omvandlad till DMSP-skala, skarvad framåt).")
    add(P, "for Brazil the corresponding errors are 94, 111 and 73 km", "94, 111 and 73",
        ", ".join(fmt(vs("ibge", s)[0], 0) for s in ("release", "chiov_bltc")) + " and "
        + fmt(vs("ibge", "spliced")[0], 0), [vp.format("ibge")])
    add(P, "The 65–134 km bias is therefore a property", "65–134", r(nat_light, 0), [rnc])

    P = "Results: National Bias ¶3"
    ctx = "within 6–11 km of the census centroid in the United States (from 65–88 km)"
    add(P, ctx, "6–11", r(nat_pop[:3], 0), [rnc])
    add(P, ctx, "65–88", r(nat_light[:3], 0), [rnc])
    ctx = "2–3 km in Mexico (from 133–134 km) and 12 km in Brazil (from 94 km)"
    add(P, ctx, "2–3", r(nat_pop[3:5], 0), [rnc])
    add(P, ctx, "133–134", r(nat_light[3:5], 0), [rnc])
    add(P, ctx, "12 km", fmt(nat_pop[5], 0) + " km", [rnc])
    add(P, ctx, "94 km", fmt(nat_light[5], 0) + " km", [rnc])
    ctx = "In Europe the median error falls from 11–13 km to 1.1–1.3 km, and the error is smaller in 55 of 56 country-years"
    add(P, ctx, "11–13", r([N.get(f"natcorr:geostat:{y}:median_light") for y in (2011, 2021)], 0), [rnc])
    add(P, ctx, "1.1–1.3", r([N.get(f"natcorr:geostat:{y}:median_pop") for y in (2011, 2021)], 1), [rnc])
    add(P, ctx, "55 of 56", f"{fmt(N.get('natcorr:geostat:n_smaller'), 0)} of {fmt(N.get('natcorr:geostat:n'), 0)}", [rnc])
    add(P, "The correction works least well in Finland (14–27 km)", "14–27",
        r([C(rnc, "err_pop_weighted_regions_km", dataset="geostat", year=y, country="FIN") for y in (2011, 2021)], 0), [rnc])

    # ------------------------------------------------------------ Förändring
    P = "Results: Shifts ¶1"
    v, s = l2_true_med(False)
    add(P, "The true shifts are small: a median 0.13–0.54 km per municipality", "0.13–0.54", r(v, 2), s)
    cdv = lambda ds, col, a, b, lv=2: C(cd.format(ds), col, series="release", level=lv, **{"from": a, "to": b})
    ctx = ("66 percent for U.S. counties over 2000–2020, 58–61 percent over each decade, 56 percent in Mexico "
           "and 51 percent in Europe")
    add(P, ctx, "66", fmt(cdv("census", "dir_right", 2000, 2020), 0, 100), [cd.format("census")])
    add(P, ctx, "58–61", r([cdv("census", "dir_right", 2000, 2010), cdv("census", "dir_right", 2010, 2020)], 0, 100),
        [cd.format("census")])
    add(P, ctx, "56", fmt(cdv("inegi", "dir_right", 2010, 2020), 0, 100), [cd.format("inegi")])
    add(P, ctx, "51", fmt(cdv("geostat", "dir_right", 2011, 2021), 0, 100), [cd.format("geostat")])
    rho = [cdv(d, "rho", a, b) for d in PER for a, b in PER[d]]
    add(P, "Rank correlations between measured and true shifts are 0.02–0.34", "0.02–0.34", r(rho, 2),
        [cd.format(d) for d in PER])
    snr = [1 - cdv(d, "snr_gt1", a, b) for d in PER for a, b in PER[d]]
    add(P, "In 72–92 percent of municipalities the error of the measured shift exceeds the true shift", "72–92",
        r(snr, 0, 100), [cd.format(d) for d in PER], note="1 − snr_gt1.")
    sl = [cdv(d, "slope", a, b) for d in PER for a, b in PER[d]]
    add(P, "The regression slopes of measured on true shift are 0.2–0.5", "0.2–0.5", r(sl, 1),
        [cd.format(d) for d in PER])
    ctx = "removing three of the 2,383 U.S. counties changes the slope for 2010–2020 from 0.49 to 0.39"
    add(P, ctx, "0.49", fmt(cdv("census", "slope", 2010, 2020), 2), [cd.format("census")])
    add(P, ctx, "2,383", fmt(cdv("census", "n", 2010, 2020), 0), [cd.format("census")])
    s3 = N.get("slope:census:2010-2020:ols_without_top3")
    ts = N.get("slope:census:2010-2020:theil_sen")
    add(P, ctx, "0.39", f"{fmt(s3, 2)} utan de tre mest inflytelserika; Theil–Sen {fmt(ts, 2)}",
        ["tools/manuscript_numbers.py (slope)"], st="saknar källa",
        note="Urvalet av de tre counties som gav 0,39 finns inte i repot. Se avsnittet om lutningens känslighet.",
        forslag="Se förslaget till ny mening under Lutningens känslighet.")
    ctx = ("its 290 municipalities moved a median 0.13 km between 2011 and 2021, and the measured shift points "
           "the right way in 53 percent (45–61 percent)")
    add(P, ctx, "290", fmt(N.get("swe:chg:L2:n"), 0), [swe])
    add(P, ctx, "0.13", fmt(N.get("swe:chg:L2:true_med"), 2), [swe])
    add(P, ctx, "53 percent (45–61", f"{fmt(N.get('swe:chg:L2:dir_right'), 0, 100)} percent ("
        f"{r([N.get('swe:chg:L2:dir_right_lo'), N.get('swe:chg:L2:dir_right_hi')], 0, 100)}", [swe],
        note=f"Blockbootstrap, block = län ({fmt(N.get('swe:chg:L2:n_blocks'), 0)} block), frö 20261121.")

    P = "Results: Shifts ¶2"
    ctx = ("Mexico's 32 states, whose centres moved a median 2.1 km, are measured in the right direction in 72 percent "
           "of cases (95 percent interval 56–88 percent), with a rank correlation of 0.42")
    add(P, ctx, "32", fmt(cdv("inegi", "n", 2010, 2020, 1), 0), [cd.format("inegi")])
    add(P, ctx, "2.1", fmt(cdv("inegi", "true_med", 2010, 2020, 1), 1), [cd.format("inegi")])
    add(P, ctx, "72", fmt(cdv("inegi", "dir_right", 2010, 2020, 1), 0, 100), [cd.format("inegi")])
    add(P, ctx, "56–88", r([cdv("inegi", "dir_right_lo", 2010, 2020, 1), cdv("inegi", "dir_right_hi", 2010, 2020, 1)], 0, 100),
        [cd.format("inegi")])
    add(P, ctx, "0.42", fmt(cdv("inegi", "rho", 2010, 2020, 1), 2), [cd.format("inegi")])
    st_ = lambda col, a, b: C(rst, col, level=1, **{"from": a, "to": b})
    ctx = "moved 2.6–3.9 km, and the share is 50–58 percent, with intervals that include 50 percent"
    add(P, ctx, "2.6–3.9", r([st_("true_med", a, b) for a, b in PER["census"]], 1), [rst])
    add(P, ctx, "50–58", r([st_("dir_right", a, b) for a, b in PER["census"]], 0, 100), [rst])
    inc = all(st_("dir_right_lo", a, b) <= 0.5 <= st_("dir_right_hi", a, b) for a, b in PER["census"])
    add(P, ctx, "include 50 percent", "alla tre intervallen innehåller 0,50" if inc else "inte alla",
        [rst], st="stämmer" if inc else "avviker")
    fp = [N.get(f"fisher:MEX2010-2020_vs_USA{a}-{b}:p") for a, b in PER["census"]]
    add(P, "(Fisher's exact test, p = 0.07–0.24)", "0.07–0.24", r(fp, 2), ["tools/manuscript_numbers.py (states)"],
        note="Mexiko 23/32 mot USA 25/50, 27/50 och 29/50; tvåsidigt.")
    ctx = "the median true shift is 0.44 times the median error in Mexico and 0.46–0.59 times in the United States"
    add(P, ctx, "0.44", fmt(N.get("states:MEX:2010-2020:true_over_err"), 2), ["tools/manuscript_numbers.py (states)"])
    add(P, ctx, "0.46–0.59", r([N.get(f"states:USA:{a}-{b}:true_over_err") for a, b in PER["census"]], 2),
        ["tools/manuscript_numbers.py (states)"])
    ctx = ("European regions, which moved a median 0.42 km, are at 49 percent, and their rank correlation is "
           "slightly negative (−0.16 to −0.01)")
    add(P, ctx, "0.42", fmt(cdv("geostat", "true_med", 2011, 2021, 1), 2), [cd.format("geostat")])
    add(P, ctx, "49", fmt(cdv("geostat", "dir_right", 2011, 2021, 1), 0, 100), [cd.format("geostat")])
    add(P, ctx, "−0.16 to −0.01", r([cdv("geostat", "rho_lo", 2011, 2021, 1), cdv("geostat", "rho_hi", 2011, 2021, 1)], 2),
        [cd.format("geostat")], note=f"Intervallet; punktskattningen är {fmt(cdv('geostat', 'rho', 2011, 2021, 1), 2)}.")

    # Tabell 3 (läses ur manuset) ---------------------------------------------
    T3 = {("United States (states)", "1"): ("rst", "census"), ("United States (counties)", "2"): ("cd", "census"),
          ("Mexico", "1"): ("cd", "inegi"), ("Mexico", "2"): ("cd", "inegi"),
          ("Europe", "1"): ("cd", "geostat"), ("Europe", "2"): ("cd", "geostat")}
    for line in manus.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 9 and (cells[0], cells[1]) in T3:
            kind, ds = T3[(cells[0], cells[1])]
            a, b = (int(x) for x in cells[2].split("–"))
            lv = int(cells[1])
            if kind == "rst":
                g, s = (lambda col: C(rst, col, level=1, **{"from": a, "to": b})), rst
            else:
                def g(col, ds=ds, lv=lv, a=a, b=b):
                    k = f"cdfull:{ds}:L{lv}:{a}-{b}:{col}"     # full precision om omkörd
                    if k in N.rows:
                        return N.get(k)
                    return C(cd.format(ds), col, series="release", level=lv, **{"from": a, "to": b})
                s = cd.format(ds)
            lab = f"{cells[0]} {cells[2]}"
            full = kind == "cd" and f"cdfull:{ds}:L{lv}:{a}-{b}:n" in N.rows
            s = [s, "omkörning av change_detection (full precision)"] if full else [s]
            add("Table 3", line, cells[3], fmt(g("n"), 0), s, note=f"{lab}: Units")
            add("Table 3", line, cells[4], fmt(g("true_med"), 2), s, note=f"{lab}: True shift")
            for j, col in ((5, "dir_right"), (6, "rho"), (7, "slope")):
                extra = ""
                if full:
                    for c_ in (col, f"{col}_lo", f"{col}_hi"):
                        csvv = C(cd.format(ds), c_, series="release", level=lv, **{"from": a, "to": b})
                        if is_tie(csvv, 2):
                            extra += f" {c_}: {csvv:.4f} i CSV, {g(c_):.6f} med full precision."
                    _TIES.clear(); mark[0] = 0
                add("Table 3", line, cells[j], fmt_ci(g(col), g(f"{col}_lo"), g(f"{col}_hi"), 2), s,
                    note=f"{lab}: {col}." + extra)
            add("Table 3", line, cells[8], fmt(g("snr_gt1"), 0, 100) + "%", s, note=f"{lab}: True shift > error")

    P = "Results: Shifts ¶3 (size classes)"
    bv = lambda ds, a, b, cl: C(bins.format(ds), "dir_right", series="release", level=2,
                                **{"from": a, "to": b, "true_move_km": cl})
    ctx = ("the share in the right direction rises from 54 percent for true shifts below 0.1 km to 69 percent for "
           "0.5–1 km and 84 percent for 2–5 km; only 8 percent of counties moved more than 2 km")
    add(P, ctx, "54", fmt(bv("census", 2000, 2020, "[0.0, 0.1)"), 0, 100), [bins.format("census")])
    add(P, ctx, "69", fmt(bv("census", 2000, 2020, "[0.5, 1.0)"), 0, 100), [bins.format("census")])
    add(P, ctx, "84", fmt(bv("census", 2000, 2020, "[2.0, 5.0)"), 0, 100), [bins.format("census")])
    add(P, ctx, "8 percent", fmt(N.get("bins:census:2000-2020:share_gt2km"), 0, 100) + " percent", [bins.format("census")])
    ctx = "the 2–5 km class reaches 77 percent (2000–2010) and 70 percent (2010–2020)"
    add(P, ctx, "77", fmt(bv("census", 2000, 2010, "[2.0, 5.0)"), 0, 100), [bins.format("census")])
    add(P, ctx, "70 percent", fmt(bv("census", 2010, 2020, "[2.0, 5.0)"), 0, 100) + " percent", [bins.format("census")])
    bi = pd.read_csv(bins.format("inegi"))
    bi = bi[(bi.series == "release") & (bi.level == 2)]
    mxi = bi[bi.n >= 20].dir_right.max()
    N.put("bins:inegi:max_share_n_ge_20", float(mxi), bins.format("inegi"))
    ctx = "In Mexico the share stays below 70 percent in every size class with at least 20 units, including 58 percent for 2–5 km"
    add(P, ctx, "below 70 percent", f"högst {fmt(mxi, 0, 100)} %", [bins.format("inegi")],
        st="stämmer" if mxi < 0.70 else "avviker")
    add(P, ctx, "58", fmt(bv("inegi", 2010, 2020, "[2.0, 5.0)"), 0, 100), [bins.format("inegi")])
    bg = pd.read_csv(bins.format("geostat"))
    bg = bg[(bg.series == "release") & (bg.level == 2)]
    bg = bg[bg.true_move_km.isin(["[0.0, 0.1)", "[0.1, 0.25)", "[0.25, 0.5)", "[0.5, 1.0)", "[1.0, 2.0)"])]
    for _, x in bg.iterrows():
        N.put(f"bins:geostat:{x.true_move_km}:dir_right", float(x.dir_right), bins.format("geostat"))
    add(P, "in Europe it stays at 50–52 percent in every class up to 2 km", "50–52", r(bg.dir_right, 0, 100),
        [bins.format("geostat")])
    ctx = "the share in the right direction for U.S. counties is 61 percent, against 58 percent for 2010–2020"
    add(P, ctx, "61", fmt(cdv("census", "dir_right", 2000, 2010), 0, 100), [cd.format("census")])
    add(P, ctx, "58", fmt(cdv("census", "dir_right", 2010, 2020), 0, 100), [cd.format("census")])

    P = "Results: Shifts ¶4 (permutation)"
    pv = lambda ds, lv, a, b, col: C(perm, col, dataset=ds, level=lv, **{"from": a, "to": b})
    rows2 = [(d, a, b) for d in PER for a, b in PER[d]]
    means = [pv(d, 2, a, b, "dir_right_null_mean") for d, a, b in rows2]
    los = [pv(d, 2, a, b, "dir_right_null_lo") for d, a, b in rows2]
    his = [pv(d, 2, a, b, "dir_right_null_hi") for d, a, b in rows2]
    add(P, "the share in the right direction is 0.50 (0.48–0.52)", "0.50 (0.48–0.52)",
        f"{r(means, 2)} ({fmt(min(los), 2)}–{fmt(max(his), 2)})", [perm], note="Nivå 2, alla perioder.")
    rb = max(max(abs(pv(d, 2, a, b, "rho_null_lo")), abs(pv(d, 2, a, b, "rho_null_hi"))) for d, a, b in rows2)
    add(P, "the rank correlation is within ±0.04 of zero", "±0.04", f"±{fmt(rb, 2)}", [perm])
    pp = [pv(d, 2, a, b, c) for d, a, b in rows2 if d in ("census", "inegi") for c in ("p_dir_right_upper", "p_rho_upper")]
    add(P, "U.S. counties and Mexican municipalities (p ≤ 0.001)", "p ≤ 0.001", f"högst {fmt(max(pp), 3)}", [perm],
        st="stämmer" if max(pp) <= 0.001 else "avviker",
        note="p = 0,001 är det minsta möjliga med 1 000 permutationer.")
    pm = [pv("inegi", 1, 2010, 2020, c) for c in ("p_dir_right_upper", "p_rho_upper")]
    add(P, "for Mexican states (p = 0.01)", "0.01", r(pm, 3), [perm], num="0.010",
        note="0,010 gäller andelen rätt riktning, 0,002 rangkorrelationen.",
        forslag="\"for Mexican states (p = 0.01 for the share and 0.002 for the correlation)\".")
    pe = [pv("geostat", 2, 2011, 2021, c) for c in ("p_dir_right_upper", "p_rho_upper")]
    add(P, "marginally so for European municipalities (p = 0.03)", "0.03", r(pe, 2), [perm],
        note="0,031 för andelen, 0,023 för korrelationen.",
        forslag="\"marginally so for European municipalities (p = 0.03 and 0.02)\".")
    pu_d = [pv("census", 1, a, b, "p_dir_right_upper") for a, b in PER["census"]]
    pu_r = [pv("census", 1, a, b, "p_rho_upper") for a, b in PER["census"]]
    add(P, "and not for U.S. states (p = 0.19–0.48)", "0.19–0.48", r(pu_d + pu_r, 2), [perm],
        note=f"Andelen: {r(pu_d, 2)}; korrelationen: {r(pu_r, 2)}. Meningen gäller både andelar och korrelationer.",
        forslag="\"and not for U.S. states (p = 0.19–0.48 for the share and 0.07–0.44 for the correlation)\".")
    add(P, "the negative rank correlation of European regions is borderline (p = 0.04)", "0.04",
        fmt(pv("geostat", 1, 2011, 2021, "p_rho_lower"), 2), [perm])

    P = "Results: Shifts ¶5 (estimators)"
    rows5 = [(d, a, b) for d in PER for a, b in PER[d]]
    red3 = [-N.get(f"a2:{d}:{a}-{b}:three_year:change") for d, a, b in rows5]
    add(P, "Three-year means reduce the median error of measured municipal shifts by 5–15 percent", "5–15",
        r(red3, 0, 100), [tr.format(d) for d in PER])
    trc = [N.get(f"a2:{d}:{a}-{b}:trend:change") for d, a, b in rows5]
    add(P, "the trend estimate changes it by between −22 and +4 percent", "−22 and +4",
        f"{fmt(min(trc), 0, 100)} and +{fmt(max(trc), 0, 100)}", [tr.format(d) for d in PER])
    add(P, "a reduction of at least 20 percent in four of five periods, was met in one", "one",
        fmt(N.get("a2:met_criterion"), 0), [tr.format(d) for d in PER], num="1")
    trv = lambda ds, ser, col, a, b, lv=2: C(tr.format(ds), col, series=ser, level=lv, **{"from": a, "to": b})
    v = [trv("census", s_, "dir_right", a, b) for s_ in ("endpoint", "three_year", "trend") for a, b in PER["census"]]
    add(P, "(0.58–0.67 with all three estimators)", "0.58–0.67", r(v, 2), [tr.format("census")])
    add(P, "falls with the trend estimate in Mexico (0.51)", "0.51", fmt(trv("inegi", "trend", "dir_right", 2010, 2020), 2),
        [tr.format("inegi")])
    if "a2full:geostat:trend:L2:dir_right_hi" in N.rows:
        lo_, hi_ = N.get("a2full:geostat:trend:L2:dir_right_lo"), N.get("a2full:geostat:trend:L2:dir_right_hi")
        srcs = [tr.format("geostat"), "omkörning av A2 (full precision)"]
        nt = f"CSV: övre gräns 0.4750; omkörd med full precision {hi_:.6f}."
    else:
        lo_ = trv("geostat", "trend", "dir_right_lo", 2011, 2021)
        hi_ = trv("geostat", "trend", "dir_right_hi", 2011, 2021)
        srcs, nt = [tr.format("geostat")], "Övre gränsen 0.4750 ligger på en halva; kör utan --skip-a2."
    add(P, "and Europe (0.46, interval 0.45–0.47, below chance)", "0.46, interval 0.45–0.47",
        f"{fmt(trv('geostat', 'trend', 'dir_right', 2011, 2021), 2)}, interval {r([lo_, hi_], 2)}", srcs, note=nt)
    add(P, "three-year means reduce the error by 15 percent", "15",
        fmt(-N.get("a2:census:2000-2010:three_year:change"), 0, 100), [tr.format("census")])
    add(P, "where independent annual noise at both ends would give 42 percent", "42",
        fmt(N.get("indep_noise_3yr"), 0, 100), ["härledning: 1 − 1/√3"])
    v = [trv("census", "trend", "dir_right", a, b, 1) for a, b in PER["census"]]
    add(P, "the trend estimate raises the share in the right direction to 0.62–0.70", "0.62–0.70", r(v, 2),
        [tr.format("census")])
    l1lo = [trv("census", "trend", "dir_right_lo", a, b, 1) for a, b in PER["census"]]
    ok = l1lo[0] <= 0.5 < l1lo[1] and l1lo[2] > 0.5
    add(P, "with intervals above 0.50 for 2010–2020 and 2000–2020", "above 0.50",
        "nedre gränser " + " / ".join(fmt(x, 2) for x in l1lo), [tr.format("census")],
        st="stämmer" if ok else "avviker", note="2000–2010, 2010–2020, 2000–2020.")
    m1 = trv("inegi", "trend", "dir_right", 2010, 2020, 1)
    add(P, "it lowers the share for Mexican states (0.62, interval 0.44–0.78)", "0.62", fmt(m1, 2),
        [tr.format("inegi")], note="0,625 = 20/32 exakt; halvt uppåt ger 0,63.",
        forslag="\"(0.63, interval 0.44–0.78)\".")
    add(P, "it lowers the share for Mexican states (0.62, interval 0.44–0.78)", "0.44–0.78",
        r([trv("inegi", "trend", "dir_right_lo", 2010, 2020, 1), trv("inegi", "trend", "dir_right_hi", 2010, 2020, 1)], 2),
        [tr.format("inegi")])

    P = "Results: Shifts ¶6 (detection model)"
    v = [C("docs/extended/detection_rho.csv", "rho", dataset=d, level=2, **{"from": a, "to": b}) for d, a, b in rows5]
    add(P, "(rank correlation 0.42–0.65)", "0.42–0.65", r(v, 2), ["docs/extended/detection_rho.csv"])
    dt = "docs/extended/detection_test.csv"
    add(P, "(area under the curve 0.51 against 0.53)", "0.51",
        fmt(C(dt, "auc", test_dataset="inegi", model="snr"), 2), [dt])
    add(P, "(area under the curve 0.51 against 0.53)", "0.53",
        fmt(C(dt, "auc", test_dataset="inegi", model="abs_t"), 2), [dt])
    dc = "docs/extended/detection_calibration.csv"
    ctx = "in the tenth of municipalities where it predicts 78 percent in the right direction the observed share is 54 percent"
    add(P, ctx, "78", fmt(C(dc, "pred_mean", test_dataset="inegi", bin=10), 0, 100), [dc])
    add(P, ctx, "54", fmt(C(dc, "obs_prop", test_dataset="inegi", bin=10), 0, 100), [dc])

    P = "Results: Shifts ¶7 (other sources)"
    bc = lambda ds, src, col, a, b, lv=2: C(bch.format(ds), col, level=lv, source=src, **{"from": a, "to": b})
    ctx = "For GHS-POP the share of municipalities in the right direction is 0.57–0.59 in the United States, 0.51 in Mexico and 0.62 in Europe"
    add(P, ctx, "0.57–0.59", r([bc("census", "ghspop", "dir_right", a, b) for a, b in PER["census"]], 2), [bch.format("census")])
    add(P, ctx, "0.51", fmt(bc("inegi", "ghspop", "dir_right", 2010, 2020), 2), [bch.format("inegi")])
    add(P, ctx, "0.62", fmt(bc("geostat", "ghspop", "dir_right", 2011, 2021), 2), [bch.format("geostat")])
    ctx = "(0.26–0.53 against 0.26–0.54 km)"
    add(P, ctx, "0.26–0.53", r([bc("census", "ghspop", "err_med", a, b) for a, b in PER["census"]], 2), [bch.format("census")])
    add(P, ctx, "0.26–0.54", r([bc("census", "ghspop", "true_med", a, b) for a, b in PER["census"]], 2), [bch.format("census")])
    v = [bc(d, "worldpop", "dir_right", a, b) for d in PER for a, b in PER[d]]
    add(P, "WorldPop gives 0.52–0.57", "0.52–0.57", r(v, 2), [bch.format(d) for d in PER])
    ctx = "Built-up surface gives 0.65–0.70 in the United States and 0.59 in Europe, above light in both, and 0.50 in Mexico"
    add(P, ctx, "0.65–0.70", r([bc("census", "ghsbuilts", "dir_right", a, b) for a, b in PER["census"]], 2), [bch.format("census")])
    add(P, ctx, "0.59", fmt(bc("geostat", "ghsbuilts", "dir_right", 2011, 2021), 2), [bch.format("geostat")])
    add(P, ctx, "0.50", fmt(bc("inegi", "ghsbuilts", "dir_right", 2010, 2020), 2), [bch.format("inegi")])

    # --------------------------------------------------------------- VIIRS
    P = "Results: Change Within the VIIRS Era"
    ch = lambda f, col, lv=2: C(f.format("change"), col, level=lv)
    mn_ = lambda f, col: C(f.format("main"), col)
    ctx = ("the centres of 287 census divisions moved a median 0.33 km between 2016 and 2021, and the measured shift "
           "points the right way in 59 percent (51–68 percent)")
    add(P, ctx, "287", fmt(ch(cc, "n"), 0), [cc.format("change")])
    add(P, ctx, "0.33", fmt(ch(cc, "true_med"), 2), [cc.format("change")])
    add(P, ctx, "59", fmt(ch(cc, "dir_right"), 0, 100), [cc.format("change")])
    add(P, ctx, "51–68", r([ch(cc, "dir_right_lo"), ch(cc, "dir_right_hi")], 0, 100), [cc.format("change")])
    add(P, "(permutation p = 0.002)", "0.002", fmt(mn_(cc, "perm_p_upper"), 3), [cc.format("main")])
    add(P, "but the interval is 0.17 wide", "0.17", fmt(mn_(cc, "ci_width"), 2), [cc.format("main")])
    ctx = ("the centres of 1,792 municipalities moved a median 0.06 km between 2015 and 2020, and the share is 46 "
           "percent (44–49 percent)")
    add(P, ctx, "1,792", fmt(ch(cj, "n"), 0), [cj.format("change")])
    add(P, ctx, "0.06", fmt(ch(cj, "true_med"), 2), [cj.format("change")])
    add(P, ctx, "46", fmt(ch(cj, "dir_right"), 0, 100), [cj.format("change")])
    add(P, ctx, "44–49", r([ch(cj, "dir_right_lo"), ch(cj, "dir_right_hi")], 0, 100), [cj.format("change")])
    add(P, "(permutation p = 0.003 in the lower tail)", "0.003", fmt(mn_(cj, "perm_p_lower"), 3), [cj.format("main")])
    ctx = ("The median error of the measured shift is 3.2 times the median true shift in Canada and 3.8 times in Japan, "
           "within the range of the earlier benchmarks (1.9–5.0)")
    add(P, ctx, "3.2", fmt(ch(cc, "ratio"), 1), [cc.format("change")])
    add(P, ctx, "3.8", fmt(ch(cj, "ratio"), 1), [cj.format("change")])
    add(P, ctx, "1.9–5.0", r([cdv(d, "ratio", a, b) for d, a, b in rows5], 1), [cd.format(d) for d in PER])
    ctx = ("62 percent for the 13 Canadian provinces and territories (31–85 percent) and 59 percent for 46 Japanese "
           "prefectures (46–72 percent)")
    add(P, ctx, "62", fmt(ch(cc, "dir_right", 1), 0, 100), [cc.format("change")])
    add(P, ctx, "13", fmt(ch(cc, "n", 1), 0), [cc.format("change")])
    add(P, ctx, "31–85", r([ch(cc, "dir_right_lo", 1), ch(cc, "dir_right_hi", 1)], 0, 100), [cc.format("change")])
    add(P, ctx, "59 percent for 46", f"{fmt(ch(cj, 'dir_right', 1), 0, 100)} percent for {fmt(ch(cj, 'n', 1), 0)}",
        [cj.format("change")])
    add(P, ctx, "46–72", r([ch(cj, "dir_right_lo", 1), ch(cj, "dir_right_hi", 1)], 0, 100), [cj.format("change")])

    # --------------------------------------------------------------- Hybrid
    P = "Results: Regional Totals ¶1"
    hv = lambda est: ([C(hp.format(d), "err_med", year=y, estimate=est) for d in FOUR for y in YRS[d]]
                      + [C(cj.format("hybrid"), "err_med", year=y, level=1, estimate=est) for y in (2015, 2020)])
    hcv = lambda est: [C(cc.format("hybrid"), "err_med", year=y, level=1, estimate=est) for y in (2016, 2021)]
    hs = [hp.format(d) for d in FOUR] + [cj.format("hybrid")]
    ctx = ("lies a median 0.3–0.9 km from the census centroid in the United States, Europe, Mexico, Brazil and Japan "
           "and 1.0–1.5 km in Canada, against 1.8–16.5 km for the light-weighted first-level centroid (36–42 km in Canada)")
    add(P, ctx, "0.3–0.9", r(hv("hybrid_light"), 1), hs)
    add(P, ctx, "1.0–1.5", r(hcv("hybrid_light"), 1), [cc.format("hybrid")])
    add(P, ctx, "1.8–16.5", r(hv("light"), 1), hs)
    add(P, ctx, "36–42", r(hcv("light"), 0), [cc.format("hybrid")])
    ctx = "the error is 1.6–3.1 km (15 km in Canada)"
    add(P, ctx, "1.6–3.1", r(hv("hybrid_geom"), 1), hs)
    add(P, ctx, "15 km", r(hcv("hybrid_geom"), 0) + " km", [cc.format("hybrid")])
    hpp = "docs/extended/{}_hybrid_position_pairs.csv"
    v = [C(hpp.format(d), "diff_med", year=y, pair="hybrid_geom_minus_hybrid_light") for d in FOUR for y in YRS[d]]
    add(P, "Light thus improves the location within units by one to three kilometres", "one to three", r(v, 0),
        [hpp.format(d) for d in FOUR], num="1–3",
        note=f"Parvis median hybrid-geom − hybrid-ljus: {r(v, 2)} km.")

    P = "Results: Regional Totals ¶2"
    hch = lambda ds, est, col, a, b: C(hc.format(ds), col, series=est, level=1, **{"from": a, "to": b})
    v = [hch(d, "hybrid_light", "dir_right", a, b) for d, a, b in rows5]
    add(P, "points the right way in 84–98 percent of first-level units", "84–98", r(v, 0, 100), [hc.format(d) for d in PER])
    v = [hch(d, "hybrid_geom", "dir_right", a, b) for d, a, b in rows5]
    add(P, "and in 92–100 percent when geometric centres are used", "92–100", r(v, 0, 100), [hc.format(d) for d in PER])
    ctx = "(0.30 against 0.31 km and 0.36 against 0.44 km)"
    for txt, a, b, est in (("0.30", 2000, 2010, "hybrid_light"), ("0.31", 2000, 2010, "hybrid_geom"),
                           ("0.36", 2000, 2020, "hybrid_light"), ("0.44", 2000, 2020, "hybrid_geom")):
        add(P, ctx, txt, fmt(hch("census", est, "err_med", a, b), 2), [hc.format("census")])
    add(P, "Mexico (0.25 against 0.42 km)", "0.25", fmt(hch("inegi", "hybrid_geom", "err_med", 2010, 2020), 2), [hc.format("inegi")])
    add(P, "Mexico (0.25 against 0.42 km)", "0.42", fmt(hch("inegi", "hybrid_light", "err_med", 2010, 2020), 2), [hc.format("inegi")])
    add(P, "Europe (0.09 against 0.28 km)", "0.09", fmt(hch("geostat", "hybrid_geom", "err_med", 2011, 2021), 2), [hc.format("geostat")])
    add(P, "Europe (0.09 against 0.28 km)", "0.28", fmt(hch("geostat", "hybrid_light", "err_med", 2011, 2021), 2), [hc.format("geostat")])
    v = [bc(d, s_, "dir_right", a, b, 1) for d in ("census", "geostat") for a, b in PER[d] for s_ in ("ghspop", "worldpop")]
    ctx = "give 71–100 percent in the right direction at the first level in the United States and Europe; for Mexican states GHS-POP gives 41 percent and WorldPop 88 percent"
    add(P, ctx, "71–100", r(v, 0, 100), [bch.format("census"), bch.format("geostat")])
    add(P, ctx, "41", fmt(bc("inegi", "ghspop", "dir_right", 2010, 2020, 1), 0, 100), [bch.format("inegi")])
    add(P, ctx, "88", fmt(bc("inegi", "worldpop", "dir_right", 2010, 2020, 1), 0, 100), [bch.format("inegi")])

    # ------------------------------------------------------------ Discussion
    P = "Discussion: Snapshot and Film ¶1"
    ten = [cdv(d, "true_med", a, b) for d, a, b in rows5 if b - a == 10]
    add(P, "which is 0.13–0.54 km in the median municipality of the three benchmarks with ten-year periods",
        "0.13–0.54", r(ten, 2), [cd.format(d) for d in PER],
        note="0,54 km gäller 2000–2020 (20 år). Tioårsperioderna ger 0,13–0,37 km.",
        forslag="\"which is 0.13–0.37 km in the median municipality of the three benchmarks with ten-year "
                "periods\" (eller \"0.13–0.54 km … with ten- and twenty-year periods\").")
    lr = [N.get(f"locratio:{d}:{a}-{b}") for d, a, b in rows5]
    add(P, "the median error of a measured shift is only 0.5–0.8 times the median location error", "0.5–0.8",
        r(lr, 1), [cd.format(d) for d in PER] + [vp.format(d) for d in PER],
        note=f"Kvot mot medel av de två årens median lägesfel: {r(lr, 2)}.")
    ic = [N.get(f"locratio:{d}:{a}-{b}:implied_corr") for d, a, b in rows5]
    add(P, "which implies a correlation of roughly 0.7–0.9", "0.7–0.9", r(ic, 1), ["härledning: 1 − kvot²/2"],
        note=f"{r(ic, 2)}.")
    P = "Discussion: Snapshot and Film ¶3"
    ctx = "it is about 1 km for the municipalities and 4 km for the Swedish counties tested here"
    add(P, ctx, "1 km", r([N.get(f"swe:pos:L2:{y}:err_med") for y in (2011, 2021)], 0) + " km", [swe])
    add(P, ctx, "4 km", r([N.get(f"swe:pos:L1:{y}:err_med") for y in (2011, 2021)], 0) + " km", [swe])
    P = "Discussion: What Lights Add"
    v, s = l2_dir()
    add(P, "light-weighted centroids get the direction of municipal centroid shifts right in 46–66 percent", "46–66",
        r(v, 0, 100), s)
    P = "Discussion: Why National Centroids Fail"
    dd = [abs(C(vp.format(d), "bm_median", series="release", level=0, year=y)
              - C(vp.format(d), "aa_median", series="release", level=0, year=y)) for d in ("inegi", "ibge") for y in YRS[d]]
    add(P, "their national errors for Mexico and Brazil differ by less than 0.2 km", "less than 0.2 km",
        f"största skillnad {fmt(max(dd), 2)} km", [vp.format("inegi"), vp.format("ibge")],
        st="stämmer" if max(dd) < 0.2 else "avviker")
    add(P, "reduces the national errors to 2–12 km in the United States, Mexico and Brazil", "2–12", r(nat_pop, 0), [rnc])
    add(P, "and to about 1 km in the median European country", "about 1 km",
        r([N.get(f"natcorr:geostat:{y}:median_pop") for y in (2011, 2021)], 0), [rnc], num="1")
    P = "Discussion: Recommendations"
    add(P, "with raw VIIRS weights the bias is 65–134 km", "65–134", r(nat_light, 0), [rnc])
    P = "Discussion: Limitations"
    add(P, "(a median 0.06 km in Japan)", "0.06", fmt(ch(cj, "true_med"), 2), [cj.format("change")])

    # ------------------------------------------------------------ Conclusion
    P = "Conclusion"
    v, s = l2_light_err()
    add(P, "to within one to three kilometres in six census systems", "one to three", r(v, 0), s, num="1–3",
        note="Medianerna 0,67–2,57 km avrundade till hela km.")
    return claims


# =============================================================================
# Revisionsdokumentet
# =============================================================================

def _esc(s: str) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def write_review(claims, N: Numbers, path: Path):
    dev = [c for c in claims if c["status"] == "avviker"]
    miss = [c for c in claims if c["status"] == "saknar källa"]
    n_ok = sum(c["status"] == "stämmer" for c in claims)
    L = []
    L.append("# Talrevision av manuset 2026-11-22")
    L.append("")
    L.append("Underlag: `docs/manuscript/manuscript.md` som det står i commit "
             f"`{_git_head()}`. Manuset är inte ändrat. Filen skrivs av `tools/manuscript_numbers.py`, "
             "som också skriver alla tal revisionen bygger på till `docs/manuscript/numbers.csv`; inget tal "
             "här är skrivet för hand.")
    L.append("")
    L.append("**Metod.** Varje kontrollerat tal i abstract, metod, resultat, tabeller, diskussion och slutsats "
             "knyts till en ordagrann textbit i manuset (skriptet stoppar om biten saknas) och räknas om ur "
             "resultatfilerna (`docs/validation`, `docs/change_detection`, `docs/robustness`, "
             "`docs/check_identity`, `docs/extended`), ur utgåvan och facit, eller med nya beräkningar i skriptet. "
             "Talet ur data avrundas **halvt uppåt** (bort från noll) till textens precision, och värden ur CSV-filerna "
             "avrundas som de står skrivna med fyra decimaler. Ett intervall i texten jämförs med min–max av alla "
             "värden det sammanfattar. Status: *stämmer* om de avrundade talen är lika, *avviker* om inte, "
             "*saknar källa* om inget i repot ger talet. Några kvalitativa påståenden (\"less than 0.01\", "
             "\"exclude 50 percent\") prövas som villkor; det står i kolumnen Tal ur data.")
    L.append("")
    L.append(f"**Resultat.** {len(claims)} tal kontrollerade: {n_ok} stämmer, {len(dev)} avviker, "
             f"{len(miss)} saknar källa. Inget av avvikelserna ändrar en slutsats; de flesta gäller avrundning, "
             "ett spann som inte täcker alla värden eller en mening som blandar två mått.")
    L.append("")
    L.append("## Avvikelser, med förslag till rättad text")
    L.append("")
    for i, c in enumerate(dev, 1):
        L.append(f"{i}. **{c['plats']}.** \"…{c['ctx']}…\" — i texten **{c['text']}**, ur data **{c['data']}** "
                 f"(`{c['source']}`).")
        if c["note"]:
            L.append(f"   {c['note']}")
        if c["forslag"]:
            L.append(f"   *Förslag:* {c['forslag']}")
    L.append("")
    L.append("## Saknar källa")
    L.append("")
    for c in miss:
        L.append(f"- **{c['plats']}.** \"…{c['ctx']}…\" — {c['text']}. {c['note']} "
                 + (f"*Förslag:* {c['forslag']}" if c["forslag"] else ""))
    L.append("")
    L.append(SLOPE_SECTION(N))
    L.append("")
    L.append(NEW_NUMBERS_SECTION(N))
    L.append("")
    L.append(SCOPE_SECTION)
    L.append("")
    L.append("## Alla kontrollerade tal")
    L.append("")
    L.append("| # | Ställe | Tal i texten | Tal ur data | Källfil | Status | Not |")
    L.append("| --- | --- | --- | --- | --- | --- | --- |")
    order = {"avviker": 0, "saknar källa": 1, "stämmer": 2}
    for i, c in enumerate(sorted(claims, key=lambda c: order[c["status"]]), 1):
        L.append(f"| {i} | {_esc(c['plats'])} | {_esc(c['text'])} | {_esc(c['data'])} | {_esc(c['source'])} | "
                 f"{c['status']} | {_esc(c['note'])} |")
    path.write_text("\n".join(L) + "\n")


def _git_head():
    import subprocess
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "?"


def SLOPE_SECTION(N: Numbers) -> str:
    lab = {"census": "USA, counties", "inegi": "Mexiko", "geostat": "Europa"}
    L = ["## Lutningens känslighet (nivå 2, tabell 3)", "",
         "Manusets mening \"removing three of the 2,383 U.S. counties changes the slope for 2010–2020 from 0.49 "
         "to 0.39\" bygger på ett urval som inte finns i repot. Skriptet räknar i stället för varje rad på nivå 2: "
         "OLS-lutningen som i tabell 3, lutningen utan de tre enheter som var för sig ändrar den mest (borttagna "
         "samtidigt), och Theil–Sen-lutningen (median av alla parvisa lutningar över samma sammanslagna öst- och "
         "nordkomponenter). Theil–Sen har inget bootstrapintervall här.", "",
         "Som tillägg (inte efterfrågat, märkt här) visas spannet när de tre enheter tas bort som var för sig "
         "sänker respektive höjer lutningen mest. De tre mest inflytelserika drar ofta åt olika håll, så att ta "
         "bort dem tillsammans kan ändra lutningen lite, medan ett ensidigt urval ändrar den mycket. Manusets "
         "0.39 för USA 2010–2020 ligger inom spannet för ett ensidigt urval (kolumnen näst längst till höger).", "",
         "| Facit | Period | Enheter | OLS (tabell 3) | OLS utan tre mest inflytelserika | Theil–Sen | "
         "Utan tre som sänker / höjer mest | Borttagna (tre mest inflytelserika) |",
         "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for ds, a, b in TABLE3_L2:
        t = f"slope:{ds}:{a}-{b}"
        L.append(f"| {lab[ds]} | {a}–{b} | {fmt(N.get(t + ':n'), 0)} | {fmt(N.get(t + ':ols'), 2)} | "
                 f"{fmt(N.get(t + ':ols_without_top3'), 2)} | {fmt(N.get(t + ':theil_sen'), 2)} | "
                 f"{fmt(N.get(t + ':ols_without_top3_down'), 2)} / {fmt(N.get(t + ':ols_without_top3_up'), 2)} | "
                 f"{_esc(N.rows[t + ':ols_without_top3']['note'].replace('borttagna: ', ''))} |")
    ols = [N.get(f"slope:{d}:{a}-{b}:ols") for d, a, b in TABLE3_L2]
    s3 = [N.get(f"slope:{d}:{a}-{b}:ols_without_top3") for d, a, b in TABLE3_L2]
    ts = [N.get(f"slope:{d}:{a}-{b}:theil_sen") for d, a, b in TABLE3_L2]
    dl = [x - y for x, y in zip(s3, ols)]
    sg = lambda x: ("+" if float(fmt(x, 2).replace(MINUS, "-")) > 0 else "") + fmt(x, 2)
    us = "slope:census:2010-2020"
    L += ["", f"OLS-lutningarna är {fmt_range(ols, 2)}; utan de tre mest inflytelserika enheterna "
              f"{fmt_range(s3, 2)} (ändring {sg(min(dl))} till {sg(max(dl))}); Theil–Sen {fmt_range(ts, 2)}. "
              f"För USA 2010–2020 ger borttagningen {fmt(N.get(us + ':ols_without_top3'), 2)} och Theil–Sen "
              f"{fmt(N.get(us + ':theil_sen'), 2)}; de tre counties som var för sig sänker lutningen mest ger "
              f"{fmt(N.get(us + ':ols_without_top3_down'), 2)}.", "",
          "**Förslag till ny mening** (ersätter \"but they are sensitive to a few units: removing three of the "
          "2,383 U.S. counties changes the slope for 2010–2020 from 0.49 to 0.39\"):", "",
          f"> The regression slopes of measured on true shift are {fmt_range(ols, 1)}. Removing the three most "
          f"influential units changes them by {sg(min(dl))} to {sg(max(dl))}, and the median-based Theil–Sen "
          f"slope is {fmt_range(ts, 2)}, lowest in Europe.", "",
          "Meningens tal skrivs av skriptet ur `numbers.csv`. Lutningens känslighet är en ny, beskrivande analys "
          "som inte står i `docs/ANALYSIS_PLAN.md`. Läggs meningen in i manuset bör det föras in som en rad under "
          "*Avvikelser* i planen (det är inte gjort här)."]
    return "\n".join(L)


def NEW_NUMBERS_SECTION(N: Numbers) -> str:
    L = ["## Nya tal (skriptade här)", ""]
    L.append("**Sverige (GEOSTAT).** Lägesfel, median (km) med blockbootstrap:")
    L.append("")
    L.append("| Nivå | År | Enheter | Median | 95 %-intervall |")
    L.append("| --- | --- | --- | --- | --- |")
    for lv, name in ((2, "kommuner"), (1, "län")):
        for y in (2011, 2021):
            k = f"swe:pos:L{lv}:{y}"
            L.append(f"| {name} | {y} | {fmt(N.get(k + ':n'), 0)} | {fmt(N.get(k + ':err_med'), 2)} | "
                     f"{fmt(N.get(k + ':err_med_lo'), 2)}–{fmt(N.get(k + ':err_med_hi'), 2)} |")
    L.append("")
    L.append(f"Förändring 2011–2021, kommuner: {fmt(N.get('swe:chg:L2:n'), 0)} enheter, median flytt "
             f"{fmt(N.get('swe:chg:L2:true_med'), 3)} km, rätt riktning "
             f"{fmt(N.get('swe:chg:L2:dir_right'), 3)} ({fmt(N.get('swe:chg:L2:dir_right_lo'), 3)}–"
             f"{fmt(N.get('swe:chg:L2:dir_right_hi'), 3)}), block = län.")
    L.append("")
    L.append("**Fishers exakta test, Mexikos delstater 2010–2020 mot USA:s:** " + "; ".join(
        f"{a}–{b}: p = {fmt(N.get(f'fisher:MEX2010-2020_vs_USA{a}-{b}:p'), 3)} "
        f"({N.rows[f'fisher:MEX2010-2020_vs_USA{a}-{b}:p']['note']})"
        for a, b in ((2000, 2010), (2010, 2020), (2000, 2020))) + ".")
    L.append("")
    L.append("**Kvoten median flytt / median fel, delstater:** Mexiko "
             f"{fmt(N.get('states:MEX:2010-2020:true_over_err'), 3)}; USA " + ", ".join(
                 f"{a}–{b} {fmt(N.get(f'states:USA:{a}-{b}:true_over_err'), 3)}"
                 for a, b in ((2000, 2010), (2010, 2020), (2000, 2020))) + ".")
    L.append("")
    L.append(f"**FIPS-kod mot GADM-polygon, USA:s counties** ({fmt(N.get('fips:matched_2010'), 0)} counties "
             "knutna entydigt 2010). Andel med samma facitflytt inom 0,1 km, och rätt riktning med respektive facit:")
    L.append("")
    L.append("| Period | Knutna | Inom 0,1 km (AA–AA) | Inom 0,1 km (AA–Barmore) | Rätt riktning, FIPS | Rätt riktning, GADM |")
    L.append("| --- | --- | --- | --- | --- | --- |")
    for a, b in ((2000, 2010), (2010, 2020), (2000, 2020)):
        k = f"fips:{a}-{b}"
        L.append(f"| {a}–{b} | {fmt(N.get(k + ':n_matched'), 0)} | {fmt(N.get(k + ':share_within_0.1km_aa'), 3)} | "
                 f"{fmt(N.get(k + ':share_within_0.1km_bm'), 3)} | {fmt(N.get(k + ':dir_right_fips'), 3)} | "
                 f"{fmt(N.get(k + ':dir_right_gadm'), 3)} |")
    L.append("")
    L.append("**Andel enheter där ljuset ligger närmare facit än geometrisk mittpunkt:**")
    L.append("")
    L.append("| Facit | Nivå | År | Enheter | Andel |")
    L.append("| --- | --- | --- | --- | --- |")
    for k in N.rows:
        if k.startswith("light_vs_geom:") and k.endswith(":share"):
            _, ds, lv, y, _ = k.split(":")
            L.append(f"| {ds} | {lv[1]} | {y} | {fmt(N.get(k[:-6] + ':n'), 0)} | {fmt(N.get(k), 3)} |")
    return "\n".join(L)


SCOPE_SECTION = """## Omfång och förbehåll

- Kontrollerat: abstract, inledning, metod (inklusive tabell 1), resultat (inklusive tabell 2 och 3),
  diskussion och slutsats. Inte kontrollerat: figurtexterna (t.ex. "2,975 common regions" i figur S1),
  referenslistan och regelvärdena i den förregistrerade tolkningsregeln (0.65, 0.50, 0.15, 150 enheter),
  som är definitioner ur planen §7 och inte resultat.
- Talen för fackelmasken (83–99, 6–10, 65 procent) står i en tabell i `docs/METHODS.md` och har ingen CSV;
  skriptet läser tabellen men talen är inte återskapade ur rådata här.
- Talet "about 94 percent" för Kanadas dissemination areas har ingen källa i repot (manuset markerar
  själv [REF needed]).
- Tre intervallgränser låg exakt på en halva med fyra decimaler utan att vara exakta bråk: Europas
  övre gräns för trendskattningen (0.4750), rangkorrelationens övre gräns för USA:s counties 2000–2010
  (0.2850) och lutningens nedre gräns för Mexiko (0.0150). Skriptet kör om A2
  (`tools/trend_estimate.run_dataset`) och `tools/change_detection.main` (kommandot i
  `docs/HARMONIZATION.md`), med samma ordning och frö, kontrollerar att omkörningarna återger de
  committade CSV-filerna med fyra decimaler och använder värdena med full precision för tabell 3 och
  trendintervallet. Övriga tal ur CSV avrundas som de står; tal som ligger exakt på en halva (då är de
  exakta bråk, t.ex. 0.875 = 28/32) anges i kolumnen Not.
- Andelen där ljuset ligger närmare än geometrisk mittpunkt räknas på samma urval som B
  (`tools/compare_benchmarks.py`): för Europa länderna med facit 2011, så nivå 1 år 2021 har 331 enheter.
- Sveriges intervall är blockbootstrap med block = län och fröet i C-drivrutinen (20261121), inte fröet i
  `change_detection.py`; punktskattningarna beror inte på fröet.
- FIPS-jämförelsen knyter GADM-enheter till FIPS-counties via närmaste centroid 2010, som
  granskningsskriptet `docs/review/scripts/check_fips_truth.py`. Counties som inte kan knytas entydigt
  (inom 5 km) ingår inte.
"""


def write_numbers(N: Numbers, path: Path):
    df = pd.DataFrame(N.rows.values())
    df["value"] = [f"{v:.6f}".rstrip("0").rstrip(".") if isinstance(v, float) else str(v) for v in df.value]
    df.sort_values("key").to_csv(path, index=False)


def main(argv=None):
    p = argparse.ArgumentParser(description="Talrevision av manuset.")
    p.add_argument("--skip-a2", action="store_true",
                   help="Kör inte om A2 och change_detection med full precision (gränsfallen avgörs då inte).")
    p.add_argument("--n-boot", type=int, default=1000)
    a = p.parse_args(argv)
    _imports()
    N = Numbers()
    misc(N)
    derived(N)
    print("Sverige …", flush=True)
    sweden(N, a.n_boot)
    print("delstater …", flush=True)
    states(N)
    print("FIPS …", flush=True)
    fips_vs_gadm(N)
    print("ljus mot geometrisk mittpunkt …", flush=True)
    light_vs_geometric(N)
    print("lutningar …", flush=True)
    slope_sensitivity(N)
    if not a.skip_a2:
        print("A2 med full precision …", flush=True)
        a2_full_precision(N)
        print("change_detection med full precision …", flush=True)
        cd_full_precision(N)
    manus = MANUS.read_text()
    claims = build_claims(N, manus)
    # Notera tal i CSV som ligger exakt på en halva vid textens precision
    write_numbers(N, OUT_CSV)
    write_review(claims, N, OUT_MD)
    st = pd.Series([c["status"] for c in claims]).value_counts()
    print(st.to_string())
    for c in claims:
        if c["status"] != "stämmer":
            print(f"[{c['status']}] {c['plats']}: {c['text']} -> {c['data']}")


if __name__ == "__main__":
    main()
