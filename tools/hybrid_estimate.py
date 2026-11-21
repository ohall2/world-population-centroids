#!/usr/bin/env python3
"""A1 – Hybridskattning (docs/ANALYSIS_PLAN.md §3).

Fråga: tillför ljuset något när befolkningstal per underenhet (nivå 2) är kända?

För en nivå 1-enhet och ett år jämförs tre skattningar mot facit:
  ljus         utgåvans ljusviktade centroid på nivå 1 (nuvarande).
  hybrid-ljus  3D-medel (Aboufadel–Austin) av nivå 2-enheternas ljuscentroider,
               viktat med nivå 2-enheternas befolkning i facit samma år.
  hybrid-geom  samma vikter, men nivå 2-enheternas geometriska ytcentroider
               (ytriktig projektion). Skiljer ut vad ljuset bidrar med.

Alla nivå 2-enheter med befolkning > 0 ingår. En nivå 1-enhet används om de
nivå 2-enheter som har ljuscentroid täcker minst 99 % av dess facitbefolkning.

Mått:
  Läge       median great circle-fel mot facit per facit och år, med
             blockbootstrap, samt parvis felskillnad per enhet (ljus − hybrid-ljus
             och hybrid-geom − hybrid-ljus) med intervall.
  Förändring de frysta måtten i tools/change_detection.py (metrics, en, block_ids,
             bootstrap), på nivå 1, för perioderna i planen §2.

Hypoteser: H1a hybrid-ljus har mindre lägesfel än ljus; H1b hybrid-ljus ger rätt
riktning oftare än ljus; H1c hybrid-ljus har mindre lägesfel än hybrid-geom.

Förbehåll (ska stå i texten): hybriden använder facits egna befolkningstal per
nivå 2-enhet. Den prövar hur mycket fel som återstår från läget inom enheterna,
inte en oberoende skattning.

Facit: USA (census, utan Puerto Rico), Mexiko (inegi), Brasilien (ibge),
Europa (geostat). Kanadas facit rörs inte.

Utdata i --out:
  <ds>_hybrid_position.csv       lägesfel per år och skattning, med intervall
  <ds>_hybrid_position_pairs.csv parvis felskillnad per år, med intervall
  <ds>_hybrid_change.csv         förändringsmått per period och skattning
  <ds>_hybrid_coverage.csv       täckning och urval per nivå 1-enhet och år

Exempel
  python tools/hybrid_estimate.py --datasets census inegi ibge geostat
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import pyogrio
import shapely

from wpc.centroids import aboufadel_austin, great_circle_m
from tools.change_detection import en, metrics, block_ids, bootstrap, truth, release

SPHERE_R = 6371000.0
RNG = np.random.default_rng(20261121)

# Perioder och länder per facit (planen §1–2). iso=None betyder alla länder i
# facit (Europa). Brasilien har bara ett år: läge men ingen förändring.
DATASETS = {
    "census":  {"iso": ["USA"], "years": [2000, 2010, 2020],
                "pairs": [(2000, 2010), (2010, 2020), (2000, 2020)]},
    "inegi":   {"iso": ["MEX"], "years": [2010, 2020], "pairs": [(2010, 2020)]},
    "ibge":    {"iso": ["BRA"], "years": [2022], "pairs": []},
    "geostat": {"iso": None,    "years": [2011, 2021], "pairs": [(2011, 2021)]},
}


# --- rena funktioner (prövas i tests/test_hybrid_estimate.py) --------------

def parent_adm1(unit_id: str) -> str:
    """GADM nivå 2 -> nivå 1, t.ex. 'USA.5.13_1' -> 'USA.5_1'."""
    iso, a1, _a2 = unit_id.split("_", 1)[0].split(".")
    return f"{iso}.{a1}_1"


def mean_centroid_3d(lon, lat, w) -> tuple[float, float]:
    """Viktat medel av 3D-enhetsvektorer på sfären (Aboufadel–Austin)."""
    return aboufadel_austin(np.asarray(lon, float), np.asarray(lat, float),
                            np.asarray(w, float))


def great_circle_km(lon1, lat1, lon2, lat2):
    """Storcirkelavstånd i km på sfären med R = 6371 km."""
    return great_circle_m(lon1, lat1, lon2, lat2) / 1000.0


def covered_fraction(pop, has_light) -> float:
    """Andel av befolkningen som ligger i nivå 2-enheter med ljuscentroid."""
    pop = np.asarray(pop, float)
    has = np.asarray(has_light, bool)
    tot = pop.sum()
    return float(pop[has].sum() / tot) if tot > 0 else 0.0


def geometric_centroids(gdf: gpd.GeoDataFrame, iso: str) -> tuple[np.ndarray, np.ndarray]:
    """Geometrisk ytcentroid per enhet i ytriktig projektion, i enhetens ordning.

    USA: longitud skalad med cos(φ) per enhet (cylindrisk ytriktig approximation;
    landet är för stort för en gemensam projektion). Övriga länder: LAEA kring
    landets mittpunkt. Samma metod som granskningens geometriska baslinje
    (docs/review/scripts/baseline_geometric_centroid.py)."""
    g = gdf.geometry
    if iso == "USA":
        b = g.bounds
        lat0 = ((b.miny + b.maxy) / 2).to_numpy()
        cx, cy = [], []
        for geom, l0 in zip(g.to_numpy(), lat0):
            k = np.cos(np.radians(l0))
            gg = shapely.transform(geom, lambda xy: np.column_stack([xy[:, 0] * k, xy[:, 1]]))
            c = gg.centroid
            cx.append(c.x / k)
            cy.append(c.y)
        return np.asarray(cx), np.asarray(cy)
    c0 = g.union_all().centroid
    crs = f"+proj=laea +lat_0={c0.y} +lon_0={c0.x} +R={SPHERE_R} +units=m +no_defs"
    cen = g.to_crs(crs).centroid.to_crs(4326)
    return cen.x.to_numpy(), cen.y.to_numpy()


def hybrid_level1(adm2: pd.DataFrame) -> pd.DataFrame:
    """Aggregera nivå 2 till nivå 1.

    adm2 har kolumnerna unit_id, adm1, pop (>0), has_light, L_lon, L_lat
    (ljuscentroid), g_lon, g_lat (geometrisk centroid). Per nivå 1-enhet ges
    hybrid-ljus (hyl_*) och hybrid-geometri (hyg_*), räknade på de nivå 2-enheter
    som har ljuscentroid, viktade med befolkningen; cover = täckningsandel; used =
    cover >= 0.99; n_adm2, n_light = antal nivå 2-enheter totalt och med ljus.
    Saknar någon av de ljussatta enheterna geometri blir hyg_* NaN."""
    rows = []
    for a1, grp in adm2.groupby("adm1"):
        cover = covered_fraction(grp["pop"].to_numpy(), grp["has_light"].to_numpy())
        lit = grp[grp["has_light"].astype(bool)]
        row = {"adm1": a1, "cover": cover, "used": cover >= 0.99,
               "n_adm2": len(grp), "n_light": len(lit)}
        if len(lit):
            w = lit["pop"].to_numpy()
            hyl = mean_centroid_3d(lit["L_lon"].to_numpy(), lit["L_lat"].to_numpy(), w)
            hyg = mean_centroid_3d(lit["g_lon"].to_numpy(), lit["g_lat"].to_numpy(), w)
            row.update(hyl_lon=hyl[0], hyl_lat=hyl[1], hyg_lon=hyg[0], hyg_lat=hyg[1])
        else:
            row.update(hyl_lon=np.nan, hyl_lat=np.nan, hyg_lon=np.nan, hyg_lat=np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


# --- inläsning -------------------------------------------------------------

def level2_pop(truth_root: Path, ds: str, year: int, isos: list[str] | None) -> pd.DataFrame:
    """Nivå 2-enheter med befolkning > 0 i facit, med nivå 1-förälder."""
    d = pd.read_csv(truth_root / ds / str(year) / "adm2.csv", keep_default_na=False)
    d["pop"] = pd.to_numeric(d["pop"], errors="coerce")
    d = d[d["pop"] > 0].copy()
    if isos is not None:
        d = d[d.unit_id.str[:3].isin(isos)]
    d["adm1"] = d.unit_id.map(parent_adm1)
    return d[["unit_id", "adm1", "pop"]]


def level1_truth(truth_root: Path, ds: str, year: int, isos: list[str] | None,
                 min_points: int) -> pd.DataFrame:
    """Nivå 1-facit (complete, minst min_points punkter), iso-filtrerat."""
    t = truth(truth_root, ds, year, 1, min_points)
    if isos is not None:
        t = t[t.index.str[:3].isin(isos)]
    return t[["bm_lon", "bm_lat"]]


def geometry_table(gadm: Path, isos: list[str]) -> pd.DataFrame:
    """Geometrisk ytcentroid per nivå 2-enhet för angivna länder."""
    where = "GID_0 IN (" + ",".join(f"'{i}'" for i in isos) + ")"
    gdf = pyogrio.read_dataframe(gadm, columns=["GID_0", "GID_2"], where=where)
    parts = []
    for iso, g in gdf.groupby("GID_0"):
        glon, glat = geometric_centroids(g, iso)
        parts.append(pd.DataFrame({"unit_id": g["GID_2"].to_numpy(),
                                   "g_lon": glon, "g_lat": glat}))
    return pd.concat(parts, ignore_index=True)


# --- blockbootstrap för lägesmåtten ---------------------------------------

def _members(blocks: np.ndarray):
    ub, inv = np.unique(blocks, return_inverse=True)
    return [np.flatnonzero(inv == k) for k in range(len(ub))], len(ub)


def boot_ci(fn, blocks: np.ndarray, n_boot: int) -> tuple[float, float]:
    """95 %-intervall för en statistik fn(index) med blockbootstrap."""
    members, nb = _members(blocks)
    vals = []
    for _ in range(n_boot):
        pick = RNG.integers(0, nb, nb)
        idx = np.concatenate([members[k] for k in pick])
        vals.append(fn(idx))
    return float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


# --- per facit -------------------------------------------------------------

EST = ["light", "hybrid_light", "hybrid_geom"]
COORD = {"light": ("L_lon", "L_lat"), "hybrid_light": ("hyl_lon", "hyl_lat"),
         "hybrid_geom": ("hyg_lon", "hyg_lat")}


def per_year_table(cfg, ds, year, truth_root, release_path, geo, min_points):
    """DataFrame indexerad på nivå 1-enhet med facit- och skattningskoordinater
    för de enheter som ingår i A1 det året (complete nivå 1-facit, ljus på nivå 1,
    och en hybrid som används enligt 99 %-regeln)."""
    isos = cfg["iso"]
    t1 = level1_truth(truth_root, ds, year, isos, min_points)          # facit nivå 1
    rel1 = release(release_path, 1, {year})
    rel1 = rel1[rel1.year == year].set_index("unit_id")[["bm_lon", "bm_lat"]]  # ljus nivå 1

    pop = level2_pop(truth_root, ds, year, isos)                        # nivå 2-befolkning
    rel2 = release(release_path, 2, {year})
    rel2 = rel2[rel2.year == year].set_index("unit_id")[["bm_lon", "bm_lat"]]
    adm2 = pop.join(rel2.rename(columns={"bm_lon": "L_lon", "bm_lat": "L_lat"}), on="unit_id")
    adm2 = adm2.join(geo.set_index("unit_id"), on="unit_id")
    adm2["has_light"] = adm2["L_lon"].notna()

    hyb = hybrid_level1(adm2).set_index("adm1")

    E = t1.rename(columns={"bm_lon": "tru_lon", "bm_lat": "tru_lat"}).copy()
    E = E.join(rel1.rename(columns={"bm_lon": "L_lon", "bm_lat": "L_lat"}))
    E = E.join(hyb[["cover", "used", "n_adm2", "n_light",
                    "hyl_lon", "hyl_lat", "hyg_lon", "hyg_lat"]])
    E["used"] = E["used"].fillna(False).astype(bool)
    # Täckning rapporteras för de nivå 1-enheter som har facit (kandidaterna som
    # 99 %-regeln gallrar), så urvalet stämmer med analysmängden.
    cov = hyb.reindex(t1.index)
    cov.index.name = "adm1"
    cov["used"] = cov["used"].fillna(False).astype(bool)
    cov["cover"] = cov["cover"].fillna(0.0)
    cov[["n_adm2", "n_light"]] = cov[["n_adm2", "n_light"]].fillna(0).astype(int)
    # A1-mängden: facit på nivå 1 (finns i E), ljus på nivå 1, hybrid som används.
    E = E[E["used"] & E["L_lon"].notna() & E["hyl_lon"].notna()]
    return E, cov


def err_to_truth(E: pd.DataFrame, est: str) -> np.ndarray:
    lon, lat = COORD[est]
    return great_circle_km(E[lon].to_numpy(), E[lat].to_numpy(),
                           E["tru_lon"].to_numpy(), E["tru_lat"].to_numpy())


def run_dataset(ds, cfg, truth_root, release_path, gadm, min_points, n_boot):
    # Länderna (och därmed geometrin) härleds ur alla års nivå 2-facit, så att
    # även länder som bara finns med vissa år (t.ex. i geostat) kommer med.
    need = set()
    for y in cfg["years"]:
        need |= set(level2_pop(truth_root, ds, y, cfg["iso"]).unit_id)
    isos = sorted(set(u[:3] for u in need))
    geo = geometry_table(gadm, isos)
    missing_iso = sorted(set(u[:3] for u in need) - set(geo.unit_id.str[:3]))
    if missing_iso:
        print(f"  OBS {ds}: ingen geometri för {missing_iso} (hybrid-geom utgår där)")

    years = {y: per_year_table(cfg, ds, y, truth_root, release_path, geo, min_points)
             for y in cfg["years"]}

    pos_rows, pair_rows, cov_rows, chg_rows = [], [], [], []

    for y, (E, hyb) in years.items():
        cov = hyb.reset_index()
        cov.insert(0, "year", y)
        cov.insert(0, "dataset", ds)
        cov_rows.append(cov[["dataset", "year", "adm1", "cover", "used", "n_adm2", "n_light"]])

        multi = E.index.str[:3].nunique() >= 5
        blk = block_ids(E.index, 1, multi)
        for est in EST:
            err = err_to_truth(E, est)
            m = np.isfinite(err)
            if m.sum() == 0:
                continue
            lo, hi = boot_ci(lambda idx, e=err: np.median(e[idx]), blk[m], n_boot) \
                if m.all() else boot_ci(lambda idx, e=err[m]: np.median(e[idx]),
                                        block_ids(E.index[m], 1, multi), n_boot)
            pos_rows.append({"dataset": ds, "year": y, "estimate": est, "n": int(m.sum()),
                             "err_med": float(np.median(err[m])), "err_med_lo": lo, "err_med_hi": hi})

        # parvis felskillnad: ljus − hybrid-ljus, hybrid-geom − hybrid-ljus
        for a, b in [("light", "hybrid_light"), ("hybrid_geom", "hybrid_light")]:
            ea, eb = err_to_truth(E, a), err_to_truth(E, b)
            m = np.isfinite(ea) & np.isfinite(eb)
            if m.sum() == 0:
                continue
            diff = ea[m] - eb[m]
            blkm = block_ids(E.index[m], 1, multi)
            lo, hi = boot_ci(lambda idx, d=diff: np.median(d[idx]), blkm, n_boot)
            pair_rows.append({"dataset": ds, "year": y, "pair": f"{a}_minus_{b}",
                              "n": int(m.sum()), "diff_med": float(np.median(diff)),
                              "diff_lo": lo, "diff_hi": hi})

    # förändring: frysta måtten per period och skattning
    for y1, y2 in cfg["pairs"]:
        E1, E2 = years[y1][0], years[y2][0]
        ids = E1.index.intersection(E2.index)
        if len(ids) == 0:
            continue
        A, B = E1.loc[ids], E2.loc[ids]
        multi = ids.str[:3].nunique() >= 5
        t = np.column_stack(en(A.tru_lon, A.tru_lat, B.tru_lon, B.tru_lat))
        for est in EST:
            lon, lat = COORD[est]
            ok = A[lon].notna().to_numpy() & B[lon].notna().to_numpy()
            if ok.sum() == 0:
                continue
            idk = ids[ok]
            e = np.column_stack(en(A.loc[idk, lon], A.loc[idk, lat],
                                   B.loc[idk, lon], B.loc[idk, lat]))
            tk = t[ok]
            m = metrics(tk, e)
            ci = bootstrap(tk, e, block_ids(idk, 1, multi), n_boot)
            chg_rows.append({"dataset": ds, "series": est, "level": 1, "from": y1, "to": y2,
                             "n": len(idk), "n_blocks": len(np.unique(block_ids(idk, 1, multi))),
                             **m, **ci})

    return (pd.DataFrame(pos_rows), pd.DataFrame(pair_rows),
            pd.concat(cov_rows, ignore_index=True), pd.DataFrame(chg_rows))


def main(argv=None):
    p = argparse.ArgumentParser(description="A1 – hybridskattning av nivå 1-centroider.")
    p.add_argument("--datasets", nargs="+", default=list(DATASETS),
                   choices=list(DATASETS))
    p.add_argument("--release", type=Path, default=Path("data/release/v2.0.0"))
    p.add_argument("--truth", type=Path, default=Path("data/validation"))
    p.add_argument("--gadm", type=Path,
                   default=Path("data/processed/boundaries/gadm41_adm2.gpkg"))
    p.add_argument("--min-points", type=int, default=10)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)

    a.out.mkdir(parents=True, exist_ok=True)
    for ds in a.datasets:
        print(f"== {ds} ==")
        pos, pair, cov, chg = run_dataset(ds, DATASETS[ds], a.truth, a.release,
                                          a.gadm, a.min_points, a.n_boot)
        pos.to_csv(a.out / f"{ds}_hybrid_position.csv", index=False, float_format="%.4f")
        pair.to_csv(a.out / f"{ds}_hybrid_position_pairs.csv", index=False, float_format="%.4f")
        cov.to_csv(a.out / f"{ds}_hybrid_coverage.csv", index=False, float_format="%.4f")
        chg.to_csv(a.out / f"{ds}_hybrid_change.csv", index=False, float_format="%.4f")
        if len(pos):
            print(pos.round(3).to_string(index=False))
        if len(pair):
            print(pair.round(3).to_string(index=False))
        if len(chg):
            cols = ["series", "from", "to", "n", "err_med", "ratio", "rho",
                    "slope", "dir_right", "dir_right_lo", "dir_right_hi", "snr_gt1"]
            print(chg[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
