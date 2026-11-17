#!/usr/bin/env python3
"""
World Population Centroids v2 – skarva DMSP- och VIIRS-centroider per enhet.

Den skarvade serien behåller DMSP-centroiderna till och med skarvåret (2013)
och flyttar varje enhets VIIRS-centroider för åren efter skarvåret med
enhetens skillnad mellan källorna skarvåret:

    Δ = DMSP(2013) − VIIRS(2013),  uttryckt i km österut och norrut vid
                                    VIIRS-punkten 2013,
    skarvad(t) = VIIRS(t) + Δ       för t > 2013.

Förflyttningen görs på en sfär (R = 6371 km) i det lokala planet vid
punkten: breddgraden ändras med Δnorr / R och längdgraden med
Δöst / (R cos φ). Skillnaden Δ mäts på samma sätt vid VIIRS-punkten 2013,
så att VIIRS(2013) + Δ = DMSP(2013) exakt.
Aboufadel–Austin- och Barmore-centroiderna skarvas var för sig.

Enheter som saknar centroid i någon av källorna skarvåret kan inte skarvas.
Deras VIIRS-år skrivs ändå, oförändrade, med splice_ok = 0, så att de kan
väljas bort. Flaggorna (aa_flag, bm_flag: 0 egen enhet, 1 utanför alla,
2 annan enhet) beräknas om för de flyttade punkterna.

Resultatet skrivs i samma format som pipelinen, i
<out>/<gränser>/<namn>/adm<nivå>/<år>.csv, med kolumnerna
splice_year, splice_ok, d_east_km_bm, d_north_km_bm, d_east_km_aa och
d_north_km_aa tillagda (0 för DMSP-åren).

Två sätt att uppskatta Δ (--method):
  point  skillnaden skarvåret, som ovan. Ett slumpmässigt hopp i DMSP eller
         VIIRS skarvåret följer med till alla senare år.
  trend  en gemensam modell per enhet och koordinat (km öster/norr): en
         rät linje för båda källorna och en konstant skillnad Δ mellan dem,
         anpassad med viktad minsta kvadrat till DMSP de --dmsp-window åren
         till och med skarvåret (standard 2004–2013) och VIIRS de
         --viirs-window åren från och med skarvåret (standard 2013–2025).
         Vikterna är 1/σ², med brusnivåernas kvot skattad ur medianerna av
         årshoppen i respektive källa. Det motsvarar hur
         brott korrigeras vid homogenisering av klimatserier: nivåskillnad
         med gemensam trend, skattad från flera år på var sida av brottet.
         Kräver minst tre år i varje källa. Med DMSP-brus dubbelt så stort
         som VIIRS-brus blir osäkerheten i Δ ungefär 45 % av den med ett enda
         år (1,01 mot 2,24 i enheter av VIIRS-bruset). Kortare fönster ger
         mindre vinst (5 + 5 år: 1,51), eftersom trenden och nivåskillnaden
         då är svåra att skilja åt; separata linjer per källa utvärderade
         skarvåret behåller 60 % av skarvårets slumpfel.

Exempel
  python tools/splice_series.py --dmsp dmsp-irqr-reg --viirs viirs-dmsp-reg \\
      --name spliced-reg --results data/results_mask2km
  python tools/splice_series.py --method trend --name spliced-reg-trend
  python tools/splice_series.py --direction backward --viirs viirs \\
      --viirs-results data/results --name spliced-back

Riktning (--direction): forward (standard förut) behåller DMSP och flyttar
VIIRS; backward behåller VIIRS från och med skarvåret och flyttar DMSP-åren
med Δ = VIIRS(skarvår) − DMSP(skarvår). Valideringen mot folkräkningar
(avsnitt 10 i docs/HARMONIZATION.md) visade att backward med rå VIIRS ger
klart mindre lägesfel alla år.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.progress import log  # noqa: E402

R_KM = 6371.0


def read_results(path) -> pd.DataFrame:
    """Läs en resultatfil. ID-kolumner läses som text utan NA-tolkning: GADM
    4.1 har en enhet med ID:t "NA", som pandas annars gör till ett saknat värde."""
    df = pd.read_csv(path, keep_default_na=False, na_values={})
    for c in df.columns:
        if c in ("unit_id", "source", "COUNTRY") or c.startswith(("GID_", "NAME_", "TYPE_", "ENGTYPE_")):
            continue
        if not pd.api.types.is_numeric_dtype(df[c]):
            try:
                df[c] = pd.to_numeric(df[c].replace("", np.nan))
            except (ValueError, TypeError):
                pass
    return df


def east_north_km(lon0, lat0, lon1, lat1):
    """Förflyttning från punkt 0 till punkt 1 i km (öster, norr), lokalt plan vid punkt 0."""
    dlon = (np.asarray(lon1) - np.asarray(lon0) + 180.0) % 360.0 - 180.0
    de = R_KM * np.radians(dlon) * np.cos(np.radians(lat0))
    dn = R_KM * np.radians(np.asarray(lat1) - np.asarray(lat0))
    return de, dn


def move(lon, lat, de, dn):
    """Flytta punkter de km österut och dn km norrut, i det lokala planet vid
    startpunkten (samma plan som east_north_km, så att move(p0, Δ(p0, p1)) = p1)."""
    lat = np.asarray(lat, dtype=float)
    lat2 = np.clip(lat + np.degrees(np.asarray(dn) / R_KM), -89.999, 89.999)
    lon2 = np.asarray(lon) + np.degrees(np.asarray(de) / (R_KM * np.cos(np.radians(lat))))
    lon2 = (lon2 + 180.0) % 360.0 - 180.0
    return lon2, lat2


def _coords(years: dict[int, pd.DataFrame], use: list[int], ref: pd.DataFrame, m: str):
    """Lägen i km (öster, norr) från referenspunkten, matris enheter × år."""
    E = np.full((len(ref), len(use)), np.nan)
    N = np.full((len(ref), len(use)), np.nan)
    for j, y in enumerate(use):
        d = years[y].set_index("unit_id").reindex(ref.index)
        E[:, j], N[:, j] = east_north_km(ref[f"{m}_lon"], ref[f"{m}_lat"], d[f"{m}_lon"], d[f"{m}_lat"])
    return E, N


def noise_ratio(dmsp: dict, viirs: dict, m: str = "bm") -> float:
    """Förhållandet mellan brusnivåerna (DMSP / VIIRS), skattat som kvoten
    mellan medianerna av alla årshopp i respektive källa. Används som vikt i
    joint_shift (vikter 1/σ²)."""
    def med(years):
        ys = sorted(years)
        out = []
        for a, b in zip(ys, ys[1:]):
            if b != a + 1:
                continue
            A = years[a].set_index("unit_id")
            B = years[b].set_index("unit_id").reindex(A.index)
            d = np.hypot(*east_north_km(A[f"{m}_lon"], A[f"{m}_lat"], B[f"{m}_lon"], B[f"{m}_lat"]))
            out.append(np.asarray(d)[np.isfinite(d)])
        return float(np.median(np.concatenate(out)))
    return med(dmsp) / med(viirs)


def joint_shift(dmsp, viirs, dy, vy, ref, m, year0, min_years, ratio: float = 1.0):
    """Δ per enhet ur en gemensam modell för båda källorna.

    För varje enhet och koordinat (km öster respektive norr):
        DMSP(t) = a + b (t − year0) + e,      t i dy
        VIIRS(t) = a − Δ + b (t − year0) + e, t i vy
    dvs. en gemensam linjär trend och en konstant skillnad Δ mellan
    sensorerna, skattad med viktad minsta kvadrat över alla år (vikt 1 för
    VIIRS och 1/ratio² för DMSP, där ratio är DMSP-brusets storlek i
    förhållande till VIIRS-brusets). Det är samma modell som vid
    homogenisering av klimatserier med ett känt brott: nivåskillnad med
    gemensam trend, skattad från flera år på var sida.
    """
    Ed, Nd = _coords(dmsp, dy, ref, m)
    Ev, Nv = _coords(viirs, vy, ref, m)
    td = np.array(dy, float) - year0
    tv = np.array(vy, float) - year0
    out = np.full((len(ref), 2), np.nan)
    for i in range(len(ref)):
        kd = np.isfinite(Ed[i]) & np.isfinite(Nd[i])
        kv = np.isfinite(Ev[i]) & np.isfinite(Nv[i])
        if kd.sum() < min_years or kv.sum() < min_years:
            continue
        A = np.vstack([np.column_stack([np.ones(kd.sum()), td[kd], np.zeros(kd.sum())]),
                       np.column_stack([np.ones(kv.sum()), tv[kv], -np.ones(kv.sum())])])
        sw = np.concatenate([np.full(kd.sum(), 1.0 / ratio), np.ones(kv.sum())])   # sqrt(vikt)
        for c, (Y1, Y2) in enumerate(((Ed[i], Ev[i]), (Nd[i], Nv[i]))):
            y = np.concatenate([Y1[kd], Y2[kv]])
            out[i, c] = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)[0][2]
    return out[:, 0], out[:, 1]


def splice_level(dmsp: dict[int, pd.DataFrame], viirs: dict[int, pd.DataFrame], year0: int,
                 method: str = "point", dmsp_window: int = 10, viirs_window: int = 13,
                 min_years: int = 3, ratio: float | None = None):
    """dmsp/viirs: år -> resultat-DataFrame. Returnerar år -> skarvad DataFrame.

    method="point": Δ = DMSP(year0) − VIIRS(year0).
    method="trend": Δ ur en gemensam modell (joint_shift) med DMSP
    year0−dmsp_window+1 … year0 och VIIRS year0 … year0+viirs_window−1 (minst
    `min_years` år med centroid i varje källa). ratio = DMSP-brus / VIIRS-
    brus; None = skattas ur data (noise_ratio).
    """
    if method == "trend" and ratio is None:
        ratio = noise_ratio(dmsp, viirs)
    d0 = dmsp[year0].set_index("unit_id")
    v0 = viirs[year0].set_index("unit_id")
    delta = pd.DataFrame(index=v0.index)
    for m in ("bm", "aa"):
        if method == "point":
            ok = d0[f"{m}_lon"].reindex(v0.index).notna() & v0[f"{m}_lon"].notna()
            de, dn = east_north_km(v0[f"{m}_lon"], v0[f"{m}_lat"],
                                   d0[f"{m}_lon"].reindex(v0.index), d0[f"{m}_lat"].reindex(v0.index))
        elif method == "trend":
            dy = [y for y in range(year0 - dmsp_window + 1, year0 + 1) if y in dmsp]
            vy = [y for y in range(year0, year0 + viirs_window) if y in viirs]
            de, dn = joint_shift(dmsp, viirs, dy, vy, v0, m, year0, min_years, ratio)
            ok = np.isfinite(de) & np.isfinite(dn) & v0[f"{m}_lon"].notna().values
        else:
            raise ValueError(method)
        delta[f"d_east_km_{m}"] = np.where(ok, de, np.nan)
        delta[f"d_north_km_{m}"] = np.where(ok, dn, np.nan)
    delta["splice_ok"] = delta[["d_east_km_bm", "d_east_km_aa"]].notna().all(axis=1).astype(int)

    out = {}
    for y, df in dmsp.items():
        if y > year0:
            continue
        df = df.copy()
        for c in ("d_east_km_bm", "d_north_km_bm", "d_east_km_aa", "d_north_km_aa"):
            df[c] = 0.0
        df["splice_ok"] = 1
        out[y] = df
    for y, df in viirs.items():
        if y <= year0:
            continue
        df = df.set_index("unit_id").join(delta, how="left")
        df["splice_ok"] = df["splice_ok"].fillna(0).astype(int)
        for m in ("bm", "aa"):
            de = df[f"d_east_km_{m}"].fillna(0.0).values
            dn = df[f"d_north_km_{m}"].fillna(0.0).values
            lon, lat = move(df[f"{m}_lon"].values, df[f"{m}_lat"].values, de, dn)
            keep = df[f"{m}_lon"].isna().values
            df[f"{m}_lon"] = np.where(keep, np.nan, lon)
            df[f"{m}_lat"] = np.where(keep, np.nan, lat)
        out[y] = df.reset_index()
    return out, delta


def splice_backward(dmsp: dict[int, pd.DataFrame], viirs: dict[int, pd.DataFrame], year0: int):
    """Omvänd skarvning: VIIRS behålls från och med year0, och varje enhets
    DMSP-centroider före year0 flyttas med Δ = VIIRS(year0) − DMSP(year0),
    mätt vid DMSP-punkten year0. Valideringen mot folkräkningar visade att
    DMSP:s fel per enhet till stor del är systematiskt och bestående, så att
    VIIRS noggrannhet på det sättet förs över till DMSP-åren
    (docs/HARMONIZATION.md, avsnitt 10)."""
    d0 = dmsp[year0].set_index("unit_id")
    v0 = viirs[year0].set_index("unit_id")
    delta = pd.DataFrame(index=d0.index)
    for m in ("bm", "aa"):
        vv = v0.reindex(d0.index)
        ok = d0[f"{m}_lon"].notna() & vv[f"{m}_lon"].notna()
        de, dn = east_north_km(d0[f"{m}_lon"], d0[f"{m}_lat"], vv[f"{m}_lon"], vv[f"{m}_lat"])
        delta[f"d_east_km_{m}"] = np.where(ok, de, np.nan)
        delta[f"d_north_km_{m}"] = np.where(ok, dn, np.nan)
    delta["splice_ok"] = delta[["d_east_km_bm", "d_east_km_aa"]].notna().all(axis=1).astype(int)
    out = {}
    for y, df in viirs.items():
        if y < year0:
            continue
        df = df.copy()
        for c in ("d_east_km_bm", "d_north_km_bm", "d_east_km_aa", "d_north_km_aa"):
            df[c] = 0.0
        df["splice_ok"] = 1
        out[y] = df
    for y, df in dmsp.items():
        if y >= year0:
            continue
        df = df.set_index("unit_id").join(delta, how="left")
        df["splice_ok"] = df["splice_ok"].fillna(0).astype(int)
        for m in ("bm", "aa"):
            lon, lat = move(df[f"{m}_lon"].values, df[f"{m}_lat"].values,
                            df[f"d_east_km_{m}"].fillna(0.0).values, df[f"d_north_km_{m}"].fillna(0.0).values)
            keep = df[f"{m}_lon"].isna().values
            df[f"{m}_lon"] = np.where(keep, np.nan, lon)
            df[f"{m}_lat"] = np.where(keep, np.nan, lat)
        out[y] = df.reset_index()
    return out, delta


def main(argv=None):
    p = argparse.ArgumentParser(description="Skarva DMSP- och VIIRS-centroider per enhet.")
    p.add_argument("--dmsp", default="dmsp-irqr-reg")
    p.add_argument("--viirs", default="viirs-dmsp-reg")
    p.add_argument("--name", default="spliced-reg")
    p.add_argument("--year", type=int, default=2013, help="Skarvår (överlappsår).")
    p.add_argument("--method", choices=["point", "trend"], default="point",
                   help="point: skillnaden skarvåret; trend: skillnaden mellan linjära trender "
                        "(DMSP före, VIIRS efter skarvåret), utvärderade skarvåret.")
    p.add_argument("--dmsp-window", type=int, default=10,
                   help="DMSP-år till och med skarvåret för --method trend (standard 10: 2004–2013).")
    p.add_argument("--viirs-window", type=int, default=13,
                   help="VIIRS-år från och med skarvåret för --method trend (standard 13: 2013–2025).")
    p.add_argument("--direction", choices=["forward", "backward"], default="forward",
                   help="forward: VIIRS flyttas till DMSP (DMSP behålls); backward: DMSP flyttas "
                        "till VIIRS (VIIRS behålls från skarvåret).")
    p.add_argument("--results", type=Path, default=Path("data/results_mask2km"))
    p.add_argument("--viirs-results", type=Path, default=None,
                   help="Resultatmapp för VIIRS-källan, om den inte är --results.")
    p.add_argument("--out", type=Path, default=None, help="Standard: samma som --results.")
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--boundaries", type=Path, default=Path("data/raw/gadm41/gadm_410-levels.gpkg"))
    p.add_argument("--cache-dir", type=Path, default=Path("data/processed/boundaries"))
    p.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--no-flags", action="store_true", help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    out_root = (args.out or args.results) / args.boundary_set / args.name

    for level in args.levels:
        t0 = time.time()
        src = {}
        for name, base in ((args.dmsp, args.results), (args.viirs, args.viirs_results or args.results)):
            d = base / args.boundary_set / name / f"adm{level}"
            src[name] = {int(f.stem): read_results(f) for f in sorted(d.glob("*.csv"))}
            if args.year not in src[name]:
                sys.exit(f"{name} adm{level} saknar {args.year} i {d}")
        ratio = noise_ratio(src[args.dmsp], src[args.viirs]) if args.method == "trend" else None
        if args.direction == "backward":
            if args.method != "point":
                sys.exit("--direction backward finns bara med --method point")
            out, delta = splice_backward(src[args.dmsp], src[args.viirs], args.year)
            moved = {y for y in out if y < args.year}
        else:
            out, delta = splice_level(src[args.dmsp], src[args.viirs], args.year, method=args.method,
                                      dmsp_window=args.dmsp_window, viirs_window=args.viirs_window,
                                      ratio=ratio)
            moved = {y for y in out if y > args.year}
        if not args.no_flags:
            from wpc.run import flags, load_units
            units = load_units(args.boundaries, args.boundary_set, level, args.cache_dir)
            pos = {u: i for i, u in enumerate(units["unit_id"])}
            for y, df in out.items():
                if y not in moved:
                    continue
                idx = df["unit_id"].map(pos).values
                for m in ("aa", "bm"):
                    lon = np.full(len(units), np.nan)
                    lat = np.full(len(units), np.nan)
                    lon[idx] = df[f"{m}_lon"].values
                    lat[idx] = df[f"{m}_lat"].values
                    df[f"{m}_flag"] = flags(units, lon, lat)[idx]
        od = out_root / f"adm{level}"
        od.mkdir(parents=True, exist_ok=True)
        for y, df in sorted(out.items()):
            df = df.assign(source=args.name, splice_year=args.year)
            df.to_csv(od / f"{y}.csv", index=False, float_format="%.8f")
        dsum = delta.dropna(subset=["d_east_km_bm"])
        dist = np.hypot(dsum.d_east_km_bm, dsum.d_north_km_bm)
        meta = {"name": args.name, "dmsp": args.dmsp, "viirs": args.viirs, "splice_year": args.year,
                "direction": args.direction, "method": args.method,
                "dmsp_window": args.dmsp_window if args.method == "trend" else None,
                "viirs_window": args.viirs_window if args.method == "trend" else None,
                "noise_ratio_dmsp_viirs": ratio,
                "level": level, "units": int(len(delta)), "spliced": int(delta.splice_ok.sum()),
                "shift_km_bm": {"median": float(dist.median()), "p90": float(dist.quantile(0.9)),
                                "max": float(dist.max())},
                "years": sorted(out)}
        (od / "splice.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
        log(f"adm{level}: {meta['spliced']} av {meta['units']} enheter skarvade, "
            f"flytt median {meta['shift_km_bm']['median']:.2f} km -> {od} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
