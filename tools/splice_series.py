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

Exempel
  python tools/splice_series.py --dmsp dmsp-irqr-reg --viirs viirs-dmsp-reg \\
      --name spliced-reg --results data/results_mask2km
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


def splice_level(dmsp: dict[int, pd.DataFrame], viirs: dict[int, pd.DataFrame], year0: int):
    """dmsp/viirs: år -> resultat-DataFrame. Returnerar år -> skarvad DataFrame."""
    d0 = dmsp[year0].set_index("unit_id")
    v0 = viirs[year0].set_index("unit_id")
    delta = pd.DataFrame(index=v0.index)
    for m in ("bm", "aa"):
        ok = d0[f"{m}_lon"].reindex(v0.index).notna() & v0[f"{m}_lon"].notna()
        de, dn = east_north_km(v0[f"{m}_lon"], v0[f"{m}_lat"],
                               d0[f"{m}_lon"].reindex(v0.index), d0[f"{m}_lat"].reindex(v0.index))
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


def main(argv=None):
    p = argparse.ArgumentParser(description="Skarva DMSP- och VIIRS-centroider per enhet.")
    p.add_argument("--dmsp", default="dmsp-irqr-reg")
    p.add_argument("--viirs", default="viirs-dmsp-reg")
    p.add_argument("--name", default="spliced-reg")
    p.add_argument("--year", type=int, default=2013, help="Skarvår (överlappsår).")
    p.add_argument("--results", type=Path, default=Path("data/results_mask2km"))
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
        for name in (args.dmsp, args.viirs):
            d = args.results / args.boundary_set / name / f"adm{level}"
            src[name] = {int(f.stem): read_results(f) for f in sorted(d.glob("*.csv"))}
            if args.year not in src[name]:
                sys.exit(f"{name} adm{level} saknar {args.year} i {d}")
        out, delta = splice_level(src[args.dmsp], src[args.viirs], args.year)
        if not args.no_flags:
            from wpc.run import flags, load_units
            units = load_units(args.boundaries, args.boundary_set, level, args.cache_dir)
            pos = {u: i for i, u in enumerate(units["unit_id"])}
            for y, df in out.items():
                if y <= args.year:
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
                "level": level, "units": int(len(delta)), "spliced": int(delta.splice_ok.sum()),
                "shift_km_bm": {"median": float(dist.median()), "p90": float(dist.quantile(0.9)),
                                "max": float(dist.max())},
                "years": sorted(out)}
        (od / "splice.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
        log(f"adm{level}: {meta['spliced']} av {meta['units']} enheter skarvade, "
            f"flytt median {meta['shift_km_bm']['median']:.2f} km -> {od} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
