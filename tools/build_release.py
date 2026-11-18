#!/usr/bin/env python3
"""
World Population Centroids v2 – bygg den publicerbara datamängden.

Huvudserien är den bakåt skarvade (`spliced-back`, tools/splice_series.py
--direction backward): samregistrerad DMSP-IRQR 1992–2012 flyttad per
enhet till rå VIIRS 2013, och rå VIIRS VNL v2 2013–2025. Valet och
valideringen beskrivs i docs/HARMONIZATION.md (avsnitt 10).

Nivå 0 (länder) tas från körningen där ljuset norr om 65° N är maskat, som i
v1 (`wpc.run --max-lat 65`, data/results_m65); nivå 1 och 2 från körningen
utan den maskningen (data/results_mask2km). Maskningen tar bort ett
nord–sydligt fel i ljus per invånare som bara märks när en enhet spänner över
gradienten (de nordiska ländernas centroider låg 50–150 km för långt
norrut); på nivå 1 och 2 skulle den i stället ta bort hela enheter och ge
stora fel i enheter som skärs av 65° N (avsnitt 11).

En CSV-fil per nivå, i långt format (en rad per enhet och år):

  unit_id, GID_0, COUNTRY[, GID_1, NAME_1[, GID_2, NAME_2]], year, sensor,
  bm_lon, bm_lat, aa_lon, aa_lat, bm_flag, aa_flag, splice_ok,
  splice_shift_km, noise_dmsp_km, noise_viirs_km, n_cells, area_km2, max_lat

  sensor          DMSP (1992–2012) eller VIIRS (2013–2025)
  bm_*, aa_*      Barmore- och Aboufadel–Austin-centroid (grader, WGS84)
  *_flag          0 i egen enhet, 1 utanför alla enheter, 2 i en annan enhet,
                  -1 ingen centroid
  splice_ok       1 om enheten kunde skarvas (centroid i båda källorna 2013);
                  0 betyder att DMSP-åren inte är flyttade och inte är
                  jämförbara med VIIRS-åren
  splice_shift_km hur långt DMSP-åren flyttades (Barmore), 0 för VIIRS-år
  noise_*_km      enhetens typiska årshopp (median av avståndet mellan på
                  varandra följande år, Barmore) under DMSP- respektive
                  VIIRS-åren – ett mått på brus, inte ett konfidensintervall
  n_cells         antal rutnätsceller som enheten täcker i källan
  area_km2        enhetens yta enligt rutnätet
  max_lat         nordligaste latitud med ljus som används (65 på nivå 0,
                  tom på nivå 1–2 där inget är maskat)

Filerna skrivs gzip-komprimerade (`.csv.gz`; pandas, R och de flesta
verktyg läser dem direkt). En GeoPackage per nivå
(`wpc_v2_gadm41_adm<N>.gpkg`) har samma rader och kolumner i lagret
`centroids`, med Barmore-centroiden som punkt (WGS84, EPSG:4326);
Aboufadel–Austin finns som kolumnerna aa_lon och aa_lat. Enheter och år utan
centroid har tom geometri. GeoPackage ersätter v1:s shapefiler (en fil per
år, nivå och metod).

Dessutom skrivs metadata.json (källor, versioner, valideringens lägesfel
per nivå) och en kort README.md i utdatamappen.

CSV-filerna kräver bara pandas och numpy; GeoPackage kräver geopandas
(`--formats csv` hoppar över den).

Exempel
  python tools/build_release.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_series import read_results  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc import __version__  # noqa: E402
from wpc.progress import Progress, log  # noqa: E402

R_KM = 6371.0
ID_COLS = {0: ["GID_0", "COUNTRY"], 1: ["GID_0", "COUNTRY", "GID_1", "NAME_1"],
           2: ["GID_0", "COUNTRY", "GID_1", "NAME_1", "GID_2", "NAME_2"]}


def gc_km(lon1, lat1, lon2, lat2):
    l1, p1, l2, p2 = map(np.radians, (lon1, lat1, lon2, lat2))
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin((l2 - l1) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def build_level(src_dir: Path, level: int, year0: int, max_lat: float | None = None) -> pd.DataFrame:
    files = sorted(src_dir.glob("*.csv"))
    d = pd.concat([read_results(f) for f in files], ignore_index=True)
    d["sensor"] = np.where(d.year >= year0, "VIIRS", "DMSP")
    for c in ("d_east_km_bm", "d_north_km_bm"):
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0.0)
    d["splice_shift_km"] = np.hypot(d.d_east_km_bm, d.d_north_km_bm)
    d = d.sort_values(["unit_id", "year"]).reset_index(drop=True)
    g = d.groupby("unit_id")
    prev_lon, prev_lat, prev_year = g.bm_lon.shift(), g.bm_lat.shift(), g.year.shift()
    jump = gc_km(prev_lon, prev_lat, d.bm_lon, d.bm_lat)
    same = (d.year - prev_year == 1) & (g.sensor.shift() == d.sensor)
    d["_jump"] = np.where(same, jump, np.nan)
    noise = d.groupby(["unit_id", "sensor"])["_jump"].median().unstack()
    d["noise_dmsp_km"] = d.unit_id.map(noise.get("DMSP", pd.Series(dtype=float)))
    d["noise_viirs_km"] = d.unit_id.map(noise.get("VIIRS", pd.Series(dtype=float)))
    d["max_lat"] = max_lat if max_lat is not None else np.nan
    cols = (["unit_id"] + [c for c in ID_COLS[level] if c in d.columns and c != f"GID_{level}"]
            + ["year", "sensor", "bm_lon", "bm_lat", "aa_lon", "aa_lat", "bm_flag", "aa_flag",
               "splice_ok", "splice_shift_km", "noise_dmsp_km", "noise_viirs_km", "n_cells", "area_km2", "max_lat"])
    return d[[c for c in cols if c in d.columns]]


GPKG_LAYER = "centroids"


def write_csv(d: pd.DataFrame, path: Path) -> None:
    # mtime=0 gör filen bitidentisk mellan körningar med samma innehåll.
    d.to_csv(path, index=False, float_format="%.6f", compression={"method": "gzip", "mtime": 0})


def write_gpkg(d: pd.DataFrame, path: Path) -> None:
    import geopandas as gpd  # lazy: krävs bara för GeoPackage

    if path.exists():
        path.unlink()
    ok = d.bm_lon.notna() & d.bm_lat.notna()
    geom = gpd.GeoSeries.from_xy(d.bm_lon, d.bm_lat, crs="EPSG:4326")
    geom[~ok] = None
    gpd.GeoDataFrame(d, geometry=geom, crs="EPSG:4326").to_file(
        path, layer=GPKG_LAYER, driver="GPKG", engine="pyogrio")


def main(argv=None):
    p = argparse.ArgumentParser(description="Bygg den publicerbara datamängden.")
    p.add_argument("--source", default="spliced-back")
    p.add_argument("--results", type=Path, default=Path("data/results_mask2km"),
                   help="Körning för nivå 1 och 2 (utan maskning norr om 65° N).")
    p.add_argument("--results-adm0", type=Path, default=Path("data/results_m65"),
                   help="Körning för nivå 0 (wpc.run --max-lat 65).")
    p.add_argument("--adm0-max-lat", type=float, default=65.0)
    p.add_argument("--boundary-set", default="gadm41")
    p.add_argument("--validation", type=Path, default=Path("docs/validation"))
    p.add_argument("--out", type=Path, default=None, help="Standard: data/release/v<version>.")
    p.add_argument("--formats", nargs="+", choices=["csv", "gpkg"], default=["csv", "gpkg"])
    args = p.parse_args(argv)
    out = args.out or Path("data/release") / f"v{__version__}"
    out.mkdir(parents=True, exist_ok=True)

    meta = {"title": "World Population Centroids v2",
            "version": __version__, "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "git_commit": git_commit(), "source_series": args.source, "boundaries": "GADM 4.1",
            "method": "docs/METHODS.md, docs/HARMONIZATION.md", "files": {}}
    for level in (0, 1, 2):
        root, max_lat = (args.results_adm0, args.adm0_max_lat) if level == 0 else (args.results, None)
        src = root / args.boundary_set / args.source / f"adm{level}"
        splice = json.loads((src / "splice.json").read_text())
        year0 = splice["splice_year"]
        with Progress(f"adm{level}: läser {src}") as prog:
            d = build_level(src, level, year0, max_lat)
            # Samma precision i CSV och GeoPackage (6 decimaler ≈ 0,1 m).
            fl = d.select_dtypes("float").columns
            d[fl] = d[fl].round(6)
            f = out / f"wpc_v2_gadm41_adm{level}.csv.gz"
            prog.set_text(f"adm{level}: skriver {f.name}")
            write_csv(d, f)
            stale = f.with_suffix("")  # okomprimerad CSV från en äldre körning
            if stale.exists():
                stale.unlink()
            prog.set_text(f"adm{level}: {len(d)} rader, {d.unit_id.nunique()} enheter -> {f.name}")
        files = {f.name: sha256(f)}
        if "gpkg" in args.formats:
            g = out / f"wpc_v2_gadm41_adm{level}.gpkg"
            with Progress(f"adm{level}: skriver {g.name}") as prog:
                write_gpkg(d, g)
            files[g.name] = sha256(g)
        meta["files"][f"adm{level}"] = {
            "level": level, "results": str(root), "max_lat": max_lat, "rows": int(len(d)), "units": int(d.unit_id.nunique()),
            "years": [int(d.year.min()), int(d.year.max())], "splice": splice,
            "units_not_spliced": int(d[d.year == year0 - 1].splice_ok.eq(0).sum()),
            "noise_km_median": {"DMSP": float(d.drop_duplicates("unit_id").noise_dmsp_km.median()),
                                "VIIRS": float(d.drop_duplicates("unit_id").noise_viirs_km.median())},
            "sha256": files}
    val = {}
    for ds in ("census", "geostat", "inegi"):
        vp = args.validation / f"{ds}_position.csv"
        if vp.exists():
            v = pd.read_csv(vp)
            v = v[v.series == "release"]
            val[ds] = {f"adm{int(r.level)}_{int(r.year)}": {"n": int(r.n), "bm_median_km": round(r.bm_median, 3),
                                                           "bm_p90_km": round(r.bm_p90, 3)}
                       for r in v.itertuples()}
    meta["validation_position_error"] = val
    (out / "metadata.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n")
    doc = __doc__.split("En CSV-fil per nivå", 1)[1].split("Dessutom skrivs", 1)[0]
    vtab = "\n".join(f"| {ds} | {k} | {v['n']} | {v['bm_median_km']} | {v['bm_p90_km']} |"
                     for ds, d in val.items() for k, v in d.items())
    (out / "README.md").write_text(
        f"# World Population Centroids v2 ({__version__})\n\n"
        "Befolkningscentroider ur nattljus för GADM 4.1 nivå 0–2, 1992–2025. DMSP-OLS 1992–2012 "
        "(samregistrerad, IRQR-kalibrerad) är flyttad per enhet till VIIRS VNL v2, som används 2013–2025. "
        "Metod: docs/METHODS.md och docs/HARMONIZATION.md i "
        "https://github.com/ohall2/world-population-centroids (gren v2).\n\n"
        "## Filer\n\nEn CSV-fil per nivå" + doc +
        "\n## Lägesfel mot folkräkningar (Barmore, km)\n\n"
        "census = USA (block groups), geostat = Europa (1 km-rutor), inegi = Mexiko (orter). "
        "Nivå 0 mäts med maskning norr om 65° N, nivå 1–2 utan.\n\n"
        "| Facit | Nivå och år | Enheter | Median | 90 % |\n| --- | --- | --- | --- | --- |\n" + vtab + "\n\n"
        "## Begränsningar\n\n"
        "- Förändring över tio år kan inte mätas på subnationell nivå i USA, Europa och Mexiko: felet i "
        "förändringen är större än befolkningens faktiska förflyttning.\n"
        "- Mexikos nationella centroid ligger omkring 130 km fel: ljus per invånare följer inkomst, "
        "som är högre i norr. Liknande fel kan finnas i andra länder med stora regionala "
        "inkomstskillnader.\n"
        "- Nivå 0 använder inget ljus norr om 65° N (som v1); nivå 1–2 använder allt ljus. "
        "Enheter helt norr om 65° N, som Svalbard, saknar därför nationell centroid.\n"
        "- Enheter med splice_ok = 0 har DMSP-år som inte är jämförbara med VIIRS-åren.\n"
        "- Enhets-ID och namn följer GADM 4.1, inklusive dess egna koder som \"?\" för okända enheter.\n")
    log(f"metadata och README -> {out}")
    # Slutlig validering, förändringstest och README på engelska (ersätter texten ovan).
    from write_release_docs import main as write_docs  # noqa: E402
    write_docs(["--release", str(out)])


if __name__ == "__main__":
    main()
