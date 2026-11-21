#!/usr/bin/env python3
"""Hämta Japans folkräknings-rutnätsstatistik (e-Stat statistiska GIS), 500 m-mesh.

Bekräftande facit för analys C (docs/ANALYSIS_PLAN.md §7). 500 m-mesh, total
befolkning: 2015 = statsId T000847, 2020 = T001101 (e-Stats definitionsfiler).
Hämtas per prefektur (47 filer per år) från den statistiska GIS-endpointen; varje
prefektur blir en zip med en cp932-CSV. En JSON-manifest med SHA-256 skrivs per år.

Före facit certifieras den nationella totalen mot folkräkningen (2015 = 127 094 745,
2020 = 126 146 099). Stämmer den inte avbryts körningen.

Källa/licens: Statistics Bureau of Japan / e-Stat; Government of Japan Standard
Terms of Use v2.0 (förenlig med CC BY 4.0). Se docs/DATA.md.

Exempel
  python tools/download_estat_jp.py --years 2015 2020
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path

from tools.validation_centroids import read_estat_jp

STATS = {2015: "T000847", 2020: "T001101"}          # 500 m-mesh, total befolkning
EXPECTED_TOTAL = {2015: 127_094_745, 2020: 126_146_099}
ENDPOINT = "https://www.e-stat.go.jp/gis/statmap-search/data"
UA = {"User-Agent": "Mozilla/5.0 (Macintosh) AppleWebKit/537.36 Chrome/120 Safari/537.36"}


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def download_year(year: int, out: Path) -> None:
    stats_id = STATS[year]
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"source": "Statistics Bureau of Japan / e-Stat",
                "product": f"Population Census grid square statistics, 500 m mesh, total population ({year})",
                "statsId": stats_id, "license": "Government of Japan Standard Terms of Use v2.0 (CC BY 4.0-compatible)",
                "endpoint": ENDPOINT, "year": year, "prefectures": {}}
    for p in range(1, 48):
        code = f"{p:02d}"
        fn = out / f"tbl{stats_id}_{code}.zip"
        if fn.exists():
            manifest["prefectures"][code] = {"file": fn.name, "sha256": sha256_bytes(fn.read_bytes())}
            continue
        url = f"{ENDPOINT}?statsId={stats_id}&code={code}&downloadType=2&aggregateUnit=S"
        for attempt in range(3):
            try:
                data = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120).read()
                break
            except Exception as e:
                if attempt == 2:
                    raise
                time.sleep(2)
        if data[:2] != b"PK":
            raise RuntimeError(f"{year} prefektur {code}: svar är inte en zip ({len(data)} byte)")
        fn.write_bytes(data)
        manifest["prefectures"][code] = {"file": fn.name, "sha256": sha256_bytes(data)}
        print(f"{year} {code}: {len(data)/1e6:.2f} MB")
    (out / f"manifest_{stats_id}.json").write_text(json.dumps(manifest, indent=1))


def certify(year: int, root: Path) -> float:
    _df, info = read_estat_jp(root, year)
    total = info["pop_total"]
    exp = EXPECTED_TOTAL[year]
    ok = abs(total - exp) <= 1
    print(f"{year}: nationell total {total:,.0f} (väntat {exp:,.0f}) "
          f"{'OK' if ok else 'STÄMMER INTE'}; mesh_len={info['mesh_len']}, "
          f"rutor={info['meshes']:,}, icke-numeriska={info['meshes_nonnumeric_pop']}")
    return total


def main(argv=None):
    p = argparse.ArgumentParser(description="Hämta Japans e-Stat 500 m-mesh (folkräkning).")
    p.add_argument("--years", type=int, nargs="+", default=[2015, 2020])
    p.add_argument("--out", type=Path, default=Path("data/raw/validation/estat_jp"))
    a = p.parse_args(argv)
    for year in a.years:
        download_year(year, a.out / str(year))
    print("--- certifiering mot nationell total ---")
    bad = []
    for year in a.years:
        total = certify(year, a.out)
        if abs(total - EXPECTED_TOTAL[year]) > 1:
            bad.append(year)
    if bad:
        raise SystemExit(f"Nationell total stämmer inte för {bad}; avbryter (facit räknas inte).")
    print("Alla år certifierade mot folkräkningens nationella total.")


if __name__ == "__main__":
    main()
