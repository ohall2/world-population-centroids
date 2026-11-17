#!/usr/bin/env python3
"""
World Population Centroids v2 – packa upp Chiovelli m.fl. (2026) och skriv
metadatafiler.

Zenodo-posten (doi:10.5281/zenodo.20558766, hämtad med
tools/download_zenodo.py --product chiovelli) har fyra zip-filer med en
GeoTIFF per år 1992–2025: clean_dmsp, topcode_fix, bloom_fix och
bloomtopcode_fix. Skriptet packar upp dem (filer som redan finns hoppas
över) och skriver en JSON-fil bredvid varje GeoTIFF med år, variant och
zip-filens md5, så att pipelinen hittar filerna (wpc.lights.annual_files).

Kräver bara Pythons standardbibliotek.

Exempel
  python tools/prepare_chiovelli.py
"""

from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--dir", type=Path, default=Path("data/raw/chiovelli"))
    args = p.parse_args(argv)
    for z in sorted(args.dir.glob("*.zip")):
        meta = json.loads(z.with_suffix(".json").read_text()) if z.with_suffix(".json").exists() else {}
        with zipfile.ZipFile(z) as zf:
            names = [n for n in zf.namelist() if n.endswith(".tif")]
            for n in names:
                out = args.dir / n
                if not out.exists():
                    zf.extract(n, args.dir)
                year = int(re.search(r"DMSP(\d{4})_", Path(n).name).group(1))
                out.with_suffix(".json").write_text(json.dumps({
                    "source": "Chiovelli, Michalopoulos, Papaioannou, Regan (2026), zenodo 20558766",
                    "variant": z.stem, "year": year, "zip": z.name,
                    "zip_md5": meta.get("md5")}, indent=2) + "\n")
        print(f"{z.name}: {len(names)} filer")


if __name__ == "__main__":
    main()
