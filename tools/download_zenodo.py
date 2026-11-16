#!/usr/bin/env python3
"""
World Population Centroids v2 – hämta filer från en Zenodo-post.

Används för fackelkatalogen MYC25Q4 (--product flares-viirs) och för CCNL v1 (Zhao m.fl. 2022), som i originalet har 30 bågsekunders
upplösning. Versionen i Earth Engine (BNU/FGS/CCNL/v1) är omsamplad till ett
rutnät på 0,008983° (1 km vid ekvatorn) och används därför inte.

Hur det fungerar
  1. Fillista, storlekar och md5-summor läses från Zenodos API.
  2. Varje fil hämtas till <namn>.part. Avbryts hämtningen fortsätter nästa
     körning där den slutade (HTTP Range). Varje läsning har en tidsgräns.
  3. När filen är komplett kontrolleras md5-summan mot Zenodo; först då får
     filen sitt riktiga namn.
  4. Bredvid varje fil skrivs en JSON-fil med post, DOI, md5, år och tid.
     Färdiga filer (fil + JSON) hoppas över.

Exempel
  python tools/download_zenodo.py --product ccnl
  python tools/download_zenodo.py --record 6644980 --out data/raw/ccnl --dry-run
"""

from __future__ import annotations

import argparse
import datetime as dt
import fnmatch
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc.progress import Progress, fmt_seconds, log  # noqa: E402

PRODUCTS = {
    "ccnl": {"record": "6644980", "pattern": "CCNL_DMSP_*_V1.tif",
             "out": "data/raw/ccnl",
             "citation": "Zhao, C., Cao, X., Chen, X., Cui, X. (2022). A consistent and "
                         "corrected nighttime light dataset (CCNL 1992–2013) from DMSP-OLS "
                         "data. Scientific Data 9, 424. doi:10.1038/s41597-022-01540-x"},
    # Fackelplatser för VIIRS-perioden (CC BY-NC 4.0).
    "flares-viirs": {"record": "23091822", "pattern": "*", "out": "data/raw/flares_viirs",
                     "citation": "Zhizhin, M., Ghosh, T., Gleason, G., Elvidge, C. D., Bazilian, M. "
                                 "(2026). VIIRS Nightfire upstream gas-flaring volumes, 2012–2025 "
                                 "(MYC25Q4), version 1.0. Zenodo. doi:10.5281/zenodo.23091822"},
}
API = "https://zenodo.org/api/records/{record}"
TIMEOUT = (30, 120)          # (anslutning, läsning) i sekunder
CHUNK = 1 << 20
MB = 1e6


def record_files(record: str, pattern: str) -> tuple[dict, list[dict]]:
    r = requests.get(API.format(record=record), timeout=TIMEOUT)
    r.raise_for_status()
    rec = r.json()
    files = []
    for f in rec["files"]:
        if not fnmatch.fnmatch(f["key"], pattern):
            continue
        algo, _, digest = f["checksum"].partition(":")
        if algo != "md5":
            raise ValueError(f"{f['key']}: oväntad checksumma {f['checksum']}")
        files.append({"key": f["key"], "size": int(f["size"]), "md5": digest,
                      "url": f["links"]["self"]})
    files.sort(key=lambda f: f["key"])
    return rec, files


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def year_of(name: str) -> int | None:
    m = re.search(r"(19|20)\d{2}", name)
    return int(m.group(0)) if m else None


def fetch(f: dict, part: Path, prog: Progress, retries: int = 8) -> None:
    """Hämta (eller fortsätt hämta) en fil till `part`."""
    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0
        if have >= f["size"]:
            return
        try:
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(f["url"], headers=headers, stream=True, timeout=TIMEOUT) as r:
                if have and r.status_code != 206:      # servern ignorerade Range
                    have = 0
                    part.unlink(missing_ok=True)
                r.raise_for_status()
                with open(part, "ab") as out:
                    for chunk in r.iter_content(CHUNK):
                        out.write(chunk)
                        have += len(chunk)
                        prog.set_done(int(have / MB))
            if part.stat().st_size >= f["size"]:
                return
        except (requests.RequestException, OSError) as exc:
            if attempt == retries:
                raise RuntimeError(f"{f['key']}: misslyckades efter {retries} försök: {exc}") from exc
            wait = min(60, 2 ** attempt)
            log(f"  {f['key']}: {str(exc)[:120]} – fortsätter om {wait} s")
            time.sleep(wait)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Hämta filer från en Zenodo-post.")
    p.add_argument("--product", choices=sorted(PRODUCTS), default=None)
    p.add_argument("--record", default=None, help="Zenodo-postens nummer.")
    p.add_argument("--pattern", default=None, help="Filnamnsmönster, t.ex. '*.tif'.")
    p.add_argument("--out", default=None)
    p.add_argument("--dry-run", action="store_true", help="Visa filerna utan att hämta.")
    args = p.parse_args(argv)

    prod = PRODUCTS.get(args.product, {})
    record = args.record or prod.get("record")
    pattern = args.pattern or prod.get("pattern", "*")
    if not record:
        sys.exit("Ange --product eller --record.")
    out_dir = Path(args.out or prod.get("out", f"data/raw/zenodo_{record}"))

    rec, files = record_files(record, pattern)
    total = sum(f["size"] for f in files)
    log(f"Zenodo {record} ({rec.get('doi', '')}): {len(files)} filer, {total / 1e9:.1f} GB")
    if args.dry_run:
        for f in files:
            print(f"  {f['key']}  {f['size'] / MB:8.1f} MB  md5 {f['md5']}")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [f for f in files
            if not ((out_dir / f["key"]).exists()
                    and (out_dir / f["key"]).with_suffix(".json").exists())]
    log(f"{len(files) - len(todo)} filer finns redan, {len(todo)} hämtas till {out_dir.resolve()}")
    t0 = time.time()
    for i, f in enumerate(todo, start=1):
        dst = out_dir / f["key"]
        part = dst.with_name(dst.name + ".part")
        have = part.stat().st_size if part.exists() else 0
        with Progress(f"[{i}/{len(todo)}] {f['key']} MB",
                      total=int(f["size"] / MB)) as prog:
            prog.set_done(int(have / MB))
            fetch(f, part, prog)
        with Progress(f"[{i}/{len(todo)}] {f['key']}: kontrollerar md5"):
            got = md5(part)
        if got != f["md5"]:
            part.unlink()
            sys.exit(f"{f['key']}: fel md5 ({got}, väntade {f['md5']}). Filen är borttagen; "
                     "kör samma kommando igen.")
        part.replace(dst)
        meta = {"source": "Zenodo", "record": record, "doi": rec.get("doi"),
                "title": rec.get("metadata", {}).get("title"),
                "license": (rec.get("metadata", {}).get("license") or {}).get("id"),
                "citation": prod.get("citation"),
                "file": f["key"], "url": f["url"], "bytes": f["size"], "md5": f["md5"],
                "year": year_of(f["key"]),
                "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        dst.with_suffix(".json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    log(f"Klart på {fmt_seconds(time.time() - t0)}.")


if __name__ == "__main__":
    main()
