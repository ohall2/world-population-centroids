#!/usr/bin/env python3
"""
Ladda upp datamängden till Zenodo som ett utkast (publiceras inte).

Skapar en deposition med metadata (titel, beskrivning ur README, skapare med
ORCID, licens CC BY-NC 4.0, nyckelord, länk till koden), reserverar en DOI och
laddar upp filerna i utgåvans mapp. Publiceringen görs sedan för hand på
zenodo.org efter granskning, eftersom en publicerad post inte kan tas bort.

Nyckeln (personal access token med behörigheterna deposit:write och
deposit:actions) läses från miljövariabeln ZENODO_TOKEN eller från filen
~/.zenodo_token. Den skrivs aldrig ut och sparas inte i repot.

Depositionens nummer och den reserverade DOI:n sparas i <utgåva>/zenodo.json.
Körs skriptet igen fortsätter det med samma deposition och laddar bara upp
filer som saknas eller har ändrats (jämförs med md5).

Exempel
  python tools/zenodo_upload.py                       # data/release/v2.0.0
  python tools/zenodo_upload.py --sandbox             # provkörning mot sandbox.zenodo.org
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wpc import __version__  # noqa: E402
from wpc.progress import Progress, log  # noqa: E402

FILES = ["wpc_v2_gadm41_adm0.csv.gz", "wpc_v2_gadm41_adm1.csv.gz", "wpc_v2_gadm41_adm2.csv.gz",
         "wpc_v2_gadm41_adm0.gpkg", "wpc_v2_gadm41_adm1.gpkg", "wpc_v2_gadm41_adm2.gpkg",
         "README.md", "metadata.json"]
CHUNK = 1 << 20


def token() -> str:
    t = os.environ.get("ZENODO_TOKEN")
    f = Path.home() / ".zenodo_token"
    if not t and f.exists():
        t = f.read_text().strip()
    if not t:
        sys.exit("Ingen nyckel: sätt ZENODO_TOKEN eller lägg nyckeln i ~/.zenodo_token.")
    return t


def api(method: str, url: str, tok: str, body=None, raw=None, length=None, ctype="application/json"):
    data = json.dumps(body).encode() if body is not None else raw
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {tok}")
    if data is not None:
        req.add_header("Content-Type", ctype)
    if length is not None:
        req.add_header("Content-Length", str(length))
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            txt = r.read().decode()
            return json.loads(txt) if txt else {}
    except urllib.error.HTTPError as e:
        msg = e.read().decode(errors="replace")[:500]
        sys.exit(f"Zenodo svarade {e.code} på {method} {url.split('?')[0]}: {msg}")


class Reader:
    """Fil som läses i bitar och flyttar förloppsraden."""

    def __init__(self, path: Path, prog: Progress):
        self.f, self.prog = open(path, "rb"), prog

    def read(self, n=-1):
        b = self.f.read(CHUNK if n is None or n < 0 else min(n, CHUNK))
        self.prog.advance(len(b))
        return b


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def metadata(rel: Path) -> dict:
    readme = (rel / "README.md").read_text()
    intro = readme.split("\n## ", 1)[0].split("\n", 2)[-1].strip()
    good = readme.split("## What the data are good for", 1)[1].split("\n## ", 1)[0].strip()
    items = "".join(f"<li>{' '.join(x.strip().lstrip('- ').replace('**', '').split())}</li>" for x in good.split("\n- ") if x.strip())
    desc = (f"<p>{' '.join(intro.split())}</p><p>What the data are good for:</p><ul>{items}</ul>"
            "<p>Files: one gzip-compressed CSV and one GeoPackage per GADM level (0 countries, 1 regions, "
            "2 municipalities), one row per unit and year, with Barmore and Aboufadel–Austin centroids, flags and "
            "noise levels. README.md describes the columns, the validation against censuses in the United States, "
            "Europe, Mexico and Brazil, and the limitations; metadata.json holds checksums and the validation "
            "tables.</p>")
    return {"metadata": {
        "upload_type": "dataset",
        "title": "World Population Centroids v2: population centroids of the world's administrative units "
                 "from nighttime lights, 1992–2025",
        "publication_date": dt.date.today().isoformat(),
        "version": __version__,
        "creators": [{"name": "Hall, Ola", "affiliation": "Lund University", "orcid": "0000-0002-9231-4028"}],
        "description": desc,
        "access_right": "open",
        "license": "cc-by-nc-4.0",
        "keywords": ["population centroids", "centre of population", "nighttime lights", "DMSP-OLS", "VIIRS",
                     "GADM"],
        "related_identifiers": [{"identifier": "https://github.com/ohall2/world-population-centroids",
                                 "relation": "isSupplementedBy", "resource_type": "software"},
                                {"identifier": "10.1038/s41597-019-0250-z", "relation": "isNewVersionOf",
                                 "resource_type": "publication-article"}],
        "prereserve_doi": True,
    }}


def main(argv=None):
    p = argparse.ArgumentParser(description="Ladda upp datamängden till Zenodo som utkast.")
    p.add_argument("--release", type=Path, default=Path("data/release") / f"v{__version__}")
    p.add_argument("--sandbox", action="store_true", help="Använd sandbox.zenodo.org (provkörning).")
    a = p.parse_args(argv)
    base = "https://sandbox.zenodo.org/api" if a.sandbox else "https://zenodo.org/api"
    tok = token()
    missing = [f for f in FILES if not (a.release / f).exists()]
    if missing:
        sys.exit(f"Saknas i {a.release}: {', '.join(missing)}")

    state_f = a.release / ("zenodo_sandbox.json" if a.sandbox else "zenodo.json")
    state = json.loads(state_f.read_text()) if state_f.exists() else {}
    if state.get("id"):
        dep = api("GET", f"{base}/deposit/depositions/{state['id']}", tok)
        log(f"Fortsätter med deposition {dep['id']}")
    else:
        dep = api("POST", f"{base}/deposit/depositions", tok, body={})
        log(f"Skapade deposition {dep['id']}")
    dep = api("PUT", f"{base}/deposit/depositions/{dep['id']}", tok, body=metadata(a.release))
    doi = dep["metadata"].get("prereserve_doi", {}).get("doi")
    state = {"id": dep["id"], "doi": doi, "edit": dep["links"].get("html"), "sandbox": a.sandbox}
    state_f.write_text(json.dumps(state, indent=1) + "\n")
    log(f"Metadata satt. Reserverad DOI: {doi}")

    have = {f["filename"]: f["checksum"] for f in api("GET", f"{base}/deposit/depositions/{dep['id']}/files", tok)}
    bucket = dep["links"]["bucket"]
    todo = [f for f in FILES if have.get(f) != md5(a.release / f)]
    total = sum((a.release / f).stat().st_size for f in todo)
    with Progress(f"laddar upp {len(todo)} filer", total=total) as prog:
        for f in todo:
            prog.set_text(f"laddar upp {f}")
            path = a.release / f
            api("PUT", f"{bucket}/{f}", tok, raw=Reader(path, prog), length=path.stat().st_size,
                ctype="application/octet-stream")
    if not todo:
        log("Alla filer finns redan och stämmer.")
    log(f"Klart. Utkastet är inte publicerat. Granska och publicera här: {state['edit']}")


if __name__ == "__main__":
    main()
