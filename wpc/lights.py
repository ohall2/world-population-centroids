"""Ljuskällor för centroidberäkningen.

En ljuskälla levererar, för ett givet geografiskt område, ett fönster ur ett
globalt rutnät med ljusvärden (float32) och fönstrets affina transform.
Allt annat i pipelinen är oberoende av vilken källa som används.

Implementerat:
  DmspIrqr  DMSP-OLS v4 Stable Lights, kalibrerade med IRQR-koefficienterna i
            Hall m.fl. (2019, tabell 2), medelvärde per år, gasfacklor maskade.
  Annual    En färdig årsprodukt med en fil per år, hämtad med
            tools/download_ee.py: CCNL (Zhao m.fl. 2022) och VIIRS VNL v2.
            Värdena används som de är (negativa sätts till 0), gasfacklor
            maskas om en fackelfil anges.

make_source(name, year, ...) skapar rätt källa; available_years(name, ...)
anger vilka år som finns på disk.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio.features import rasterize
from rasterio.windows import Window

# DMSP-OLS v4: 43 201 x 16 801 pixlar à 30 bågsekunder. Pixelmittpunkterna
# ligger på hela 30"-steg från -180 till 180 och från 75 till -65, så hörnet
# ligger en halv pixel utanför. Earth Engine anger samma hörn för 32 av 35
# kompositer (skalan avrundad till 0.0083333333); för F162009, F182010 och
# F182011 anger Earth Engine hörnet (-180, 75), ett känt problem med dessa tre
# raster (docs/HARMONIZATION.md 7.1). Vi använder EOG:s rutnät för alla; den
# verkliga förskjutningen mäts och korrigeras av samregistreringen.
DMSP_WIDTH, DMSP_HEIGHT = 43201, 16801
DMSP_TRANSFORM = Affine(1 / 120, 0.0, -180.0 - 1 / 240, 0.0, -1 / 120, 75.0 + 1 / 240)


@dataclass
class Window2D:
    """Ett utsnitt ur ett globalt rutnät."""
    data: np.ndarray          # float32, rader x kolumner
    transform: Affine         # fönstrets transform

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        h, w = self.data.shape
        xmin, ymax = self.transform.c, self.transform.f
        return xmin, ymax + self.transform.e * h, xmin + self.transform.a * w, ymax


def grid_window(bounds, transform: Affine, width: int, height: int) -> Window | None:
    """Alla pixlar som rektangeln berör, begränsat till rutnätet."""
    minx, miny, maxx, maxy = bounds
    inv = ~transform
    c0, r0 = inv @ (minx, maxy)
    c1, r1 = inv @ (maxx, miny)
    col0 = max(int(np.floor(min(c0, c1))), 0)
    row0 = max(int(np.floor(min(r0, r1))), 0)
    col1 = min(int(np.ceil(max(c0, c1))), width)
    row1 = min(int(np.ceil(max(r0, r1))), height)
    if col1 <= col0 or row1 <= row0:
        return None
    return Window(col0, row0, col1 - col0, row1 - row0)


def calibrate(dn: np.ndarray, c0: float, c1: float, c2: float) -> np.ndarray:
    """IRQR-kalibrering c0 + c1*DN + c2*DN^2. Pixlar med DN = 0 förblir 0."""
    dn = dn.astype("float32")
    return np.where(dn > 0, c0 + c1 * dn + c2 * dn * dn, 0.0).astype("float32")


def read_calibration(path: Path) -> dict[tuple[str, int], tuple[float, float, float]]:
    rows = [r for r in csv.reader(l for l in path.open(encoding="utf-8") if not l.startswith("#"))]
    header, body = rows[0], rows[1:]
    idx = {k: header.index(k) for k in ("satellite", "year", "c0", "c1", "c2")}
    return {(r[idx["satellite"]], int(r[idx["year"]])):
            (float(r[idx["c0"]]), float(r[idx["c1"]]), float(r[idx["c2"]])) for r in body}


class DmspIrqr:
    """DMSP-OLS v4 Stable Lights, IRQR-kalibrerade, ett år.

    - Varje satellit kalibreras med c0 + c1*DN + c2*DN^2 (tabell 2 i Hall m.fl. 2019).
      Pixlar med DN = 0 förblir 0: utan kravet skulle varje mörk pixel få värdet
      c0 och få vikt i centroiden.
    - År med två satelliter: medelvärdet av de kalibrerade värdena.
    - Pixlar som berörs av en fackelpolygon sätts till 0.
    - Med `offsets` (satellit, år) -> (dx, dy) flyttas varje komposit efter
      kalibreringen så att den ligger rätt mot VIIRS (wpc/coregister.py,
      tools/coregister_dmsp.py). Källan heter då "dmsp-irqr-reg".
    """

    name = "dmsp-irqr"

    def __init__(self, year: int, raw_dir: Path, calibration: Path, flares: Path | None,
                 offsets: dict | None = None):
        coeffs = read_calibration(calibration)
        self.year = year
        self.satellites = sorted(sat for (sat, y) in coeffs if y == year)
        if not self.satellites:
            raise ValueError(f"Inga kalibreringskoefficienter för {year}")
        self.coeffs = [coeffs[(s, year)] for s in self.satellites]
        self.paths = [raw_dir / f"DMSP_{s}{year}_stable_lights.tif" for s in self.satellites]
        for p in self.paths:
            if not p.exists():
                raise FileNotFoundError(p)
        self._srcs = [rasterio.open(p) for p in self.paths]
        for s in self._srcs:
            if (s.width, s.height) != (DMSP_WIDTH, DMSP_HEIGHT):
                raise ValueError(f"{s.name}: {s.width}x{s.height}, väntade {DMSP_WIDTH}x{DMSP_HEIGHT}")
        self.flares = gpd.read_file(flares, engine="pyogrio") if flares else None
        self._flare_index = self.flares.sindex if self.flares is not None else None
        self.offsets = None
        if offsets is not None:
            missing = [s for s in self.satellites if (s, year) not in offsets]
            if missing:
                raise ValueError(f"Förskjutning saknas för {missing} {year}")
            self.offsets = [offsets[(s, year)] for s in self.satellites]
            self.name = "dmsp-irqr-reg"

    def describe(self) -> dict:
        d = {"source": self.name, "year": self.year, "satellites": self.satellites,
             "coefficients": dict(zip(self.satellites, self.coeffs)),
             "files": [str(p) for p in self.paths],
             "grid": "EOG DMSP v4, 30 arcsec, hörn (-180.0041667, 75.0041667)"}
        if self.offsets:
            d["offsets_px"] = dict(zip(self.satellites, [list(o) for o in self.offsets]))
        return d

    # Rutnätet som källan levererar fönster ur.
    transform, width, height = DMSP_TRANSFORM, DMSP_WIDTH, DMSP_HEIGHT

    def read(self, bounds) -> Window2D | None:
        win = grid_window(bounds, DMSP_TRANSFORM, DMSP_WIDTH, DMSP_HEIGHT)
        if win is None:
            return None
        return self.read_window(win)

    def read_window(self, win: Window) -> Window2D:
        acc = None
        offs = self.offsets or [None] * len(self._srcs)
        for src, (c0, c1, c2), off in zip(self._srcs, self.coeffs, offs):
            if off is None:
                cal = calibrate(src.read(1, window=win), c0, c1, c2)
            else:
                from .coregister import shift_bilinear
                h = int(np.ceil(max(abs(off[0]), abs(off[1])))) + 1
                big = Window(win.col_off - h, win.row_off - h, win.width + 2 * h, win.height + 2 * h)
                raw = src.read(1, window=big, boundless=True, fill_value=0)
                cal = shift_bilinear(calibrate(raw, c0, c1, c2), off[0], off[1])
                cal = cal[h:h + int(win.height), h:h + int(win.width)]
            acc = cal if acc is None else acc + cal
        data = acc / len(self._srcs)
        transform = DMSP_TRANSFORM @ Affine.translation(win.col_off, win.row_off)
        out = Window2D(data=data, transform=transform)
        self._mask_flares(out)
        return out

    def _mask_flares(self, w: Window2D) -> None:
        mask_flares(w, self.flares, self._flare_index)

    def close(self) -> None:
        for s in self._srcs:
            s.close()


def mask_flares(w: Window2D, flares, index) -> None:
    """Sätt alla pixlar som en fackelpolygon berör till 0 (all_touched)."""
    if index is None:
        return
    hits = index.query(_box(w.bounds))
    if len(hits) == 0:
        return
    mask = rasterize(((g, 1) for g in flares.geometry.iloc[hits]),
                     out_shape=w.data.shape, transform=w.transform,
                     fill=0, all_touched=True, dtype="uint8")
    w.data[mask.astype(bool)] = 0.0


# Årsprodukter: filnamn (utan .tif/.json). Varje fil har en JSON-fil bredvid
# med fältet "year", skriven av tools/download_zenodo.py eller download_ee.py.
ANNUAL = {
    "ccnl": {"glob": "CCNL_DMSP_*_V1",
             "note": "CCNL v1, 30 bågsekunder (Zhao m.fl. 2022, doi:10.5281/zenodo.6644980)"},
    "viirs": {"glob": "VIIRS_*_average_masked",
              "note": "VIIRS VNL v2 årskompositer, average_masked (nW/cm2/sr)"},
    "viirs-dmsp": {"glob": "VDMSP_*",
                   "note": "VIIRS VNL v2 harmoniserad till DMSP-IRQR (tools/harmonize_viirs.py)"},
    "viirs-dmsp-reg": {"glob": "VDMSP_*",
                       "note": "VIIRS VNL v2 harmoniserad till samregistrerad DMSP-IRQR"},
    # Publicerade harmoniserade serier, för jämförelse (docs/HARMONIZATION.md).
    "li2020-dmsp": {"glob": "LI2020_*_calDMSP_b1",
                    "note": "Li m.fl. 2020 v7, kalibrerad DMSP 1992–2013 (sat-io Harmonized_NTL/dmsp)"},
    "li2020-viirs": {"glob": "LI2020_*_simVIIRS_b1",
                     "note": "Li m.fl. 2020 v7, DMSP-lik VIIRS 2013–2021 (sat-io Harmonized_NTL/viirs)"},
    "chiovelli-cln": {"glob": "DMSP*_cln",
                      "note": "Chiovelli m.fl. 2026, clean_dmsp (zenodo 20558766)"},
    "chiovelli-tc": {"glob": "DMSP*_tcfix",
                     "note": "Chiovelli m.fl. 2026, topcode_fix (zenodo 20558766)"},
    "chiovelli-bl": {"glob": "DMSP*_blfix",
                     "note": "Chiovelli m.fl. 2026, bloom_fix (zenodo 20558766)"},
    "chiovelli-bltc": {"glob": "DMSP*_bltcfix",
                       "note": "Chiovelli m.fl. 2026, bloomtopcode_fix (zenodo 20558766)"},
    # Jämförelsekällor för analys B (docs/ANALYSIS_PLAN.md §6). Totaler per cell;
    # räknas om till täthet i DensitySource så centroidkoden ger en viktad
    # centroid (befolkning respektive bebyggd yta) utan att ändras.
    "ghspop": {"glob": "GHS_POP_E*_GLOBE_R2023A_4326_30ss_V1_0",
               "note": "GHS-POP R2023A, 30\", befolkning per cell (JRC GHSL, CC BY 4.0)"},
    "ghsbuilts": {"glob": "GHS_BUILT_S_E*_GLOBE_R2023A_4326_30ss_V1_0",
                  "note": "GHS-BUILT-S R2023A, 30\", bebyggd yta (m²) per cell (JRC GHSL, CC BY 4.0)"},
    "worldpop": {"glob": "ppp_*_1km_Aggregated",
                 "note": "WorldPop Global unconstrained, 1 km, befolkning per cell "
                         "(beror på nattljus; CC BY 4.0)"},
}

# Källor vars cellvärde är en total per cell (befolkning, bebyggd yta) och som
# därför räknas om till täthet (värde/cellyta) innan centroidkoden används.
DENSITY_SOURCES = {"ghspop", "ghsbuilts", "worldpop"}


def snap_transform(t: Affine, tol: float = 1e-8) -> Affine:
    """Rätta avrundade 30"/15"-rutnät till exakta värden.

    CCNL-filerna anger skalan 0.0083333333 och hörnet -180.00416666665, dvs.
    1/120 grad avrundat till tio decimaler (fel ca 0,2 m över hela bredden).
    Ligger skalan och hörnet inom `tol` från ett exakt rutnät med 1/120 eller
    1/240 grad används det exakta, så att rutnätet blir identiskt med DMSP:s.
    """
    for step in (1 / 120, 1 / 240):
        if abs(t.a - step) < tol and abs(t.e + step) < tol and t.b == 0 and t.d == 0:
            c = round(t.c / (step / 2)) * (step / 2)
            f = round(t.f / (step / 2)) * (step / 2)
            if abs(c - t.c) < 100 * tol and abs(f - t.f) < 100 * tol:
                return Affine(step, 0.0, c, 0.0, -step, f)
    return t


def annual_files(name: str, raw_dir: Path) -> dict[int, Path]:
    """År -> GeoTIFF, ur metadatafilerna som download_ee.py skriver."""
    spec = ANNUAL[name]
    out: dict[int, Path] = {}
    for js in sorted(raw_dir.glob(spec["glob"] + ".json")):
        meta = json.loads(js.read_text())
        tif = js.with_suffix(".tif")
        if not tif.exists():
            continue
        year = int(meta["year"])
        if year in out:
            raise ValueError(f"Två filer för {name} {year}: {out[year].name} och {tif.name}")
        out[year] = tif
    return out


class Annual:
    """En årsprodukt med en fil per år, i filens eget rutnät."""

    def __init__(self, name: str, year: int, raw_dir: Path, flares: Path | None):
        files = annual_files(name, raw_dir)
        if year not in files:
            raise FileNotFoundError(f"Ingen {name}-fil för {year} i {raw_dir}")
        self.name, self.year, self.path = name, year, files[year]
        self._src = rasterio.open(self.path)
        if self._src.crs is None or not self._src.crs.is_geographic:
            raise ValueError(f"{self.path.name}: väntade geografiska koordinater, fick {self._src.crs}")
        self.nodata = self._src.nodata
        self.transform = snap_transform(self._src.transform)
        self.width, self.height = self._src.width, self._src.height
        self.flares = gpd.read_file(flares, engine="pyogrio") if flares else None
        self._flare_index = self.flares.sindex if self.flares is not None else None

    def describe(self) -> dict:
        return {"source": self.name, "year": self.year, "file": str(self.path),
                "product": ANNUAL[self.name]["note"],
                "nodata": self.nodata, "dtype": self._src.dtypes[0],
                "grid": {"width": self.width, "height": self.height,
                         "transform": list(self.transform)[:6]},
                "flares": "maskade" if self.flares is not None else "ej maskade"}

    def read(self, bounds) -> Window2D | None:
        win = grid_window(bounds, self.transform, self.width, self.height)
        return None if win is None else self.read_window(win)

    def read_window(self, win: Window) -> Window2D:
        raw = self._src.read(1, window=win)
        data = raw.astype("float32")
        bad = ~np.isfinite(data) | (data < 0)
        if self.nodata is not None:
            bad |= raw == self.nodata
        data[bad] = 0.0
        out = Window2D(data=data,
                       transform=self.transform @ Affine.translation(win.col_off, win.row_off))
        mask_flares(out, self.flares, self._flare_index)
        return out

    def close(self) -> None:
        self._src.close()


class DensitySource(Annual):
    """Årsprodukt vars cellvärde är en total per cell (GHS-POP:s befolkning,
    GHS-BUILT-S:s bebyggda yta), omräknad till täthet per km².

    Centroidkoden viktar varje pixel med värde x andel x cellyta. För en total C
    per cell blir rätt vikt C x andel, inte C x andel x yta. Genom att leverera
    tätheten C/yta blir vikten (C/yta) x andel x yta = C x andel, exakt den
    totalviktade centroiden, utan att centroids.py ändras (planen §6,
    'samma centroidkod', preciserad under Avvikelser).
    """

    def read_window(self, win: Window) -> Window2D:
        from .centroids import cell_area_km2
        w = super().read_window(win)            # totaler, nodata/negativa -> 0
        t = w.transform
        lat = t.f + t.e * (np.arange(w.data.shape[0]) + 0.5)
        area = cell_area_km2(lat, -t.e, t.a)    # km² per rad
        w.data = (w.data / area[:, None]).astype("float32")
        return w

    def describe(self) -> dict:
        return {**super().describe(), "values": "total per cell omräknad till täthet per km²"}


SOURCES = ["dmsp-irqr", "dmsp-irqr-reg", *ANNUAL]
OFFSETS_FILE = Path("config/dmsp_offsets.csv")


def read_offsets(path: Path = OFFSETS_FILE) -> dict:
    """(satellit, år) -> (dx, dy) i pixlar, från tools/coregister_dmsp.py."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            out[(row["satellite"], int(row["year"]))] = (float(row["dx"]), float(row["dy"]))
    return out


class LatitudeMask:
    """Sätter ljuset till 0 i celler vars mittpunkt ligger norr om `max_lat`.

    Som i v1 (Hall m.fl. 2019): norr om omkring 65° N är nattljuset
    opålitligt (sena solnedgångar, norrsken, snö) och ljus per invånare
    mycket högre än söderut, vilket drar t.ex. Sveriges och Norges nationella
    centroider 100–150 km norrut (docs/HARMONIZATION.md). Enheter helt norr om
    gränsen får ingen centroid.
    """

    def __init__(self, source, max_lat: float):
        self._s, self.max_lat = source, max_lat
        self.transform, self.width, self.height = source.transform, source.width, source.height
        self.name = source.name

    def describe(self) -> dict:
        return {**self._s.describe(), "max_lat": self.max_lat}

    def read(self, bounds):
        win = grid_window(bounds, self.transform, self.width, self.height)
        return None if win is None else self.read_window(win)

    def read_window(self, win: Window) -> Window2D:
        w = self._s.read_window(win)
        t = w.transform
        lat = t.f + t.e * (np.arange(w.data.shape[0]) + 0.5)
        w.data[lat > self.max_lat, :] = 0.0
        return w

    def close(self) -> None:
        self._s.close()


def make_source(name: str, year: int, raw_dir: Path, calibration: Path | None,
                flares: Path | None, max_lat: float | None = None):
    if name == "dmsp-irqr":
        src = DmspIrqr(year, raw_dir, calibration, flares)
    elif name == "dmsp-irqr-reg":
        src = DmspIrqr(year, raw_dir, calibration, flares, offsets=read_offsets())
    elif name in ANNUAL:
        cls = DensitySource if name in DENSITY_SOURCES else Annual
        src = cls(name, year, raw_dir, flares)
    else:
        raise ValueError(f"Okänd ljuskälla: {name}")
    return LatitudeMask(src, max_lat) if max_lat is not None else src


def available_years(name: str, raw_dir: Path, calibration: Path | None) -> list[int]:
    if name in ("dmsp-irqr", "dmsp-irqr-reg"):
        return sorted({y for (_, y) in read_calibration(calibration)})
    return sorted(annual_files(name, raw_dir))


def _box(bounds):
    from shapely.geometry import box
    return box(*bounds)
