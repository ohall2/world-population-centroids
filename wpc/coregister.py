"""
Samregistrering av DMSP-kompositer (subpixel-förskjutning).

DMSP-kompositernas geografiska läge varierar med omkring en halv pixel
mellan år och satelliter (docs/HARMONIZATION.md, avsnitt 5). Här uppskattas
en förskjutning (dx, dy) i pixlar per komposit mot en referensbild och
kompositen flyttas så att den ligger rätt.

Referens: den harmoniserade VIIRS-bilden 2013 (VDMSP_2013.tif) – VIIRS har
god geolokalisering (storleksordningen 100 m), och bilden ligger på DMSP:s
rutnät och är utsmetad som DMSP. Ljuset förändras mellan åren, men en
korskorrelation påverkas av förändringen bara om den är riktad, vilket är
osannolikt i genomsnitt över många rutor.

Uppskattning per ruta (TILE × TILE pixlar):
  - båda bilderna transformeras med log(1 + x) (minskar mättade
    stadskärnors vikt), medelvärdet dras ifrån och ett Hann-fönster läggs på;
  - korskorrelationen beräknas med FFT; toppen inom ±MAX_SHIFT pixlar söks
    och läget förfinas med en parabel genom toppen och grannarna i varje led;
  - rutan används om båda bilderna har minst MIN_LIT tända pixlar.
Förskjutningen för kompositen är den viktade medianen av rutornas värden
(vikt = antal tända pixlar). Spridningen (viktad MAD) redovisas.

Teckenkonvention: (dx, dy) är sådan att referens(x, y) ≈ komposit(x + dx,
y + dy), med x åt öster (kolumner) och y åt söder (rader). Den rättade
kompositen är komposit'(x, y) = komposit(x + dx, y + dy), dvs. den samplas
dx kolumner österut och dy rader söderut. Bilinjär interpolation bevarar
summan och flyttar ljusets tyngdpunkt exakt (−dx, −dy), utom vid kanterna.
"""

from __future__ import annotations

import numpy as np

TILE = 600            # 5 grader
MAX_SHIFT = 3
MIN_LIT = 1500


def _prep(a: np.ndarray) -> np.ndarray:
    x = np.log1p(np.maximum(a.astype("float64"), 0))
    x -= x.mean()
    w = np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    return x * w


def _parabola(cm: float, c0: float, cp: float) -> float:
    den = cm - 2 * c0 + cp
    return 0.0 if den == 0 else 0.5 * (cm - cp) / den


def tile_offset(ref: np.ndarray, img: np.ndarray) -> tuple[float, float, float] | None:
    """(dx, dy, toppkorrelation) för en ruta, eller None om toppen ligger på
    kanten av sökområdet (ingen tydlig topp)."""
    a, b = _prep(ref), _prep(img)
    fa, fb = np.fft.rfft2(a), np.fft.rfft2(b)
    cc = np.fft.irfft2(np.conj(fa) * fb, s=a.shape)        # cc[s] = sum a(x) b(x+s)
    cc /= np.sqrt((a * a).sum() * (b * b).sum()) or 1.0
    m = MAX_SHIFT
    win = np.roll(np.roll(cc, m, 0), m, 1)[: 2 * m + 1, : 2 * m + 1]   # s = -m..m
    iy, ix = np.unravel_index(np.argmax(win), win.shape)
    if iy in (0, 2 * m) or ix in (0, 2 * m):
        return None
    dy = iy - m + _parabola(win[iy - 1, ix], win[iy, ix], win[iy + 1, ix])
    dx = ix - m + _parabola(win[iy, ix - 1], win[iy, ix], win[iy, ix + 1])
    return float(dx), float(dy), float(win[iy, ix])


def weighted_median(v: np.ndarray, w: np.ndarray) -> float:
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return float(v[o][np.searchsorted(cw, 0.5 * cw[-1])])


def shift_bilinear(a: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """a'(r, c) = a(r + dy, c + dx), bilinjärt; utanför bilden räknas som 0.

    Anroparen läser en marginal på minst ceil(|d|)+1 pixlar och skär bort den,
    så att kanterna inte påverkar resultatet inne i fönstret.
    """
    out = a.astype("float32")
    for axis, d in ((1, dx), (0, dy)):
        if d == 0:
            continue
        k = int(np.floor(d))
        f = np.float32(d - k)
        n = out.shape[axis]
        pad = abs(k) + 2
        p = np.pad(out, [(pad, pad) if ax == axis else (0, 0) for ax in range(2)])
        lo = np.take(p, np.arange(n) + pad + k, axis=axis)
        hi = np.take(p, np.arange(n) + pad + k + 1, axis=axis)
        out = (1 - f) * lo + f * hi
    return out
