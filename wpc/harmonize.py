"""
Harmonisering av VIIRS till DMSP-liknande värden.

Metod (i linje med Li m.fl. 2020, Scientific Data 7:168, men icke-parametrisk):

1. Aggregering till DMSP:s rutnät. VIIRS-cellernas mittpunkter ligger på
   -180 + k/240 och 75 - k/240, DMSP:s på -180 + i/120 och 75 - j/120. En
   DMSP-cell täcker alltså en hel VIIRS-cell (2j, 2i) och halva grannceller.
   Ytviktat medelvärde = faltning med [1/4, 1/2, 1/4] i båda led, tagen i
   jämna index. Cellernas olika area inom en DMSP-cell försummas (relativt
   fel < 1e-5).
2. Utsmetning. DMSP-OLS ser ljus utsmetat över flera kilometer (blooming).
   Den aggregerade VIIRS-bilden faltas med en gaussisk kärna med
   standardavvikelsen sigma (DMSP-pixlar). sigma väljs med korsvalidering.
3. Överföringsfunktion i två steg ("hurdle"). För överlappsåret 2013 och
   för klasser av log10(utsmetad VIIRS-radians) beräknas (a) andelen pixlar
   som är tända i DMSP-IRQR (F18 2013, kalibrerad) och (b) medelvärdet av
   DMSP-värdet bland de tända. Båda görs monotona (pool-adjacent-violators,
   viktad med antal pixlar). En pixel får värdet (b) om andelen (a) är minst
   1/2, annars 0. Det efterliknar DMSP:s tröskel: bakgrunden blir exakt 0,
   i stället för ett litet positivt förväntat värde överallt (som skulle dra
   centroiderna mot ytornas geometriska mitt). Värdet interpoleras linjärt
   i log10 mellan klassernas mittpunkter.
4. Tillämpning. Samma sigma och funktion används för alla VIIRS-år; VIIRS
   radianser är kalibrerade ombord, så funktionen behöver inte ändras.

Pixlar inom fackelmasken tas inte med när funktionen anpassas, och de sätts
till 0 i resultatet (masken används också när centroiderna beräknas).
"""

from __future__ import annotations

import numpy as np

# Klassgränser för log10(radians). Radians under 10^LOG_MIN räknas som 0.
LOG_MIN, LOG_MAX, LOG_STEP = -2.0, 4.6, 0.025
N_BINS = int(round((LOG_MAX - LOG_MIN) / LOG_STEP))
ZERO_BIN = 0                 # klass 0 = radians < 10^LOG_MIN (inklusive 0)


def aggregate_to_dmsp(v: np.ndarray) -> np.ndarray:
    """VIIRS-block -> DMSP-block, ytviktat.

    `v` har alla 86 401 VIIRS-kolumner och raderna 2a-1 ... 2b+1, dvs. de
    jämna rader som motsvarar DMSP-raderna a ... b plus en udda rad på var
    sida (rader utanför bilden fylls med 0 av anroparen). Resultatet har
    raderna a ... b och alla 43 201 DMSP-kolumner. Kolumner utanför bilden
    (väster om -180, öster om 180) räknas som 0.
    """
    if v.shape[0] % 2 == 0:
        raise ValueError("VIIRS-blocket ska ha ett udda antal rader (2(b-a)+3)")
    r = 0.25 * v[:-2:2] + 0.5 * v[1:-1:2] + 0.25 * v[2::2]
    p = np.pad(r, ((0, 0), (1, 1)))
    return 0.25 * p[:, 0:-2:2] + 0.5 * p[:, 1:-1:2] + 0.25 * p[:, 2::2]


def gaussian_kernel(sigma: float) -> np.ndarray:
    if sigma <= 0:
        return np.array([1.0])
    rad = int(np.ceil(3 * sigma))
    x = np.arange(-rad, rad + 1)
    k = np.exp(-0.5 * (x / sigma) ** 2)
    return k / k.sum()


def blur(a: np.ndarray, sigma: float) -> np.ndarray:
    """Separerbar gaussisk faltning; utanför blocket räknas som 0.

    Anroparen läser en marginal på minst 3*sigma rader runt blocket och
    skär bort den efteråt, så att blockkanterna inte påverkar resultatet.
    """
    k = gaussian_kernel(sigma)
    if k.size == 1:
        return a.astype("float32", copy=True)
    rad = k.size // 2
    p = np.pad(a.astype("float32"), rad)
    rows = np.zeros((a.shape[0], p.shape[1]), "float32")
    for o, w in enumerate(k):
        rows += w * p[o:o + a.shape[0], :]
    out = np.zeros(a.shape, "float32")
    for o, w in enumerate(k):
        out += w * rows[:, o:o + a.shape[1]]
    return out


def bin_index(x: np.ndarray) -> np.ndarray:
    """Klass för radians x: 0 för x < 10^LOG_MIN, annars 1..N_BINS."""
    with np.errstate(divide="ignore"):
        lx = np.log10(np.maximum(x, 1e-30))
    i = np.floor((lx - LOG_MIN) / LOG_STEP).astype(np.int64) + 1
    i[lx < LOG_MIN] = ZERO_BIN
    return np.clip(i, 0, N_BINS)


def bin_centers() -> np.ndarray:
    """log10-mittpunkt för klasserna 1..N_BINS (klass 0 saknar mittpunkt)."""
    return LOG_MIN + (np.arange(N_BINS) + 0.5) * LOG_STEP


class BinStats:
    """Löpande summor per klass av DMSP-värdet y: antal, summa, kvadratsumma,
    antal tända (y > 0) och summa bland tända."""

    def __init__(self):
        z = lambda: np.zeros(N_BINS + 1)  # noqa: E731
        self.n, self.s, self.ss, self.nlit, self.slit = z(), z(), z(), z(), z()

    def add(self, idx: np.ndarray, y: np.ndarray) -> None:
        idx = idx.ravel()
        y = y.ravel().astype("float64")
        lit = y > 0
        bc = lambda i, w=None: np.bincount(i, weights=w, minlength=N_BINS + 1)  # noqa: E731
        self.n += bc(idx)
        self.s += bc(idx, y)
        self.ss += bc(idx, y * y)
        self.nlit += bc(idx[lit])
        self.slit += bc(idx[lit], y[lit])

    def merged(self, other: "BinStats") -> "BinStats":
        out = BinStats()
        for k in ("n", "s", "ss", "nlit", "slit"):
            setattr(out, k, getattr(self, k) + getattr(other, k))
        return out

    def to_dict(self) -> dict:
        return {k: getattr(self, k).tolist() for k in ("n", "s", "ss", "nlit", "slit")}


def pav(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Viktad isoton regression (icke-avtagande), pool-adjacent-violators."""
    vals, wts, sizes = [], [], []
    for yi, wi in zip(y, w):
        vals.append(float(yi)); wts.append(float(wi)); sizes.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, w2, s2 = vals.pop(), wts.pop(), sizes.pop()
            v1, w1, s1 = vals.pop(), wts.pop(), sizes.pop()
            wt = w1 + w2
            vals.append((v1 * w1 + v2 * w2) / wt if wt > 0 else (v1 + v2) / 2)
            wts.append(wt); sizes.append(s1 + s2)
    return np.repeat(vals, sizes)


def _ffill(v: np.ndarray, n: np.ndarray) -> np.ndarray:
    out = v.copy()
    for i in range(1, out.size):
        if n[i] == 0:
            out[i] = out[i - 1]
    return out


def fit_table(st: BinStats) -> dict:
    """Tabeller per klass 0..N_BINS: p_lit, mean_lit och value (= mean_lit om
    p_lit >= 1/2, annars 0). Klass 0 (radians ~ 0) får alltid värdet 0."""
    p = np.divide(st.nlit, st.n, out=np.zeros_like(st.n), where=st.n > 0)
    m = np.divide(st.slit, st.nlit, out=np.zeros_like(st.n), where=st.nlit > 0)
    p = pav(_ffill(p, st.n), st.n)
    m = pav(_ffill(m, st.nlit), st.nlit)
    value = np.where(p >= 0.5, m, 0.0)
    value[ZERO_BIN] = 0.0
    return {"p_lit": p, "mean_lit": m, "value": value}


def apply_table(x: np.ndarray, value: np.ndarray) -> np.ndarray:
    """DMSP-liknande värde för utsmetad radians x.

    Linjär interpolation i log10 mellan klassmittpunkterna; där tabellen är 0
    (andel tända < 1/2) blir resultatet 0, även nära en tröskel.
    """
    out = np.zeros(x.shape, "float32")
    pos = x >= 10 ** LOG_MIN
    if pos.any():
        idx = bin_index(x[pos])
        lx = np.log10(x[pos])
        v = np.interp(lx, bin_centers(), value[1:]).astype("float32")
        v[value[idx] == 0] = 0.0
        out[pos] = v
    return out


def sse(st: BinStats, value: np.ndarray) -> float:
    """Kvadratfelsumma på stickprovet `st` när varje klass förutsägs med
    `value` (beräknad från klassummorna)."""
    return float(np.sum(st.ss - 2 * value * st.s + st.n * value ** 2))
