"""Förloppsvisning: snurra, förfluten tid och, när totalen är känd, ETA.

Används i alla skript så att det alltid syns att de arbetar och inte har stannat.

Raden ritas av en egen process, inte en tråd. Vissa bibliotek (t.ex. GDAL vid
läsning av stora filer) blockerar Pythons tolk i minuter åt gången; en tråd
står då still, men en separat process fortsätter att uppdatera klockan.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time

_CLEAR = "\033[K"  # rensa resten av raden, så att gammal text inte blir kvar
_TEXT_BYTES = 200


def fmt_seconds(s: float) -> str:
    s = int(round(s))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{sec:02d}" if h else f"{m:d}:{sec:02d}"


def _line(symbol: str, text: str, done: int, total: int, elapsed: float) -> str:
    if total > 0:
        pct = 100 * done / total
        tail = f"{done}/{total} ({pct:5.1f} %)  tid {fmt_seconds(elapsed)}"
        if 0 < done < total:
            tail += f"  kvar ca {fmt_seconds(elapsed / done * (total - done))}"
    else:
        tail = f"tid {fmt_seconds(elapsed)}"
    return f"\r  {symbol} {text}  {tail}{_CLEAR}"


def _draw(text, done, total, t0, stop) -> None:
    spin, i = "|/-\\", 0
    while not stop.wait(0.25):
        t = bytes(text[:]).rstrip(b"\0").decode("utf-8", "replace")
        sys.stdout.write(_line(spin[i % 4], t, done.value, total.value, time.time() - t0.value))
        sys.stdout.flush()
        i += 1


class Progress:
    """with Progress("Beräknar adm2 1992", total=47217) as p:
           for ...:
               p.advance()
    """

    def __init__(self, text: str, total: int | None = None):
        ctx = mp.get_context("spawn")
        self._text = ctx.Array("c", _TEXT_BYTES)
        self._done = ctx.Value("q", 0)
        self._total = ctx.Value("q", total or 0)
        self._t0 = ctx.Value("d", time.time())
        self._stop = ctx.Event()
        self.set_text(text)
        self._proc = ctx.Process(target=_draw, daemon=True,
                                 args=(self._text, self._done, self._total, self._t0, self._stop))

    @property
    def done(self) -> int:
        return self._done.value

    def advance(self, n: int = 1) -> None:
        with self._done.get_lock():
            self._done.value += n

    def set_done(self, n: int) -> None:
        self._done.value = n

    def set_text(self, text: str) -> None:
        raw = text.encode("utf-8")[:_TEXT_BYTES - 1]
        self._text[:] = raw + b"\0" * (_TEXT_BYTES - len(raw))

    def __enter__(self) -> "Progress":
        self._proc.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        self._proc.join(timeout=2)
        text = bytes(self._text[:]).rstrip(b"\0").decode("utf-8", "replace")
        symbol = "✗" if exc_type else "✓"
        sys.stdout.write(_line(symbol, text, self._done.value, self._total.value,
                               time.time() - self._t0.value) + "\n")
        sys.stdout.flush()


def log(msg: str) -> None:
    # \r + rensa raden: skriver över en eventuell förloppsrad i stället för
    # att hamna efter den.
    sys.stdout.write(f"\r{_CLEAR}{time.strftime('%H:%M:%S')} {msg}\n")
    sys.stdout.flush()
