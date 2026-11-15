"""Förloppsvisning: snurra, förfluten tid och, när totalen är känd, ETA.

Används i alla skript så att det alltid syns att de arbetar och inte har stannat.
"""

from __future__ import annotations

import sys
import threading
import time

_CLEAR = "\033[K"  # rensa resten av raden, så att gammal text inte blir kvar


def fmt_seconds(s: float) -> str:
    s = int(round(s))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{sec:02d}" if h else f"{m:d}:{sec:02d}"


class Progress:
    """Visar en rad som uppdateras fyra gånger per sekund.

    with Progress("Beräknar adm2 1992", total=47217) as p:
        for ...:
            p.advance()
    """

    def __init__(self, text: str, total: int | None = None, stream=sys.stdout):
        self.text, self.total, self.done = text, total, 0
        self.t0 = time.time()
        self.stream = stream
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def advance(self, n: int = 1) -> None:
        with self._lock:
            self.done += n

    def set_text(self, text: str) -> None:
        with self._lock:
            self.text = text

    def _line(self, symbol: str) -> str:
        elapsed = time.time() - self.t0
        with self._lock:
            text, done, total = self.text, self.done, self.total
        if total:
            pct = 100 * done / total
            eta = (elapsed / done * (total - done)) if done else 0
            tail = f"{done}/{total} ({pct:5.1f} %)  tid {fmt_seconds(elapsed)}"
            if done and done < total:
                tail += f"  kvar ca {fmt_seconds(eta)}"
        else:
            tail = f"tid {fmt_seconds(elapsed)}"
        return f"\r  {symbol} {text}  {tail}{_CLEAR}"

    def _run(self) -> None:
        spin, i = "|/-\\", 0
        while not self._stop.wait(0.25):
            self.stream.write(self._line(spin[i % 4]))
            self.stream.flush()
            i += 1

    def __enter__(self) -> "Progress":
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        self._thread.join()
        symbol = "✗" if exc_type else "✓"
        self.stream.write(self._line(symbol) + "\n")
        self.stream.flush()


def log(msg: str, stream=sys.stdout) -> None:
    stream.write(f"{time.strftime('%H:%M:%S')} {msg}\n")
    stream.flush()
