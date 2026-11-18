#!/usr/bin/env python3
"""
Bygg en läsversion av manuset som PDF: text, tabeller och figurerna på sina
platser, och först en ruta med länkar till data, kod och dokumentation.

Läser docs/manuscript/manuscript.md (export från manusdokumentet) och
figurerna i docs/figures/*.pdf. Kräver pandoc och xelatex.

Exempel
  python tools/build_review_pdf.py --data-link "https://zenodo.org/records/...?preview=1&token=..." \\
      --out docs/manuscript/manuscript_review.pdf
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import subprocess
import tempfile
from pathlib import Path

R = Path(__file__).resolve().parents[1]
REPO = "https://github.com/ohall2/world-population-centroids"
FIGS = {1: "fig1_location_error", 2: "fig2_map", 3: "fig3_national_bias", 4: "fig4_measured_vs_true",
        5: "fig5_direction_by_shift", 6: "fig6_seam"}
WIDTH = {1: "0.85", 2: "1.0", 3: "1.0", 4: "1.0", 5: "0.8", 6: "0.85"}


def build_md(src: str, data_link: str) -> str:
    lines = src.split("\n")
    title = re.search(r"^\*\*(Can Nighttime.*?)\*\*$", src, re.M).group(1)
    body = src.split("**Abstract.**", 1)[1]
    body = "**Abstract.**" + body
    body = re.sub(r" \\\[\d+ words\\\]", "", body)
    body = re.sub(r"```latex\n(.*?)\n```", lambda m: "$$" + m.group(1).strip() + "$$", body, flags=re.S)
    # bildtexter: plocka ut och ta bort avsnittet med dem
    caps = {int(m.group(1)): m.group(2).strip() for m in
            re.finditer(r"^\*\*Figure (\d)\.\*\* (.+)$", body, re.M)}
    body = re.sub(r"^\*\*Figure \d\.\*\* .+\n?", "", body, flags=re.M)
    body = body.replace("## Data Availability and Figure Captions", "## Data Availability")
    # intern anteckning över referenslistan
    body = re.sub(r"^\*Every entry below.*\n", "", body, flags=re.M)
    # figurerna efter stycket där de först nämns
    paras = body.split("\n\n")
    placed = set()
    out = []
    for p in paras:
        out.append(p)
        for n in sorted(FIGS):
            if n not in placed and re.search(rf"\(Figure {n}[a-z]?[\);]", p) and not p.startswith("#"):
                cap = caps[n].replace("[", "(").replace("]", ")")
                out.append(f"![**Figure {n}.** {cap}]({R / 'docs/figures' / (FIGS[n] + '.pdf')})"
                           f"{{width={WIDTH[n]}\\textwidth}}")
                placed.add(n)
    body = "\n\n".join(out)
    box = (
        "::: {.box}\n"
        "**Review copy, " + dt.date.today().strftime("%d %B %Y").lstrip("0") + ". Not for circulation.**\n\n"
        f"- **Data** (draft, not yet public): [Zenodo preview]({data_link})\n"
        f"- **Code and full documentation:** [GitHub, branch v2]({REPO}/tree/v2)\n"
        f"- **Analyses and decisions, step by step:** [docs/HARMONIZATION.md]({REPO}/blob/v2/docs/HARMONIZATION.md)\n"
        f"- **Literature with verification status:** [docs/LITERATURE.md]({REPO}/blob/v2/docs/LITERATURE.md)\n"
        f"- **Independent review of this draft and the changes made:** "
        f"[docs/review/review-manuscript-2026-11-18.md]({REPO}/blob/v2/docs/review/review-manuscript-2026-11-18.md)\n"
        f"- **Robustness checks:** [docs/robustness/]({REPO}/tree/v2/docs/robustness)\n"
        ":::\n\n")
    head = (f"---\ntitle: \"{title}\"\nauthor: \"Ola Hall, Lund University\"\n"
            f"date: \"Draft for Annals of the American Association of Geographers\"\n---\n\n")
    return head + box + body


TEX_HEADER = r"""
\usepackage{fontspec}
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\usepackage[margin=2.3cm]{geometry}
\usepackage{caption}
\captionsetup{font=small,labelformat=empty}
\usepackage{tcolorbox}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\scriptsize}
\setlength{\parskip}{0.4em}
\usepackage{xcolor}
\usepackage{xurl}
\AtBeginDocument{\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black}}
"""

LUA = r"""
function Div(el)
  if el.classes:includes('box') then
    return {pandoc.RawBlock('latex', '\\begin{tcolorbox}[colback=gray!8,colframe=gray!50]')}
      .. el.content .. {pandoc.RawBlock('latex', '\\end{tcolorbox}')}
  end
end
"""


def main(argv=None):
    p = argparse.ArgumentParser(description="Bygg läsversionen av manuset som PDF.")
    p.add_argument("--manuscript", type=Path, default=R / "docs/manuscript/manuscript.md")
    p.add_argument("--data-link", required=True)
    p.add_argument("--out", type=Path, default=R / "docs/manuscript/manuscript_review.pdf")
    a = p.parse_args(argv)
    md = build_md(a.manuscript.read_text(), a.data_link)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "m.md").write_text(md)
        (td / "h.tex").write_text(TEX_HEADER)
        (td / "box.lua").write_text(LUA)
        subprocess.run(["pandoc", str(td / "m.md"), "-o", str(a.out), "--pdf-engine=xelatex",
                        "-H", str(td / "h.tex"), "--lua-filter", str(td / "box.lua"),
                        "-V", "fontsize=10pt", "-V", "papersize=a4", "--resource-path", str(R)], check=True)
    print(f"skrev {a.out}")


if __name__ == "__main__":
    main()
