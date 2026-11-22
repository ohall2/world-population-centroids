#!/usr/bin/env python3
"""
Bygg läsexemplaret som PDF: manus med figurerna på sina platser, bilaga och
bilagetabeller, och först en ruta med länkar till data, kod och dokumentation.

Läser docs/manuscript/{manuscript,supplement,supplement_tables}.md och
figurerna i docs/figures/*.png. Kräver pandoc och xelatex (med tcolorbox och
pdflscape). Länkarna pekar på den gren som anges med --branch.

Exempel
  python tools/build_review_pdf.py
  python tools/build_review_pdf.py --data-link "https://zenodo.org/records/...?preview=1&token=..."
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
        5: "fig5_direction_by_shift", 6: "fig6_sources", 7: "fig7_change_overview"}
WIDTH = {3: "95%", 5: "80%"}
# Relativa kolumnbredder för manusets tre tabeller (rubrikradens början -> bredder).
TABLES = {"| Benchmark | Census years |": [30, 12, 30, 28],
          "| Benchmark | Units | Light |": [22, 9, 9, 13, 11, 11, 13],
          "| Benchmark | Level | Period |": [17, 6, 11, 7, 9, 15, 17, 17, 8]}


def esc(text: str) -> str:
    return text.replace("*", r"\*").replace("_", r"\_")


def shift_headings(text: str) -> str:
    out = []
    for line in text.split("\n"):
        if line.startswith("### "):
            line = "## " + line[4:]
        elif line.startswith("## "):
            line = "# " + line[3:]
        out.append(line)
    return "\n".join(out)


def table_notes(text: str) -> str:
    """Kursiva tabellnoter med filnamn: skydda *, _ och <> så att pandoc inte tolkar dem."""
    lines = text.split("\n")
    for i, l in enumerate(lines):
        if l.startswith("*") and l.endswith("*") and not l.startswith("**") and (".csv" in l or "<" in l):
            inner = re.sub(r"(?<!\\)([_<>*])", r"\\\1", l[1:-1])
            lines[i] = "*" + inner + "*"
    return "\n".join(lines)


def manuscript(src: str) -> tuple[str, str]:
    title = re.search(r"\*\*(Can Nighttime Lights Track Where People Move\?[^*]*)\*\*", src).group(1)
    m = src[src.index("**Abstract.**"):]
    m, n = re.subn(r"```latex\n(.*?)\n```", lambda g: "$$" + g.group(1).strip() + "$$", m, flags=re.S)
    assert n == 1, "formeln hittades inte"
    c0, r0 = m.index("## Data Availability and Figure Captions"), m.index("## References")
    capsec = m[c0:r0]
    caps = {}
    for k in FIGS:
        g = re.search(r"\*\*Figure %d\.\*\* (.*?)(?=\n\n|\Z)" % k, capsec, flags=re.S)
        assert g, f"bildtext {k} saknas"
        caps[k] = g.group(1).strip()
    da = re.search(r"\*\*Data availability\.\*\*.*?(?=\n\n)", capsec, flags=re.S).group(0)
    m = m[:c0] + "## Data Availability\n\n" + da + "\n\n" + m[r0:]
    for k in sorted(FIGS):  # figuren efter stycket där den först nämns
        g = re.search(r"Figure %d[ab]?\b" % k, m[:m.index("## Data Availability")])
        assert g, f"figur {k} nämns inte"
        j = m.find("\n\n", g.start())
        img = "\n\n![%s](%s){width=%s}" % (esc("Figure %d. %s" % (k, caps[k])),
                                           R / "docs/figures" / (FIGS[k] + ".png"), WIDTH.get(k, "100%"))
        m = m[:j] + img + m[j:]
    for head, w in TABLES.items():
        i = m.index(head)
        j = m.find("\n", i)
        k = m.find("\n", j + 1)
        m = m[:j + 1] + "|" + "|".join("-" * x for x in w) + "|" + m[k:]
    return title, shift_headings(m)


def supplement(src: str) -> str:
    g = re.search(r"\*\*Figure S1\.\*\* (.*?)(?=\n\n|\Z)", src, flags=re.S)
    cap = ("Figure S1. " + g.group(1).strip()).replace("*", r"\*")
    s = src[:g.start()] + "![%s](%s){width=100%%}" % (cap, R / "docs/figures/figS1_seam.png") + src[g.end():]
    s = s[:s.index("## S3. Supplementary Tables")].rstrip() + "\n"
    return s.replace('# Supplement to "Can Nighttime Lights Track Where People Move?"', "# Supplement")


def box(data_link: str | None, branch: str) -> str:
    blob, tree = f"{REPO}/blob/{branch}", f"{REPO}/tree/{branch}"
    data = (f"[Zenodo, draft of version 2.0.0]({data_link}) – a private link, please do not pass it on" if data_link
            else "Zenodo draft, reserved DOI 10.5281/zenodo.24336701; the private link is sent with this copy")
    return ("::: {.box}\n"
            "**Reading copy, " + dt.date.today().strftime("%d %B %Y").lstrip("0") + ". Not for circulation.**\n\n"
            f"- **Start here:** [guide for readers]({blob}/READERS.md) – where everything is and what can be rerun\n"
            f"- **Data** (not yet public): {data}\n"
            f"- **Code, tests and result tables:** [GitHub, branch {branch}]({tree})\n"
            f"- **Census benchmarks and comparison centroids per unit:** [replication/]({tree}/replication)\n"
            f"- **Raw data:** not redistributed; sources and download commands in [docs/DATA.md]({blob}/docs/DATA.md)\n"
            f"- **Analysis plan and deviations:** [docs/ANALYSIS_PLAN.md]({blob}/docs/ANALYSIS_PLAN.md)\n"
            f"- **Reviews and the audit of every number:** [docs/review/]({tree}/docs/review)\n"
            ":::\n\n")


TEX_HEADER = r"""
\usepackage{caption}
\captionsetup[figure]{labelformat=empty,font=small,skip=6pt}
\usepackage{float}
\floatplacement{figure}{H}
\usepackage{pdflscape}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\footnotesize}
\setlength{\LTleft}{0pt}
\usepackage{microtype}
\usepackage{tcolorbox}
\usepackage{xurl}
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
    p = argparse.ArgumentParser(description="Bygg läsexemplaret som PDF.")
    p.add_argument("--data-link", help="Zenodos förhandslänk; utan den hänvisar rutan till följebrevet.")
    p.add_argument("--branch", default="utokade-kontroller", help="Gren som länkarna pekar på.")
    p.add_argument("--font", default="TeX Gyre Pagella")
    p.add_argument("--out", type=Path, default=R / "docs/manuscript/manuscript_review.pdf")
    a = p.parse_args(argv)
    M = R / "docs/manuscript"
    title, body = manuscript((M / "manuscript.md").read_text())
    sup = supplement((M / "supplement.md").read_text())
    tabs = (M / "supplement_tables.md").read_text().replace("# Supplementary Tables S1–S6", "# Supplementary Tables")
    front = '---\ntitle: "%s"\nauthor: "Ola Hall, Lund University"\n---\n\n' % title.replace('"', '\\"')
    doc = table_notes(front + box(a.data_link, a.branch) + body.rstrip() + "\n\n\\newpage\n\n" + sup
                      + "\n\n```{=latex}\n\\newpage\n\\begin{landscape}\n```\n\n" + tabs.rstrip()
                      + "\n\n```{=latex}\n\\end{landscape}\n```\n")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "m.md").write_text(doc)
        (td / "h.tex").write_text(TEX_HEADER)
        (td / "box.lua").write_text(LUA)
        subprocess.run(["pandoc", str(td / "m.md"), "-o", str(a.out), "--pdf-engine=xelatex",
                        "-H", str(td / "h.tex"), "--lua-filter", str(td / "box.lua"), "--toc", "--toc-depth=1",
                        "-V", f"mainfont={a.font}", "-V", "monofont=DejaVu Sans Mono", "-V", "fontsize=11pt",
                        "-V", "linestretch=1.15", "-V", "geometry:a4paper", "-V", "geometry:margin=2.5cm",
                        "-V", "colorlinks=true", "-V", "linkcolor=black", "-V", "urlcolor=blue"], check=True)
    print(f"skrev {a.out}")


if __name__ == "__main__":
    main()
