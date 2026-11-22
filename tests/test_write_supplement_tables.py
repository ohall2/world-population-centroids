"""tools/write_supplement_tables.py: formatering, tolkningsregeln (planen §7) och
avvikelsetabellen, på syntetiska data med känt svar."""
from pathlib import Path

import pandas as pd
import pytest

from tools import write_supplement_tables as W

R = Path(__file__).resolve().parents[1]


def test_num_and_interval():
    assert W.num(0.6144) == "0.61"
    assert W.num(-0.08) == "−0.08"
    assert W.num(2436, 0) == "2,436"
    assert W.ci(0.6144, 0.5846, 0.6417) == "0.61 (0.58–0.64)"
    assert W.ci(-0.08, -0.16, -0.01) == "−0.08 (−0.16 to −0.01)"
    assert W.ci(0.22, -0.03, 0.47) == "0.22 (−0.03 to 0.47)"
    assert W.num(float("nan")) == "–"


def test_md_table():
    t = W.md_table(["a", "b"], [["1", "x|y"]])
    assert t.splitlines() == ["| a | b |", "| --- | --- |", "| 1 | x\\|y |"]


@pytest.mark.parametrize("est,lo,hi,n,row", [
    (0.70, 0.55, 0.80, 300, 1),     # rad 1
    (0.58, 0.52, 0.64, 300, 2),     # rad 2 (bredd 0,12)
    (0.589, 0.509, 0.680, 287, 3),  # Kanada: bredd 0,17
    (0.58, 0.55, 0.61, 120, 3),     # för få enheter
    (0.4637, 0.4398, 0.4896, 1792, 4),  # Japan: under 0,50, smalt
    (0.65, 0.58, 0.70, 300, 1),     # gränsen 0,65 hör till rad 1
])
def test_rule_outcome(est, lo, hi, n, row):
    assert W.rule_outcome(est, lo, hi, n) == row


def test_rule_outcome_uncovered_returns_none():
    # punktskattning ≥ 0,65 men nedre gräns ≤ 0,50 och smalt: täcks inte av någon rad
    assert W.rule_outcome(0.66, 0.50, 0.64, 300) is None


PLAN = """# Plan

## 10. Avvikelser

| Datum | Avsnitt | Ändring | Skäl |
| --- | --- | --- | --- |
| 2026-11-21 | §2 | a \\| b | c |
| 2026-11-22 | §6/§7 | d | e |
"""


def test_parse_deviations():
    rows = W.parse_deviations(PLAN)
    assert rows == [("2026-11-21", "§2"), ("2026-11-22", "§6/§7")]


def test_deviations_must_match_translation():
    en = [("2026-11-21", "§2", "x", "y")]
    with pytest.raises(RuntimeError):
        W.check_deviations(W.parse_deviations(PLAN), en)
    W.check_deviations(W.parse_deviations(PLAN), en + [("2026-11-22", "§6/§7", "x", "y")])


def test_translation_matches_committed_plan():
    """Den engelska S6 måste följa planens avvikelsetabell rad för rad."""
    rows = W.parse_deviations((R / "docs/ANALYSIS_PLAN.md").read_text())
    W.check_deviations(rows, W.DEVIATIONS_EN)


def test_read_empty_csv(tmp_path):
    p = tmp_path / "x.csv"; p.write_text("")
    assert W.read(p).empty
