"""Test av den parvisa läges-jämförelsen i analys B (compare_benchmarks.pairwise)
med konstgjorda data där svaret är känt. Måtten och bootstrapen återanvänds och
prövas inte här."""

import numpy as np
import pytest

from tools.compare_benchmarks import pairwise


def test_pairwise_known():
    err_light = np.array([2.0, 2.0, 2.0, 2.0])
    err_source = np.array([1.0, 3.0, 1.0, 1.0])     # närmare, längre, närmare, närmare
    closer, diff = pairwise(err_light, err_source)
    assert closer == pytest.approx(0.75)             # 3 av 4 enheter
    assert diff == pytest.approx(1.0)                # median av [1, -1, 1, 1]


def test_pairwise_ties_not_counted_as_closer():
    err_light = np.array([1.0, 1.0, 1.0])
    err_source = np.array([1.0, 1.0, 0.5])           # lika räknas inte som närmare
    closer, diff = pairwise(err_light, err_source)
    assert closer == pytest.approx(1 / 3)
    assert diff == pytest.approx(0.0)                # median av [0, 0, 0.5]


# Integrationstest (hoppas över om data saknas, som data/ inte ligger i git).
# Hade fångat felet där WorldPop använde GHS-epoken för geostat 2011, så att
# Europa 2011 föll ur: varje facitår och period ska finnas i utdata för varje
# källa som har data för det året.
from pathlib import Path                                                      # noqa: E402
from tools.compare_benchmarks import run_dataset, DS                          # noqa: E402

_GEOSTAT_DATA = (Path("data/validation/geostat/2011/adm2.csv").exists()
                 and Path("data/results/gadm41/worldpop/adm2/2011.csv").exists()
                 and Path("data/release/v2.0.0").exists())


@pytest.mark.skipif(not _GEOSTAT_DATA, reason="geostat-facit/centroider saknas (data/ ingår inte i git)")
def test_every_facit_year_and_period_present_for_each_source():
    from wpc.lights import available_years
    sources = ["ghspop", "ghsbuilts", "worldpop"]
    epochs = {s: available_years(s, Path("data/raw") / s, None) for s in sources}
    pos, chg, _pair, sel = run_dataset(
        "geostat", DS["geostat"], Path("data/validation"), Path("data/release/v2.0.0"),
        Path("data/results/gadm41"), Path("data/processed/boundaries"), "gadm41",
        sources, epochs, n_boot=50)
    names = ["light", *sources, "geometric"]
    for lv in (1, 2):
        for y in DS["geostat"]["years"]:                 # 2011 och 2021
            for src in names:
                assert ((pos.level == lv) & (pos.year == y) & (pos.source == src)).any(), \
                    f"läge saknas: nivå {lv}, år {y}, källa {src}"
        for y1, y2 in DS["geostat"]["pairs"]:            # 2011–2021
            for src in names:
                assert ((chg.level == lv) & (chg["from"] == y1) & (chg.to == y2) & (chg.source == src)).any(), \
                    f"förändring saknas: nivå {lv}, {y1}-{y2}, källa {src}"
    # geostat 2011 ska inte tappa Europa-enheter pga WorldPop (per-källa-epok)
    s2011 = sel[(sel.level == 2) & (sel.year == 2011)].iloc[0]
    assert s2011.n_used > 0 and s2011.missing_worldpop < s2011.n_facit
