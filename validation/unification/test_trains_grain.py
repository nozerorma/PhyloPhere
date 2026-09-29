"""trains_grain.py: train flags under the union of hypotheses versus per hypothesis."""
import random
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trains_grain import measure  # noqa: E402


def _table(rows):
    return pd.DataFrame(rows, columns=["gene", "caap_group", "trait", "position"])


def test_positions_that_only_form_a_train_together_are_flagged_by_the_union_alone():
    # each hypothesis has two positions (below minlen 3); together they are four in a row
    d = _table([("g", "US", "H1", 10), ("g", "US", "H1", 12), ("g", "US", "H2", 11), ("g", "US", "H2", 13)])
    m = measure(d, ["H1", "H2"], maxcaas=0.7, minlen=3)
    assert m["positions"] == 4 and m["flagged_union"] == 4 and m["flagged_any_hyp"] == 0
    assert m["only_union"] == 4 and m["only_hyp"] == 0
    assert m["records"] == 4 and m["records_removed_union"] == 4 and m["records_removed_hyp"] == 0
    assert m["units_with_train_union"] == 1 and m["units_with_train_hyp"] == 0


def test_a_train_inside_one_hypothesis_is_found_by_both_grains():
    d = _table([("g", "US", "H1", p) for p in (5, 6, 7)] + [("g", "US", "H2", 50)])
    m = measure(d, ["H1", "H2"], maxcaas=0.7, minlen=3)
    assert m["flagged_union"] == 3 and m["flagged_any_hyp"] == 3 and m["only_union"] == 0
    assert m["records_removed_union"] == 3 and m["records_removed_hyp"] == 3


def test_one_hypothesis_makes_the_grains_identical_and_hyp_flags_never_exceed_the_union():
    rng = random.Random(5)
    rows = [(f"g{rng.randint(0, 3)}", "US", f"H{rng.randint(1, 6)}", rng.randint(1, 40)) for _ in range(600)]
    d = _table(rows).drop_duplicates()
    for h in ("H1", "H3"):
        m = measure(d, [h], maxcaas=0.7, minlen=3)
        assert m["only_union"] == 0 and m["flagged_union"] == m["flagged_any_hyp"]
    m = measure(d, [f"H{i}" for i in range(1, 7)], maxcaas=0.7, minlen=3)
    assert m["only_hyp"] == 0 and m["flagged_union"] >= m["flagged_any_hyp"]
    assert m["records_removed_union"] >= m["records_removed_hyp"]
