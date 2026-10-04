"""POSENRICH without a phylogenetic permulation null: no p-value under the name of a test it did not run."""
import math

import numpy as np
import pandas as pd
import pytest

from test_null_consumers import _module

TERMS = {"T1": {"g1:1", "g1:2", "g2:1"}, "T2": {"g3:1", "g3:2", "g3:3"}}
DESCS = {"T1": "one", "T2": "two"}
OBS = {"g1:1": 1.0, "g2:1": 2.0, "g3:1": 0.5}
BACKGROUND = ["g1:1", "g1:2", "g2:1", "g3:1", "g3:2", "g3:3", "g4:1", "g4:2"]


@pytest.fixture(scope="module")
def pe():
    return _module("posenrich_enrich")


def _run(pe, **kw):
    return pe.run_permulation_for_terms(TERMS, DESCS, OBS, BACKGROUND, 1, 0, n_perms=200, seed=1998, **kw)


def _null(pe):
    rows = [(p, "global", c, s) for c in ("b_1", "b_2", "b_3", "b_4") for p, s in (("g1:1", 0.3), ("g3:2", 0.1))]
    df = pd.DataFrame(rows, columns=["pos_id", "side", "cycle", "score"])
    return df[["pos_id", "cycle", "score"]], np.array(["b_1", "b_2", "b_3", "b_4"])


def test_without_a_null_the_observed_sums_stay_and_every_null_based_value_is_undefined(pe):
    res = _run(pe)
    assert [r["pathway"] for r in res] == ["T1", "T2"]
    assert [r["obs_sum"] for r in res] == [3.0, 0.5]
    for r in res:
        assert all(math.isnan(r[k]) for k in ("null_mean", "null_sd", "perm_nes", "p_value"))
        assert r["direction"] == "" and r["_overlap"]


def test_with_a_null_the_values_are_defined(pe):
    sub, cycles = _null(pe)
    res = _run(pe, caas_null_sub=sub, caas_null_cycles=cycles)
    assert all(not math.isnan(r["p_value"]) and r["direction"] in ("enriched", "depleted") for r in res)


def test_the_label_shuffle_is_an_explicit_opt_in(pe):
    res = _run(pe, allow_label_shuffle=True)
    assert all(not math.isnan(r["p_value"]) for r in res)


def test_main_leaves_the_adjusted_p_undefined_and_nothing_significant_without_a_null():
    import inspect
    src = inspect.getsource(_module("posenrich_enrich").main)
    assert "no_null" in src and "np.full(len(res), np.nan)" in src
