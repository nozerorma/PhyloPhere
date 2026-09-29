"""Chunks of a gene arrive from the worker pool in any order; their merge must not depend on it."""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.utils.gene_wrapper import _merge_gene_chunks  # noqa: E402

CHUNKS = [
    [("b_3", ["r3a", "r3b"]), ("b_4", ["r4"])],
    [("b_1", ["r1a", "r1b"])],
    [("b_10", ["r10"]), ("b_2", ["r2"])],
]


def test_merge_is_independent_of_arrival_order():
    merged = {tuple(map(str, _merge_gene_chunks([c for chunk in order for c in chunk])))
              for order in itertools.permutations(CHUNKS)}
    assert len(merged) == 1


def test_merge_orders_by_cycle_tag_and_keeps_records_within_a_cycle():
    out = _merge_gene_chunks([c for chunk in CHUNKS for c in chunk])
    # plain tag order, the order a single-chunk replay of the sorted tags produces
    assert [tag for tag, _ in out] == sorted(tag for chunk in CHUNKS for tag, _ in chunk)
    assert dict(out)["b_3"] == ["r3a", "r3b"] and dict(out)["b_1"] == ["r1a", "r1b"]
