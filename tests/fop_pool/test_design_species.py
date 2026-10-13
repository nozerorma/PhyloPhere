#!/usr/bin/env python3
# test_design_species.py: the species count S of a design reaches the pooler the same way for the observed design and the null cycles.
# PhyloPhere | tests/fop_pool/
"""Checks core.labelings.read_design_species / design_species_count and core.driver.pool_labelings."""

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"))
from src.core.driver import pool_labelings  # noqa: E402
from src.core.labelings import design_species_count, read_design_species, trait_pairs_from  # noqa: E402
from src.convergence.disambiguate_single import PositionAxes  # noqa: E402

ROWS = [("b_1", "H1", 1, "a1", "a2"), ("b_1", "H2", 1, "a1", "a3"), ("b_2", "H1", 1, "b1", "b2")]


def axes_record(species):
    """A scored position of a labeling: the raw `sides` of two converging domains (1 and 2)."""
    side = {"domain_scores": {1: 0.9, 2: 0.9}, "domain_der": {1: "A", 2: "A"}, "domain_der_enc": {1: "A", 2: "A"},
            "domain_anc": {1: "V", 2: "V"}, "domain_der_species": species, "pair_lca": [[1, 2, 9, 0.9]]}
    return PositionAxes(position=5, caap_group="US", asr_path_score=0.0, side=None, hypothesis=None,
                        domain_scores=None, derived_agreement=None,
                        sides={"top": side, "bottom": {}, "domain_meta": {1: {}, 2: {}}},
                        n_hypotheses=None, n_species=None)


def main():
    with tempfile.TemporaryDirectory(dir=Path.home() / ".cc") as td:
        path = Path(td) / "fop_pairs.tsv"
        path.write_text("cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n"
                        + "".join(f"{c}\t{h}\t{p}\t{s1}\t{s2}\t1.0\n" for c, h, p, s1, s2 in ROWS))
        # Distinct species per cycle, the species shared by two hypotheses counted once.
        assert read_design_species(path) == {"b_1": 3, "b_2": 2}
        assert read_design_species(path.parent / "missing.tsv") == {} and read_design_species("NO_FILE") == {}
        # contrast_hypotheses_pairs.tsv has no cycle column: its rows belong to the observed cycle.
        obs = Path(td) / "contrast_hypotheses_pairs.tsv"
        obs.write_text("hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\nH1\t1\ta1\ta2\t1\nH2\t1\ta1\ta3\t1\n")
        assert list(read_design_species(obs).values()) == [3]

    assert design_species_count({1: [("a1", "a2"), ("a3", "a1")], 2: [("a2", "a4")]}) == 4

    # A null cycle with and without fop_pairs, and the observed design, give the same S and so the same score.
    sp = {1: "a1", 2: "a2"}
    with_fop = pool_labelings([("b_1~H1", [axes_record(sp)])], {"b_1": {("H1", 1): 1.0, ("H1", 2): 1.0}}, {"b_1": 20})
    labelings = {"b_1": ([f"f{i}" for i in range(10)], [f"g{i}" for i in range(10)])}   # S = 20 species
    assert design_species_count(trait_pairs_from(*labelings["b_1"])) == 20
    no_fop = pool_labelings([("b_1", [axes_record(sp)])], None, None, labelings)
    expected = 0.9 * 2 / 20                      # K = 2 domains, both at 0.9: strength 0.9; breadth 2 / 20
    for pooled in (with_fop, no_fop):
        score = pooled[0][1][0].asr_path_score
        assert math.isclose(score, expected), score
    # Without S the breadth is 1.
    assert math.isclose(pool_labelings([("b_1", [axes_record(sp)])], None)[0][1][0].asr_path_score, 0.9)
    print("test_design_species: OK")


if __name__ == "__main__":
    main()
