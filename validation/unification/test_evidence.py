"""Evidence of a position: the unpooled rows of the scorer, as a table of what each domain of each hypothesis saw."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_observed_b0 import SRC, pepc  # noqa: E402,F401  (the PEPC fixture)
from src.core.observed import analyze_observed, observed_entries, score_observed  # noqa: E402

POSITIONS = {"629", "631"}          # 631 holds a tied residue in the frozen master


def _entries(pepc):
    rows = [r for r in pepc["rows"] if r["position"] in POSITIONS]
    assert rows
    return observed_entries("PEPC", rows)


def _analyze(pepc, **kw):
    return analyze_observed(pepc["ctx"], "PEPC", _entries(pepc), pepc["trait_pairs"], pepc["pss"], 0.1, **kw)


def test_the_unpooled_rows_are_kept_only_when_asked_and_the_pooled_result_is_the_same(pepc):
    plain, diag_plain = _analyze(pepc)
    kept, diag_kept = _analyze(pepc, keep_unpooled=True)
    assert "unpooled" not in diag_plain and "unpooled" in diag_kept
    assert [vars(r) for r in plain] == [vars(r) for r in kept]
    assert [vars(r) for r in plain] == [vars(r) for r in score_observed(pepc["ctx"], "PEPC", _entries(pepc), pepc["trait_pairs"], pepc["pss"], 0.1)]


def test_there_is_one_unpooled_row_per_entry_and_each_carries_its_hypothesis_and_domains(pepc):
    _, diag = _analyze(pepc, keep_unpooled=True)
    unpooled = diag["unpooled"]
    assert len(unpooled) == len(_entries(pepc))
    assert all(r.hypothesis for r in unpooled) and any(r.pair_details for r in unpooled)


# ── the evidence table ───────────────────────────────────────────────────────────────────────────────────────────

from src.core.evidence import EVIDENCE_COLUMNS, evidence_rows, select_top_positions  # noqa: E402


@pytest.fixture(scope="module")
def evidence(pepc):
    results, diag = _analyze(pepc, keep_unpooled=True)
    return results, diag["unpooled"], evidence_rows(diag["unpooled"])


def test_one_row_per_entry_and_domain_with_the_documented_columns(pepc, evidence):
    _, unpooled, rows = evidence
    assert rows and list(rows[0]) == EVIDENCE_COLUMNS
    domains = {r["domain"] for r in rows}
    assert domains == {"1", "2", "3", "4"}                       # PEPC has four Voronoi domains
    assert len(rows) == len(unpooled) * len(domains)             # every domain of every hypothesis is listed, scored or not
    assert {r["msa_pos"] for r in rows} == POSITIONS


def test_the_rows_say_what_the_domain_saw_at_its_mrca_and_at_the_tips(evidence):
    _, _, rows = evidence
    seen = [r for r in rows if r["mrca_node"] != ""]
    assert seen, "some domain has an MRCA"
    r = seen[0]
    assert r["mrca_state"] != "" and 0.0 <= float(r["mrca_posterior"]) <= 1.0
    assert r["top_species"] and r["bottom_species"] and r["top_tip_aa"] and r["bottom_tip_aa"]
    assert all(c in "ACDEFGHIKLMNPQRSTVWY-?" for c in r["top_tip_aa"].replace(",", ""))


def test_the_rows_are_ordered_by_position_group_hypothesis_number_and_domain(evidence):
    _, _, rows = evidence
    import re
    key = lambda r: (int(r["msa_pos"]), r["caap_group"], int(re.search(r"\d+", r["hypothesis"]).group()), int(r["domain"]))
    assert [key(r) for r in rows] == sorted(key(r) for r in rows)


def test_the_domain_scores_pool_to_the_master_by_the_mean_over_hypotheses(pepc, evidence):
    """The master's domain_<d>_score is the mean over the hypotheses of the position of the domain's score on that side
    (0 where the domain did not change): the table must reproduce it from its own rows."""
    import pandas as pd
    from src.core.master import write_master_csv
    from src.core.observed import observed_master_rows
    results, _, rows = evidence
    out = Path(pepc["dir"]) / "m.csv"
    write_master_csv(observed_master_rows("PEPC", results, pepc["fields"]), out, pepc["fields"])
    master = pd.read_csv(out, keep_default_na=False)
    checked = 0
    for m in master.itertuples():
        sel = [r for r in rows if r["msa_pos"] == str(m.msa_pos) and r["caap_group"] == m.caap_group]
        hyps = {r["hypothesis"] for r in sel}
        if m.side not in ("top", "bottom"):
            continue
        for d in range(1, 5):
            vals = [float(r[f"{m.side}_domain_score"] or 0.0) for r in sel if r["domain"] == str(d)]
            assert len(vals) == len(hyps)
            expect = sum(vals) / len(vals)
            got = getattr(m, f"domain_{d}_score")
            assert (got if got != "" else 0.0) == pytest.approx(expect, abs=1e-12), (m.msa_pos, m.caap_group, m.side, d)
            checked += 1
    assert checked > 20


def test_the_rows_of_a_position_carry_the_tag_of_their_entry(pepc, evidence):
    _, unpooled, rows = evidence
    assert {r["tag"] for r in rows} == {e.tag for e in _entries(pepc)}


# ── the choice of the positions ──────────────────────────────────────────────────────────────────────────────────

def _score(gene, pos, side, score, p=None):
    return {"Gene": gene, "Position": str(pos), "side": side, "CAAS_score": repr(score), "p.emp": "NA" if p is None else repr(p)}


def test_the_top_positions_are_the_best_by_score_across_sides_with_ties_to_the_smaller_p_then_the_name():
    rows = [_score("B", 5, "top", 0.9), _score("B", 5, "bottom", 0.95),          # one position: its best side counts
            _score("A", 7, "top", 0.9, 0.2), _score("A", 8, "top", 0.9, 0.1),     # tie on 0.9: the smaller p.emp first
            _score("C", 1, "top", 0.9, 0.1),                                      # tie on score and p: gene, then position
            _score("A", 2, "top", 0.1)]
    got = select_top_positions(rows, 4)
    assert [(g, p) for g, p, _ in got] == [("B", "5"), ("A", "8"), ("C", "1"), ("A", "7")]


def test_fewer_positions_than_asked_returns_them_all_and_zero_returns_none():
    rows = [_score("A", 1, "top", 0.5), _score("A", 2, "top", 0.4)]
    assert len(select_top_positions(rows, 10)) == 2 and select_top_positions(rows, 0) == []


def test_a_position_without_a_score_is_not_a_candidate():
    rows = [_score("A", 1, "top", 0.5), {"Gene": "A", "Position": "2", "side": "top", "CAAS_score": "NA", "p.emp": "NA"}]
    assert [(g, p) for g, p, _ in select_top_positions(rows, 5)] == [("A", "1")]
