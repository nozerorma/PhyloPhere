#!/usr/bin/env python3
"""Regenerate the T0 golden fixtures (roadmap tier T0).

Writes two files next to this script:

* ``path_scores_golden.json`` — synthetic scenarios + the full frozen return
  dict of ``compute_asr_path_score`` for each. ``test_path_scores_golden.py``
  replays the scenarios and diffs against these values.
* ``fop_pool_fixture.json`` — shared hypothesis-pooling scenarios + the frozen
  ``pool_hypotheses`` output. ``test_fop_pool_golden.py`` and
  ``test_fop_pool_golden.R`` both read this file and assert their pooling matches
  the stored ``expected`` to 1e-9 — the single shared oracle that keeps the
  Python / R twins from drifting.

Run from anywhere:  python gen_golden.py
Each roadmap tier that changes the scoring maths re-runs this and commits the
diff, rather than hand-editing asserts.
"""
from __future__ import annotations

import json
from pathlib import Path

from _common import FIXTURE_DIR, _load, run_scenario, scenario_to_json

ps = _load("path_scores")
fop = _load("fop_pool")

# A(nc) dominant, V / T / L as derived candidates. Distances kept short so the
# isolation walks stay legible.
A = {"A": 0.9, "V": 0.05, "T": 0.05}
A_soft = {"A": 0.6, "V": 0.3, "T": 0.1}
V_here = {"V": 0.85, "A": 0.1, "T": 0.05}
T_here = {"T": 0.85, "A": 0.1, "V": 0.05}


def _tree_two_clades():
    """root 0 -> {1 -> {3,4}, 2 -> {5,6}}; MRCAs live at 3/4/5/6, deep enough
    that each has >=2 private nodes above it before the clades merge at 1/2/0."""
    edges = [(0, 1), (0, 2), (1, 10), (10, 3), (10, 4),
             (2, 20), (20, 5), (20, 6)]
    node_ids = [0, 1, 2, 10, 20, 3, 4, 5, 6]
    return edges, node_ids


def scenarios():
    e2, n2 = _tree_two_clades()
    base_post = {
        0: A, 1: A, 2: A, 10: A, 20: A,
        3: V_here, 4: V_here, 5: V_here, 6: T_here,
    }

    out = []

    # 1 — single changed pair: core must be 0 (< 2 observations on any side).
    out.append(dict(
        name="single_changed_pair",
        doc="one pair changes top A->V; core=0, asr=0",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[{"pair_id": 1, "node_id": 3, "focal_state": "A",
                       "top_tip_mode": "V", "bottom_tip_mode": "A"}],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 2 — two pairs change on the SAME (top) side, same residue: core > 0.
    out.append(dict(
        name="two_same_side_converge",
        doc="pairs 1 & 2 both change top A->V in separate clades",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 3 — opposite-only: one pair top-only, other bottom-only -> core=0.
    out.append(dict(
        name="opposite_sides_only",
        doc="pair 1 changes top, pair 2 changes bottom; neither side reaches >=2",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "V"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 4 — a pair that changes on BOTH sides, plus a top-only partner.
    out.append(dict(
        name="pair_changes_both_sides",
        doc="pair 1 changes top+bottom to V, pair 2 top-only to V",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "V"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 5 — two converging pairs + a conserved pair -> conservation_gate < 1.
    out.append(dict(
        name="with_conserved_pair",
        doc="pairs 1 & 2 converge top->V; pair 3 held ancestral (conserved)",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 3, "node_id": 4, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=True, conserved_pair="3",
    ))

    # 6 — contamination at hop+1: node directly above pair-1 MRCA is modally V.
    contam_post = dict(base_post)
    contam_post[10] = {"V": 0.7, "A": 0.3}   # parent of nodes 3 and 4
    out.append(dict(
        name="contaminated_hop1",
        doc="derived V already modal at the node above pair-1's MRCA",
        edges=e2, node_ids=n2, posteriors=contam_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 7 — n>2 changed pairs on the top side, mixed residues (agreement < 1).
    e3 = e2 + [(0, 30), (30, 7), (30, 8)]
    n3 = n2 + [30, 7, 8]
    post3 = dict(base_post)
    post3.update({30: A, 7: V_here, 8: T_here})
    out.append(dict(
        name="n_gt_2_mixed_residues",
        doc="4 pairs change top: 3 to V, 1 to T -> derived_agreement 0.75",
        edges=e3, node_ids=n3, posteriors=post3,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 3, "node_id": 6, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 4, "node_id": 7, "focal_state": "A",
             "top_tip_mode": "T", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 8 — GS scheme: V and I co-encode to 'l' under GS3 (A is 'n'), so a V/I
    # split still counts as a change AND agrees. Under US it would disagree.
    post_gs = dict(base_post)
    post_gs[6] = {"I": 0.85, "A": 0.1, "V": 0.05}
    out.append(dict(
        name="gs3_coencoded_agreement",
        doc="pair 1->V, pair 2->I; disagree under US, both 'l' (agree) under GS3",
        edges=e2, node_ids=n2, posteriors=post_gs,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 6, "focal_state": "A",
             "top_tip_mode": "I", "bottom_tip_mode": "A"},
        ],
        scheme="GS3", is_conserved_meta=False, conserved_pair="",
    ))

    # 9 — MRCA at the root: no background above it -> EMPTY_PATH_SCORE regime.
    out.append(dict(
        name="mrca_at_root",
        doc="both MRCAs are the root; path_to_root empty",
        edges=[(0, 1), (0, 2)], node_ids=[0, 1, 2],
        posteriors={0: A, 1: V_here, 2: V_here},
        pair_details=[
            {"pair_id": 1, "node_id": 0, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 0, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 10 — sibling merge: MRCAs are siblings, private segment empty -> 1.0.
    out.append(dict(
        name="sibling_merge",
        doc="pair MRCAs 3 and 4 share parent 10; private segment empty",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 4, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 11 — no changed pairs at all (all sides conserved) -> early zero return.
    out.append(dict(
        name="no_changed_pairs",
        doc="every pair holds ancestral on both sides",
        edges=e2, node_ids=n2, posteriors=base_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 12 — soft posteriors (real ASR uncertainty) so the axes land off 0/1.
    soft_post = {
        0: A_soft, 1: A_soft, 2: A_soft, 10: A_soft, 20: A_soft,
        3: {"V": 0.6, "A": 0.4}, 4: V_here,
        5: {"V": 0.55, "A": 0.45}, 6: V_here,
    }
    out.append(dict(
        name="soft_posteriors_midrange",
        doc="uncertain ASR; independence/core land strictly inside (0,1)",
        edges=e2, node_ids=n2, posteriors=soft_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=False, conserved_pair="",
    ))

    # 13 — T3a: both-sides split into two rows + a conserved pair in D_s of both.
    # Literal fixture from docs/scoring_v2_T3_core_pareado.md §10.1. Exercises
    # native_side_split=True: top and bottom are independent rows (no
    # 1-(1-t)(1-b)), the conserved pair counts in |D_s| on BOTH sides, and
    # conserved_pair_scores is still emitted without multiplying anything.
    #   TOP:    D={P1,P2,P4} n=3 -> core_top    = 2*0.77378/3 = 0.51585
    #   BOTTOM: D={P1,P3,P4} n=3 -> core_bottom = 2*0.81451/3 = 0.54301
    bs_edges = [(0, 1), (0, 2), (1, 10), (10, 3), (10, 4), (2, 20), (20, 5),
                (20, 6), (2, 21), (21, 9), (21, 12), (0, 30), (30, 7), (30, 8)]
    bs_nodes = [0, 1, 2, 10, 20, 21, 30, 3, 4, 5, 6, 9, 12, 7, 8]
    bs_post = {n: {"A": 0.90, "V": 0.05, "L": 0.05}
               for n in [0, 1, 2, 10, 20, 21, 30]}
    out.append(dict(
        name="both_sides_two_rows",
        doc="P1 changes both sides to different residues, P2 top-only, P3 "
            "bottom-only, P4 conserved; native_side_split -> two rows (§10.1)",
        edges=bs_edges, node_ids=bs_nodes, posteriors=bs_post,
        pair_details=[
            {"pair_id": 1, "node_id": 3, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "L"},
            {"pair_id": 2, "node_id": 5, "focal_state": "A",
             "top_tip_mode": "V", "bottom_tip_mode": "A"},
            {"pair_id": 3, "node_id": 7, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "L"},
            {"pair_id": 4, "node_id": 9, "focal_state": "A",
             "top_tip_mode": "A", "bottom_tip_mode": "A"},
        ],
        scheme="US", is_conserved_meta=True, conserved_pair="4",
        native_side_split=True,
    ))

    return out


def fop_scenarios():
    """Shared Python/R pooling oracle. Each entry: inputs + frozen pooled dict."""
    S = []

    S.append(dict(
        name="single_hypothesis_passthrough",
        hyp_records=[{"hyp": "H1", "asr_path_score": 0.42, "core": 0.6,
                      "independence": 0.9, "mrca_diversity": 0.5,
                      "derived_agreement": 1.0, "conservation_gate": 1.0,
                      "pair_scores": {"1": 0.8, "2": 0.7}}],
        pss=[], scheme="US",
    ))

    S.append(dict(
        name="two_hyp_two_domain_pss",
        hyp_records=[
            {"hyp": "H1", "asr_path_score": 0.70, "independence": 1.0,
             "mrca_diversity": 1.0, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.72,
             "pair_scores": {"1": 0.9, "2": 0.8}},
            {"hyp": "H2", "asr_path_score": 0.20, "independence": 1.0,
             "mrca_diversity": 0.5, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.24,
             "pair_scores": {"1": 0.2, "2": 0.8}},
        ],
        pss=[["H1", 1, 10.0], ["H1", 2, 1.0], ["H2", 1, 2.0], ["H2", 2, 2.0]],
        scheme="US",
    ))

    S.append(dict(
        name="two_hyp_equal_weight",
        hyp_records=S[-1]["hyp_records"],
        pss=[], scheme="US",
    ))

    S.append(dict(
        name="directional_core_opposite_sides",
        hyp_records=[
            {"hyp": "H1", "asr_path_score": 0.5, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.5,
             "pair_scores": {"1": 0.85, "2": 0.80},
             "pair_top_scores": {"1": 0.85}, "pair_bottom_scores": {"2": 0.80}},
            {"hyp": "H2", "asr_path_score": 0.5, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.5,
             "pair_scores": {"1": 0.85, "2": 0.80},
             "pair_top_scores": {"1": 0.85}, "pair_bottom_scores": {"2": 0.80}},
        ],
        pss=[], scheme="US",
    ))

    S.append(dict(
        name="conserved_gate_dedup",
        hyp_records=[
            {"hyp": "H1", "asr_path_score": 0.5, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 0.8, "core": 0.5,
             "pair_scores": {"1": 0.7, "2": 0.6},
             "conserved_pair_scores": {"5": 0.6},
             "conserved_pair_nodes": {"5": "cn_shared"}},
            {"hyp": "H2", "asr_path_score": 0.5, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 0.8, "core": 0.5,
             "pair_scores": {"1": 0.7, "2": 0.6},
             "conserved_pair_scores": {"5": 0.6},
             "conserved_pair_nodes": {"5": "cn_shared"}},
            {"hyp": "H3", "asr_path_score": 0.4, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 0.7, "core": 0.4,
             "pair_scores": {"1": 0.6, "2": 0.5},
             "conserved_pair_scores": {"5": 0.4},
             "conserved_pair_nodes": {"5": "cn_uniq"}},
        ],
        pss=[], scheme="US",
    ))

    S.append(dict(
        name="harvest_wide_da_us",
        hyp_records=[
            {"hyp": "H1", "asr_path_score": 0.5, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.5,
             "pair_scores": {"1": 0.7, "2": 0.6},
             "pair_derived_top": {}, "pair_derived_bot": {"1": "I", "2": "V"}},
            {"hyp": "H2", "asr_path_score": 0.4, "independence": 1.0,
             "mrca_diversity": 0.0, "derived_agreement": 1.0,
             "conservation_gate": 1.0, "core": 0.4,
             "pair_scores": {"1": 0.7, "2": 0.6},
             "pair_derived_top": {}, "pair_derived_bot": {"1": "V", "2": "I"}},
        ],
        pss=[], scheme="US",
    ))

    S.append(dict(
        name="harvest_wide_da_gs4",
        hyp_records=S[-1]["hyp_records"],
        pss=[], scheme="GS4",
    ))

    return S


def _pss_map(rows):
    return {(h, int(d)): float(v) for h, d, v in rows} or None


def _keys_to_int(d):
    return {int(k): v for k, v in d.items()}


def main():
    # ── path_scores golden ──────────────────────────────────────────────────
    ps_out = []
    for sc in scenarios():
        result = run_scenario(scenario_to_json(sc), ps)
        ps_out.append({"scenario": scenario_to_json(sc), "expected": result})
    (FIXTURE_DIR / "path_scores_golden.json").write_text(
        json.dumps(ps_out, indent=2, sort_keys=True) + "\n"
    )
    print(f"path_scores_golden.json: {len(ps_out)} scenarios")

    # ── fop_pool shared fixture ─────────────────────────────────────────────
    fp_out = []
    for sc in fop_scenarios():
        recs = [
            {**r,
             **{k: _keys_to_int(r[k]) for k in
                ("pair_scores", "pair_top_scores", "pair_bottom_scores",
                 "pair_derived_top", "pair_derived_bot",
                 "conserved_pair_scores", "conserved_pair_nodes") if k in r}}
            for r in sc["hyp_records"]
        ]
        expected = fop.pool_hypotheses(recs, _pss_map(sc["pss"]), scheme=sc["scheme"])
        fp_out.append({
            "name": sc["name"],
            "hyp_records": sc["hyp_records"],
            "pss": sc["pss"],
            "scheme": sc["scheme"],
            "expected": expected,
        })
    (FIXTURE_DIR / "fop_pool_fixture.json").write_text(
        json.dumps(fp_out, indent=2, sort_keys=True) + "\n"
    )
    print(f"fop_pool_fixture.json: {len(fp_out)} scenarios")


if __name__ == "__main__":
    main()
