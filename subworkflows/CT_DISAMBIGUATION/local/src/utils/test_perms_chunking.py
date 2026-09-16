#!/usr/bin/env python3
"""Regression net for Stage 2 of docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md:
splitting a gene's null-replay across multiple worker chunks
(_perms_worker_replay + _perms_worker_finalize) instead of one task per gene
(the old monolithic _perms_worker). The critical correctness property is that
chunking must never change output: a chunk boundary must never split one base
cycle's "<base>~H*" variants (_chunk_gene_cycles), and _perms_worker_finalize's
whole-gene reduction (n_detected, clustering, detail rows) must be insensitive
to how its merged input was partitioned across replay chunks.

Run:  python -m pytest test_perms_chunking.py
"""
from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # .../local
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import src.utils.gene_wrapper as gw  # noqa: E402

_chunk_gene_cycles = gw._chunk_gene_cycles
_perms_worker = gw._perms_worker
_perms_worker_finalize = gw._perms_worker_finalize
_perms_worker_replay = gw._perms_worker_replay


# ── _chunk_gene_cycles ────────────────────────────────────────────────────────

def test_chunk_never_splits_a_base_cycles_hypothesis_variants():
    tags = ["b_1~H1", "b_1~H2", "b_1~H3", "b_2~H1", "b_2~H2", "b_3~H1"]
    chunks = _chunk_gene_cycles(tags, target_chunk_size=2)
    seen = set()
    for chunk in chunks:
        bases_here = {t.split("~", 1)[0] for t in chunk}
        # every base cycle appearing in this chunk contributes ALL its tags here
        for b in bases_here:
            assert b not in seen, f"base cycle {b} split across chunks"
            seen.add(b)
    # every tag accounted for exactly once
    assert sorted(t for c in chunks for t in c) == sorted(tags)


def test_chunk_respects_target_size_when_base_cycles_are_singletons():
    tags = [f"b_{i}" for i in range(10)]
    chunks = _chunk_gene_cycles(tags, target_chunk_size=3)
    assert [len(c) for c in chunks] == [3, 3, 3, 1]


def test_chunk_target_larger_than_input_returns_one_chunk():
    tags = ["b_1~H1", "b_1~H2", "b_2~H1"]
    chunks = _chunk_gene_cycles(tags, target_chunk_size=1000)
    assert chunks == [tags]


def test_chunk_empty_input():
    assert _chunk_gene_cycles([], target_chunk_size=10) == [[]]


# ── _perms_worker_finalize: partition-invariant whole-gene reduction ──────────

def _rec(position, side, score, group="US", cyc=None):
    return types.SimpleNamespace(
        position=position, caap_group=group, side=side,
        asr_path_score=score, cycle=cyc,
    )


def test_finalize_is_invariant_to_how_chunks_are_partitioned():
    # 4 base cycles' worth of pooled records, as _perms_worker_replay would
    # return them (already FOP-pooled: one record per (base_cycle, position, side)).
    full = [
        ("b_1", [_rec(10, "top", 0.5), _rec(20, "bottom", 0.2)]),
        ("b_2", [_rec(10, "top", 0.6)]),
        ("b_3", [_rec(20, "bottom", 0.1), _rec(30, "top", 0.9)]),
        ("b_4", [_rec(10, "top", 0.4)]),
    ]
    n_cycles_total = 4

    gene_one_chunk, detail_one, pval_one = _perms_worker_finalize(
        "GENE", full, n_cycles_total,
    )

    # Same data, delivered as if two replay chunks (b_1,b_2 then b_3,b_4) had
    # been merged in the parent before finalize -- exactly what
    # process_all_genes_perms' pending_pooled accumulation does.
    chunk_a = full[:2]
    chunk_b = full[2:]
    merged = list(chunk_a) + list(chunk_b)
    gene_two_chunks, detail_two, pval_two = _perms_worker_finalize(
        "GENE", merged, n_cycles_total,
    )

    assert gene_one_chunk == gene_two_chunks == "GENE"
    assert sorted(detail_one, key=lambda r: (r["cycle"], r["Position"], r["side"])) == \
        sorted(detail_two, key=lambda r: (r["cycle"], r["Position"], r["side"]))
    assert sorted(pval_one, key=lambda r: (r["Position"], r["caap_group"])) == \
        sorted(pval_two, key=lambda r: (r["Position"], r["caap_group"]))

    # And insensitive to chunk arrival ORDER (imap_unordered gives no guarantee).
    reordered = list(chunk_b) + list(chunk_a)
    _, detail_reordered, pval_reordered = _perms_worker_finalize(
        "GENE", reordered, n_cycles_total,
    )
    assert sorted(detail_reordered, key=lambda r: (r["cycle"], r["Position"], r["side"])) == \
        sorted(detail_one, key=lambda r: (r["cycle"], r["Position"], r["side"]))
    assert sorted(pval_reordered, key=lambda r: (r["Position"], r["caap_group"])) == \
        sorted(pval_one, key=lambda r: (r["Position"], r["caap_group"]))


# ── End-to-end: unchunked _perms_worker vs. replay-in-N-chunks + finalize ─────
# Drives the REAL _perms_worker_replay (including its FOP domain-pooling), with
# _load_gene_asr_context and analyze_gene_disambiguation faked out so the test
# doesn't need a real alignment/tree/ASR fixture -- only the chunking/merge
# machinery around them is under test here.

def _fake_ctx():
    return {
        "alignment_data": types.SimpleNamespace(species_to_taxid={}),
        "tree_data": object(),
        "node_posteriors": types.SimpleNamespace(posteriors_node={}),
    }


def _fake_analyze(gene, caas_entries, trait_pairs, **kwargs):
    """One record per caas_entry position, `sides` scored deterministically from
    (position, fg species count) so different cycles/hypotheses actually differ
    and the FOP pooler has something non-trivial to average."""
    fg_n = len(trait_pairs[1])
    out = []
    for entry in caas_entries:
        pos = entry.position
        score = (pos % 7) / 10.0 + fg_n * 0.01
        sides = {
            "domain_meta": {0: "d0"},
            "top": {
                "domain_scores": {0: score},
                "domain_der_enc": {0: "A"},
                "domain_der": {0: "A"},
                "domain_anc": {0: "B"},
            },
            "bottom": {"domain_scores": {0: 0.0}},
        }
        out.append(types.SimpleNamespace(
            position=pos, caap_group="US", side=None,
            hypothesis=None, sides=sides,
        ))
    return out, None


def _make_discovery_file(tmp_path, gene, cycle_tags):
    """One perm_discovery row per (cycle tag, position) so every replayed
    labeling has >=1 CAAS entry to replay -- _perms_worker_replay skips cycles
    with no discovery hits. Columns/casing match _parse_discovery_entries'
    expectations exactly (lowercase "cycle"/"position", "gene" for the
    gene_filter path used when perm_discovery_file is a single file)."""
    disc = tmp_path / "perm_disc.tsv"
    with open(disc, "w") as f:
        f.write("gene\tcycle\tposition\n")
        for cyc in cycle_tags:
            for pos in (10, 20):
                f.write(f"{gene}\t{cyc}\t{pos}\n")
    return disc


def _cycle_labelings(base_cycles, n_hyp):
    labelings = {}
    fg = ["SpA", "SpB"]
    bg = ["SpC", "SpD"]
    for b in base_cycles:
        for h in range(1, n_hyp + 1):
            labelings[f"{b}~H{h}"] = (fg, bg)
    return labelings


def test_chunked_replay_matches_unchunked_end_to_end():
    import tempfile

    base_cycles = ["b_1", "b_2", "b_3", "b_4"]
    n_hyp = 3
    cycle_tags = [f"{b}~H{h}" for b in base_cycles for h in range(1, n_hyp + 1)]
    cycle_labelings = _cycle_labelings(base_cycles, n_hyp)
    fop_pairs = {b: {} for b in base_cycles}  # empty pss map -> equal weight

    with tempfile.TemporaryDirectory() as td:
        tmp_path = Path(td)
        disc_path = _make_discovery_file(tmp_path, "GENE1", cycle_tags)

        with patch("src.utils.gene_wrapper._load_gene_asr_context", return_value=_fake_ctx()), \
             patch("src.utils.gene_wrapper.analyze_gene_disambiguation", side_effect=_fake_analyze):

            # Unchunked: the pre-Stage-2 single-task-per-gene call.
            gene_u, detail_u, pval_u = _perms_worker(
                "GENE1", "unused_align_dir", "unused_tree", None,
                "lg", "unused_asr_cache", 0.1,
                cycle_tags, cycle_labelings, str(disc_path),
                None, fop_pairs,
            )

            # Chunked: force 3 chunks (target_chunk_size smaller than one base
            # cycle's hyp-variant count would risk splitting it -- exercise that
            # exact edge by asking for a tiny target so _chunk_gene_cycles must
            # fall back to "whole base cycle" chunks).
            chunks = _chunk_gene_cycles(cycle_tags, target_chunk_size=2)
            assert len(chunks) > 1, "test is only meaningful if it actually chunks"

            pooled_all = []
            for chunk in chunks:
                _, pooled = _perms_worker_replay(
                    "GENE1", "unused_align_dir", "unused_tree", None,
                    "lg", "unused_asr_cache", 0.1,
                    chunk, cycle_labelings, str(disc_path),
                    None, fop_pairs,
                )
                pooled_all.extend(pooled)

            n_cycles_total = len(base_cycles)
            gene_c, detail_c, pval_c = _perms_worker_finalize(
                "GENE1", pooled_all, n_cycles_total,
            )

    assert gene_u == gene_c == "GENE1"
    key_d = lambda r: (r["cycle"], r["Position"], r["caap_group"], r["side"])
    key_p = lambda r: (r["Position"], r["caap_group"])
    assert sorted(detail_u, key=key_d) == sorted(detail_c, key=key_d)
    assert sorted(pval_u, key=key_p) == sorted(pval_c, key=key_p)
    assert len(detail_u) > 0  # sanity: the fixture actually produced rows
