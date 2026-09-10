#!/usr/bin/env python3
"""Tests for residue_descriptors.add_residue_descriptors.

Layout convention: ``derived_residues`` is ``<top>/<bottom>`` (same left/right as
the ``caas`` string). The ``side``-sanctioned side shows derived residues;
the other side shows the ancestral residue.

Run: ``python3 -m pytest subworkflows/CT_POSTPROC/local/src/test_residue_descriptors.py``
"""

import pandas as pd

from residue_descriptors import (
    DESCRIPTOR_COLUMNS,
    SPECIES_TALLY_COLUMNS,
    add_residue_descriptors,
    add_species_tally,
)


def _mk_p3() -> pd.DataFrame:
    # test_fop_pool.R mk_p3: 3 hypothesis rows, one (Gene, Position), bottom-side
    # change (bot_aa I/V), no side column -> inferred from the data.
    return pd.DataFrame(
        {
            "Gene": ["BRCA1", "BRCA1", "BRCA1"],
            "Position": [96, 96, 96],
            "mrca_1_node": ["p1", "p3", "p1"],
            "mrca_2_node": ["p2", "p4", "p2"],
            "mrca_1_anc_aa": ["A", "A", "A"],
            "mrca_2_anc_aa": ["A", "A", "A"],
            "mrca_1_top_aa": ["", "", ""],
            "mrca_2_top_aa": ["", "", ""],
            "mrca_1_bot_aa": ["I", "I", "I"],
            "mrca_2_bot_aa": ["V", "V", "V"],
        }
    )


def test_point3_top_is_ancestral_bottom_is_derived():
    out = add_residue_descriptors(_mk_p3())
    row = out.iloc[0]
    assert row["derived_residues"] == "A/IV"          # top=anc A, bottom=derived IV
    # actual support = distinct CAAS pairs (2 pairs; pair1->I, pair2->V, both anc A)
    assert row["top_residue_support"] == "A:2"
    assert row["bottom_residue_support"] == "I:1,V:1"
    # _detail = distinct reconstructed nodes (pair1 nodes {p1,p3}, pair2 {p2,p4})
    assert row["top_residue_support_detail"] == "A:4"
    assert row["bottom_residue_support_detail"] == "I:2,V:2"
    assert row["n_conserved_pairs"] == ""             # no conserved block
    assert (out["derived_residues"] == "A/IV").all()  # broadcast


def test_side_top_puts_derived_left_ancestral_right():
    # PEPC:539-style — side == "top" everywhere; one row also has a stray
    # bottom residue (n90) which must NOT surface (side sanctions top only).
    df = pd.DataFrame(
        {
            "Gene": ["P", "P", "P"],
            "Position": [539, 539, 539],
            "side": ["top", "top", "top"],
            "mrca_1_node": ["n90", "n91", "n93"],
            "mrca_1_anc_aa": ["P", "P", "P"],
            "mrca_1_top_aa": ["T", "T", "T"],
            "mrca_1_bot_aa": ["S", "", ""],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "T/P"       # top derived T, bottom anc P
    assert out.iloc[0]["top_residue_support"] == "T:1"    # 1 physical pair
    assert out.iloc[0]["bottom_residue_support"] == "P:1"
    assert out.iloc[0]["top_residue_support_detail"] == "T:3"     # 3 hypothesis reconstructions
    assert out.iloc[0]["bottom_residue_support_detail"] == "P:3"


def test_side_bottom_puts_derived_right():
    df = pd.DataFrame(
        {
            "Gene": ["G"],
            "Position": [1],
            "side": ["bottom"],
            "mrca_1_node": ["a"],
            "mrca_1_anc_aa": ["A"],
            "mrca_1_top_aa": ["W"],   # ignored: side is bottom
            "mrca_1_bot_aa": ["C"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "A/C"       # top anc A, bottom derived C
    assert out.iloc[0]["top_residue_support"] == "A:1"
    assert out.iloc[0]["bottom_residue_support"] == "C:1"
    assert out.iloc[0]["top_residue_support_detail"] == "A:1"
    assert out.iloc[0]["bottom_residue_support_detail"] == "C:1"


def test_side_none_infers_both_from_data():
    # A row whose `side` is unusable ("none") falls back to inferring derived
    # sides from which clades actually substituted (here: both).
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "side": ["none", "none"],
            "mrca_1_node": ["a", "b"],
            "mrca_1_anc_aa": ["M", "M"],
            "mrca_1_top_aa": ["L", "L"],
            "mrca_1_bot_aa": ["F", "F"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "L/F"
    assert out.iloc[0]["top_residue_support"] == "L:1"           # 1 physical pair
    assert out.iloc[0]["bottom_residue_support"] == "F:1"
    assert out.iloc[0]["top_residue_support_detail"] == "L:2"    # nodes a, b
    assert out.iloc[0]["bottom_residue_support_detail"] == "F:2"


def test_multi_residue_side_sorted():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "side": ["top", "top"],
            "mrca_1_node": ["a", "b"],
            "mrca_1_anc_aa": ["I", "I"],
            "mrca_1_top_aa": ["S", "M"],
            "mrca_1_bot_aa": ["", ""],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "MS/I"
    # one pair, but hypotheses disagree S vs M -> the pair counts once for each
    assert out.iloc[0]["top_residue_support"] == "M:1,S:1"
    assert out.iloc[0]["top_residue_support_detail"] == "M:1,S:1"


def test_n_conserved_pairs_counts_distinct_nodes():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "side": ["top", "top"],
            "mrca_1_node": ["a", "a"],
            "mrca_1_anc_aa": ["I", "I"],
            "mrca_1_top_aa": ["F", "F"],
            "mrca_1_bot_aa": ["", ""],
            "conserved_1_node": ["c1", "c2"],
            "conserved_1_cons": ["0.9", "0.7"],
            "conserved_2_node": ["c1", ""],
            "conserved_2_cons": ["0.8", ""],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "F/I"
    assert out.iloc[0]["n_conserved_pairs"] == "2"   # c1, c2


def test_no_changed_pairs_all_empty_but_conserved_counted():
    df = pd.DataFrame(
        {
            "Gene": ["G"],
            "Position": [1],
            "side": ["top"],
            "mrca_1_node": ["a"],
            "mrca_1_anc_aa": ["A"],
            "mrca_1_top_aa": [""],
            "mrca_1_bot_aa": [""],
            "conserved_1_node": ["c1"],
            "conserved_1_cons": ["0.9"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == ""
    assert out.iloc[0]["top_residue_support"] == ""
    assert out.iloc[0]["bottom_residue_support"] == ""
    assert out.iloc[0]["n_conserved_pairs"] == "1"


def test_no_raw_block_stable_schema():
    df = pd.DataFrame({"Gene": ["G"], "Position": [1], "caas": ["A/B"]})
    out = add_residue_descriptors(df)
    for c in DESCRIPTOR_COLUMNS:
        assert c in out.columns
        assert out.iloc[0][c] == ""


def test_pair_vs_node_support_diverge():
    # One physical pair whose 3 hypothesis rows resolve to 2 distinct nodes.
    df = pd.DataFrame(
        {
            "Gene": ["G", "G", "G"],
            "Position": [1, 1, 1],
            "side": ["bottom", "bottom", "bottom"],
            "mrca_1_node": ["a", "a", "b"],
            "mrca_1_anc_aa": ["A", "A", "A"],
            "mrca_1_top_aa": ["", "", ""],
            "mrca_1_bot_aa": ["L", "L", "L"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["bottom_residue_support"] == "L:1"          # 1 physical pair
    assert out.iloc[0]["bottom_residue_support_detail"] == "L:2"   # nodes a, b
    assert out.iloc[0]["derived_residues"] == "A/L"


# ── add_species_tally ───────────────────────────────────────────────────────

def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def test_species_tally_counts_contrast_species(tmp_path):
    aln_dir = tmp_path / "alignments"
    aln_dir.mkdir()
    # column index 2 (0-based, = Position value): fg mostly N, bg all A
    _write_fasta(aln_dir / "GENE.filtered.fa", {
        "sp_fg1": "AANXX", "sp_fg2": "AANXX", "sp_fg3": "AAAXX",  # 2 N, 1 A
        "sp_bg1": "AAAXX", "sp_bg2": "AAAXX",                     # 2 A
        "sp_bg3": "AA-XX",                                        # gap -> not counted
        "outgroup": "AAKXX",                                     # not in either list
    })
    df = pd.DataFrame({"Gene": ["GENE", "GENE"], "Position": [2, 2]})
    out = add_species_tally(
        df, str(aln_dir),
        ["sp_fg1", "sp_fg2", "sp_fg3"],
        ["sp_bg1", "sp_bg2", "sp_bg3"],
    )
    assert out.iloc[0]["top_species_residues"] == "N:2,A:1"
    assert out.iloc[0]["bottom_species_residues"] == "A:2"
    assert out.iloc[0]["n_top_species"] == "3"
    assert out.iloc[0]["n_bottom_species"] == "2"          # bg3 is a gap
    assert (out["top_species_residues"] == "N:2,A:1").all()  # broadcast


def test_species_tally_noop_without_alignment(tmp_path):
    df = pd.DataFrame({"Gene": ["G"], "Position": [1]})
    out = add_species_tally(df, None, ["a"], ["b"])
    for c in SPECIES_TALLY_COLUMNS:
        assert c in out.columns and out.iloc[0][c] == ""


def test_species_tally_out_of_range_position(tmp_path):
    aln_dir = tmp_path / "alignments"
    aln_dir.mkdir()
    _write_fasta(aln_dir / "G.fa", {"a": "MMM", "b": "MMM"})
    df = pd.DataFrame({"Gene": ["G"], "Position": [99]})
    out = add_species_tally(df, str(aln_dir), ["a"], ["b"])
    assert out.iloc[0]["top_species_residues"] == ""
