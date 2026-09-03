#!/usr/bin/env python3
"""Tests for residue_descriptors.add_residue_descriptors.

Layout convention: ``derived_residues`` is ``<top>/<bottom>`` (same left/right as
the ``caas`` string). The ``change_side``-sanctioned side shows derived residues;
the other side shows the ancestral residue.

Run: ``python3 -m pytest subworkflows/CT_POSTPROC/local/src/test_residue_descriptors.py``
"""

import pandas as pd

from residue_descriptors import DESCRIPTOR_COLUMNS, add_residue_descriptors


def _mk_p3() -> pd.DataFrame:
    # test_fop_pool.R mk_p3: 3 hypothesis rows, one (Gene, Position), bottom-side
    # change (bot_aa I/V), no change_side column -> inferred from the data.
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
    assert row["top_residue_support"] == "A:4"
    assert row["bottom_residue_support"] == "I:2,V:2"
    assert row["n_conserved_pairs"] == ""             # no conserved block
    assert (out["derived_residues"] == "A/IV").all()  # broadcast


def test_change_side_top_puts_derived_left_ancestral_right():
    # PEPC:539-style — change_side == "top" everywhere; one row also has a stray
    # bottom residue (n90) which must NOT surface (change_side sanctions top only).
    df = pd.DataFrame(
        {
            "Gene": ["P", "P", "P"],
            "Position": [539, 539, 539],
            "change_side": ["top", "top", "top"],
            "mrca_1_node": ["n90", "n91", "n93"],
            "mrca_1_anc_aa": ["P", "P", "P"],
            "mrca_1_top_aa": ["T", "T", "T"],
            "mrca_1_bot_aa": ["S", "", ""],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "T/P"       # top derived T, bottom anc P
    assert out.iloc[0]["top_residue_support"] == "T:3"
    assert out.iloc[0]["bottom_residue_support"] == "P:3"


def test_change_side_bottom_puts_derived_right():
    df = pd.DataFrame(
        {
            "Gene": ["G"],
            "Position": [1],
            "change_side": ["bottom"],
            "mrca_1_node": ["a"],
            "mrca_1_anc_aa": ["A"],
            "mrca_1_top_aa": ["W"],   # ignored: change_side is bottom
            "mrca_1_bot_aa": ["C"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "A/C"       # top anc A, bottom derived C
    assert out.iloc[0]["top_residue_support"] == "A:1"
    assert out.iloc[0]["bottom_residue_support"] == "C:1"


def test_change_side_both_shows_derived_on_both():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "change_side": ["both", "both"],
            "mrca_1_node": ["a", "b"],
            "mrca_1_anc_aa": ["M", "M"],
            "mrca_1_top_aa": ["L", "L"],
            "mrca_1_bot_aa": ["F", "F"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "L/F"
    assert out.iloc[0]["top_residue_support"] == "L:2"
    assert out.iloc[0]["bottom_residue_support"] == "F:2"


def test_multi_residue_side_sorted():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "change_side": ["top", "top"],
            "mrca_1_node": ["a", "b"],
            "mrca_1_anc_aa": ["I", "I"],
            "mrca_1_top_aa": ["S", "M"],
            "mrca_1_bot_aa": ["", ""],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["derived_residues"] == "MS/I"
    assert out.iloc[0]["top_residue_support"] == "M:1,S:1"


def test_n_conserved_pairs_counts_distinct_nodes():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G"],
            "Position": [1, 1],
            "change_side": ["top", "top"],
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
            "change_side": ["top"],
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


def test_distinct_node_counting_dedups_repeats():
    df = pd.DataFrame(
        {
            "Gene": ["G", "G", "G"],
            "Position": [1, 1, 1],
            "change_side": ["bottom", "bottom", "bottom"],
            "mrca_1_node": ["a", "a", "b"],
            "mrca_1_anc_aa": ["A", "A", "A"],
            "mrca_1_top_aa": ["", "", ""],
            "mrca_1_bot_aa": ["L", "L", "L"],
        }
    )
    out = add_residue_descriptors(df)
    assert out.iloc[0]["bottom_residue_support"] == "L:2"  # nodes a, b
    assert out.iloc[0]["derived_residues"] == "A/L"
