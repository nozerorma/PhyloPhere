#!/usr/bin/env python3
# test_species_mapping.py — Unit tests of the tree/alignment matching by tax_id.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/phylo/
#
# Run from subworkflows/CT_DISAMBIGUATION/local: python3 -m pytest src/phylo/test_species_mapping.py

from io import StringIO

import pytest
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.phylo.species_mapping import match_tree_alignment_by_taxid


def _inputs(species):
    tree = Phylo.read(StringIO("(" + ",".join(f"{s}:1" for s in species) + ");"), "newick")
    aln = MultipleSeqAlignment([SeqRecord(Seq("ACGT"), id=s, description="") for s in species])
    return tree, aln


def test_unique_tax_ids_relabel_tree_and_alignment():
    tree, aln = _inputs(["A_sp", "B_sp", "C_sp"])
    out = match_tree_alignment_by_taxid(tree, aln, {"A_sp": "1", "B_sp": "2", "C_sp": "3"})
    assert len(out) == 4
    _, matched_aln, tree_map, aln_map = out
    assert sorted(r.id for r in matched_aln) == ["1", "2", "3"]
    assert tree_map == aln_map == {"1": "A_sp", "2": "B_sp", "3": "C_sp"}


def test_names_missing_from_the_map_are_dropped():
    tree, aln = _inputs(["A_sp", "B_sp", "Unmapped"])
    _, matched_aln, _, aln_map = match_tree_alignment_by_taxid(tree, aln, {"A_sp": "1", "B_sp": "2"})
    assert sorted(aln_map) == ["1", "2"]


def test_shared_tax_id_raises_and_names_the_curated_map():
    tree, aln = _inputs(["A_sp", "B_sp"])
    with pytest.raises(ValueError, match="species_taxid_map.tsv"):
        match_tree_alignment_by_taxid(tree, aln, {"A_sp": "1", "B_sp": "1"})
