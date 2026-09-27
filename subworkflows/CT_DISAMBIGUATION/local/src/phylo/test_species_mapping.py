"""Regression tests for match_tree_alignment_by_taxid with shared tax_ids.

Several accessions of one species share an NCBI tax_id. The alphabetically
first accession keeps the tax_id, the rest get synthetic ones, and every
accession must stay resolvable through the inverted species -> tax_id map that
asr_single.py builds, independently of PYTHONHASHSEED.
"""
import json
import os
import subprocess
import sys
from io import StringIO
from pathlib import Path

from Bio import AlignIO, Phylo

LOCAL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(LOCAL))

from src.phylo.species_mapping import match_tree_alignment_by_taxid  # noqa: E402

TREE = "((Sp_a_1:1,Sp_a_2:1):1,(Sp_a_3:1,(Sp_b:1,Sp_c:1):1):1);"
FASTA = ">Sp_a_1\nAAAA\n>Sp_a_2\nAAAC\n>Sp_a_3\nAACC\n>Sp_b\nACCC\n>Sp_c\nCCCC\n"
TAXIDS = {"Sp_a_1": "100", "Sp_a_2": "100", "Sp_a_3": "100", "Sp_b": "200", "Sp_c": "300"}


def _run():
    tree = Phylo.read(StringIO(TREE), "newick")
    aln = AlignIO.read(StringIO(FASTA), "fasta")
    _, _, tree_t2s, aln_t2s, synth = match_tree_alignment_by_taxid(tree, aln, dict(TAXIDS))
    return tree_t2s, aln_t2s, synth


def test_kept_accession_owns_shared_taxid():
    tree_t2s, aln_t2s, synth = _run()
    assert aln_t2s["100"] == "Sp_a_1"
    assert tree_t2s["100"] == "Sp_a_1"
    assert set(synth) == {"Sp_a_2", "Sp_a_3"}


def test_every_accession_resolvable():
    _, aln_t2s, _ = _run()
    species_to_taxid = {sp: t for t, sp in aln_t2s.items()}
    assert set(species_to_taxid) == set(TAXIDS)
    assert len(set(species_to_taxid.values())) == len(TAXIDS)


def test_independent_of_hash_seed():
    code = (
        "import json, sys; sys.path.insert(0, %r); "
        "from src.phylo.test_species_mapping import _run; "
        "t, a, s = _run(); print(json.dumps([t, a, sorted(s)], sort_keys=True))"
    ) % str(LOCAL)
    outs = set()
    for seed in ("0", "1", "4", "6", "99"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        res = subprocess.run([sys.executable, "-c", code], env=env, cwd=str(LOCAL),
                             capture_output=True, text=True, check=True)
        outs.add(res.stdout.strip().splitlines()[-1])
    assert len(outs) == 1, outs
