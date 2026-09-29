"""CONCAT_DISCOVERY / CONCAT_BACKGROUND: the concatenated tables do not depend on the order the tasks finished.

The scripts are taken from the process definitions themselves and run with bash on staged files, the way
Nextflow leaves them (discovery_1, discovery_2, ... in arrival order).
"""
import itertools
import re
import subprocess
from pathlib import Path

import pytest

NF = Path(__file__).resolve().parents[2] / "subworkflows/CT/ct_concat.nf"
HEADER = "gene\tmode\tcaap_group\ttrait\tposition"


def _script(process):
    text = NF.read_text()
    body = text[text.index(f"process {process} "):]
    body = body[body.index('script:\n    """') + len('script:\n    """'):]
    body = body[:body.index('\n    """')]
    # Groovy string semantics for the two escapes the scripts use: \$ -> $ and \t -> tab
    return body.replace("\\$", "$").replace("\\t", "\t")


def _run(process, staged, tmp_path, prefix):
    work = tmp_path / f"w_{process}_{abs(hash(tuple(staged)))}"
    work.mkdir()
    for i, content in enumerate(staged, start=1):
        (work / f"{prefix}_{i}").write_text(content)
    (work / "run.sh").write_text(_script(process))
    p = subprocess.run(["bash", "run.sh"], cwd=work, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-800:] + p.stderr[-800:]
    return work


# one file per batch of genes; a gene's rows keep their production order; H10 before H2 within a gene
F1 = f"{HEADER}\nSPN\tCAAP\tUS\tH10\t9\nSPN\tCAAP\tUS\tH2\t9\nSPN\tCAAP\tGS1\tH10\t100\n"
F2 = f"{HEADER}\nAIPL1\tCAAP\tUS\tH1\t5\nAIPL1\tCAAP\tUS\tH1\t3\n"
F3 = f"{HEADER}\nSHANK2\tCAAP\tUS\tH1\t7\n"
EXPECTED = (HEADER + "\n"
            "AIPL1\tCAAP\tUS\tH1\t5\nAIPL1\tCAAP\tUS\tH1\t3\n"      # production order inside a gene is kept
            "SHANK2\tCAAP\tUS\tH1\t7\n"
            "SPN\tCAAP\tUS\tH10\t9\nSPN\tCAAP\tUS\tH2\t9\nSPN\tCAAP\tGS1\tH10\t100\n")


def test_discovery_is_gene_ordered_whatever_the_arrival_order(tmp_path):
    outs = set()
    for order in itertools.permutations([F1, F2, F3]):
        work = _run("CONCAT_DISCOVERY", order, tmp_path, "discovery")
        outs.add((work / "discovery.tab").read_text())
    assert outs == {EXPECTED}


BG = {"a": "SPN\t9,100\n", "b": "AIPL1\tNULL\n", "c": "SHANK2\t7\n", "d": "ANK3\t\n"}


def test_background_is_gene_ordered_and_lists_genes_with_positions(tmp_path):
    outs = set()
    for order in itertools.permutations(list(BG.values())):
        work = _run("CONCAT_BACKGROUND", order, tmp_path, "background")
        outs.add(((work / "background.output").read_text(), (work / "background_genes.output").read_text()))
    assert outs == {("AIPL1\tNULL\nANK3\t\nSHANK2\t7\nSPN\t9,100\n", "SHANK2\nSPN\n")}


def test_data_sorts_keep_their_temporary_files_out_of_tmp():
    """Cluster policy: no /tmp. Every sort of data (those with options) runs with -T . (the work dir)."""
    for process in ("CONCAT_DISCOVERY", "CONCAT_BACKGROUND"):
        data_sorts = [l for l in _script(process).splitlines() if re.search(r"\bsort -", l)]
        assert data_sorts and all("-T ." in l and "LC_ALL=C" in l for l in data_sorts), (process, data_sorts)
