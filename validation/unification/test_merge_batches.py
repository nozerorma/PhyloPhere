"""merge_disambiguation_batches.py: the merged outputs do not depend on the order the batches finished."""
import csv
import itertools
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local/scripts/merge_disambiguation_batches.py"
HEADER = ["gene", "msa_pos", "side", "asr_path_score"]


def _batch(root, name, rows, asr_rows=(), header=HEADER):
    """A batch ct_disambiguation/ directory; rows are (gene, msa_pos, side, score) in production order."""
    d = root / name
    (d / "diagnostics").mkdir(parents=True)
    for rel in ("caas_convergence_master.csv", "diagnostics/no_change_debug.csv"):
        with open(d / rel, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(r for r in rows if rel.startswith("caas") or r[2] == "none")
    with open(d / "diagnostics/caas_hypothesis_domain_asr.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "position", "hypothesis"])
        w.writerows(asr_rows)
    genes = {}
    for r in rows:
        genes[r[0]] = genes.get(r[0], 0) + 1
    (d / "caas_convergence_summary.json").write_text(json.dumps(
        {"metadata": {"num_genes": len(genes), "total_positions": len(rows), "num_pairs": 4, "schema_version": "x"},
         "by_gene_counts": genes}))
    (d / "aggregation.sqlite3").write_text(name)
    return d


def _merge(dirs, out):
    p = subprocess.run([sys.executable, str(SCRIPT), "--batch-dirs", *map(str, dirs), "--output-dir", str(out)],
                       capture_output=True, text=True)
    return p


def _read(path, delimiter=","):
    return list(csv.reader(open(path, newline=""), delimiter=delimiter))


# gene rows in production order; positions 9 and 100 check numeric (not lexical) ordering
B1 = [("AIPL1", 100, "top", 0.5), ("AIPL1", 9, "bottom", 0.25), ("AIPL1", 9, "top", 0.75), ("ANK3", 5, "none", 0.0)]
B2 = [("SPN", 3, "top", 0.1), ("SHANK2", 7, "none", 0.0), ("SHANK2", 7, "top", 0.2)]
ASR1 = [("AIPL1", 9, "H1"), ("AIPL1", 100, "H1")]
ASR2 = [("SPN", 3, "H1"), ("SHANK2", 7, "H2")]


def test_merged_outputs_are_independent_of_batch_arrival_order(tmp_path):
    b1 = _batch(tmp_path, "batch_1", B1, ASR1)
    b2 = _batch(tmp_path, "batch_2", B2, ASR2)
    outs = {}
    for i, order in enumerate(itertools.permutations([b1, b2])):
        out = tmp_path / f"out{i}"
        assert _merge(order, out).returncode == 0
        outs[i] = {rel: (out / rel).read_bytes() for rel in
                   ("caas_convergence_master.csv", "diagnostics/no_change_debug.csv",
                    "diagnostics/caas_hypothesis_domain_asr.tsv", "caas_convergence_summary.json",
                    "aggregation_batch_00001.sqlite3", "aggregation_batch_00002.sqlite3")}
    assert outs[0] == outs[1]
    assert outs[0]["aggregation_batch_00001.sqlite3"] == b"batch_1"  # numbered by gene order, not by arrival


def test_merge_equals_the_unbatched_master(tmp_path):
    """One directory with every gene (an unbatched run: gene, msa_pos, production order) is the reference."""
    everything = _batch(tmp_path, "unbatched", B1 + B2, ASR1 + ASR2)
    ref = _merge([everything], tmp_path / "ref")
    b2, b1 = _batch(tmp_path, "batch_2", B2, ASR2), _batch(tmp_path, "batch_1", B1, ASR1)  # finished in reverse
    assert ref.returncode == 0 and _merge([b2, b1], tmp_path / "merged").returncode == 0
    merged = _read(tmp_path / "merged/caas_convergence_master.csv")
    assert merged == _read(tmp_path / "ref/caas_convergence_master.csv")
    assert [(r[0], r[1], r[2]) for r in merged[1:]] == [
        ("AIPL1", "9", "bottom"), ("AIPL1", "9", "top"), ("AIPL1", "100", "top"), ("ANK3", "5", "none"),
        ("SHANK2", "7", "none"), ("SHANK2", "7", "top"), ("SPN", "3", "top")]
    assert _read(tmp_path / "merged/diagnostics/no_change_debug.csv")[1:] == [["ANK3", "5", "none", "0.0"], ["SHANK2", "7", "none", "0.0"]]
    asr = _read(tmp_path / "merged/diagnostics/caas_hypothesis_domain_asr.tsv", "\t")
    assert [r[0] for r in asr[1:]] == ["AIPL1", "AIPL1", "SHANK2", "SPN"] and asr[1][1] == "9"  # a gene keeps its own order


def test_schema_mismatch_is_still_an_error(tmp_path):
    b1 = _batch(tmp_path, "batch_1", B1)
    b2 = _batch(tmp_path, "batch_2", B2, header=HEADER + ["domain_5_score"])
    p = _merge([b1, b2], tmp_path / "out")
    assert p.returncode != 0 and "Schema mismatch" in p.stderr
