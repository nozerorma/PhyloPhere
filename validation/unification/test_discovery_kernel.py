"""Vectorized perm-replay kernel vs the scalar CAAStools discovery, on a synthetic alignment.

`miss_pair` discards a labeling when the fg and bg thresholds are equal and the fg and bg sides
both have gapped (or missing) species, but in different pairs. The tests build alignments where
that rule changes the result and check that the scalar path and the kernel agree in both modes,
and that the rule is not vacuous (the scalar path itself differs with and without it).
"""
import subprocess
from pathlib import Path

import pandas as pd
import pytest

CT = Path(__file__).resolve().parents[2] / "subworkflows/CT/local/ct"

# three pairs: (f1,b1) (f2,b2) (f3,b3); fg residues D, bg residues K
TRAITS = "f1\t1\t1\nb1\t0\t1\nf2\t1\t2\nb2\t0\t2\nf3\t1\t3\nb3\t0\t3\n"
FG, BG = "f1,f2,f3", "b1,b2,b3"


def _fasta(path, seqs):
    path.write_text("".join(f">{s}\n{q}\n" for s, q in seqs.items()))


def _run(cmd, cwd):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr


def _discover(tmp, seqs, thresholds, miss_pair):
    """-> (scalar rows, kernel rows), each a set of (caap_group, position, caas)."""
    _fasta(tmp / "G1.fasta", seqs)
    (tmp / "traits.tab").write_text(TRAITS)
    (tmp / "b0.tab").write_text(f"b_0\t{FG}\t{BG}\n")
    flag = ["--miss_pair"] if miss_pair else []
    common = ["-a", "G1.fasta", "-t", "traits.tab", "--fmt", "fasta", "--patterns", "1,2,3", "--caap_mode",
              "--max_conserved", "1", *thresholds, *flag]
    _run([str(CT), "discovery", *common, "-o", "scalar.out", "--background_output", "scalar.bg"], tmp)
    _run([str(CT), "perm-replay", "-a", "G1.fasta", "-t", "traits.tab", "-s", "b0.tab", "-o", "kernel.out",
          "--fmt", "fasta", "--patterns", "1,2,3", "--caap_mode", "--max_conserved", "1", *thresholds, *flag,
          "--export_perm_discovery", "kernel.disc"], tmp)
    key = ["caap_group", "position", "caas"]
    # a run that keeps no position writes no output file: that is the empty set
    rows = lambda f: set(map(tuple, pd.read_csv(tmp / f, sep="\t")[key].astype(str).itertuples(index=False, name=None))) \
        if (tmp / f).is_file() else set()
    return rows("scalar.out"), rows("kernel.disc")


def _positions(rows):
    return {int(p) for _g, p, _c in rows}


# Column 0: f1 and b1 gapped (same pair)  -> kept
# Column 1: f1 (pair 1) and b2 (pair 2) gapped (different pairs) -> discarded by miss_pair
# Column 2: no gap (control); column 3: conserved (no CAAS)
GAP_ALN = {"f1": "--DA", "b1": "-KKA", "f2": "DDDA", "b2": "K-KA", "f3": "DDDA", "b3": "KKKA"}
GAP_T = ["--max_fg_gaps", "1", "--max_bg_gaps", "1", "--max_gaps", "2",
         "--max_fg_miss", "0", "--max_bg_miss", "0", "--max_miss", "0"]

# f1 (pair 1) and b2 (pair 2) are absent from the alignment: missing sets {1} vs {2}
MISS_ALN = {"f2": "DD", "b1": "KK", "f3": "DD", "b3": "KK"}
MISS_T = ["--max_fg_gaps", "0", "--max_bg_gaps", "0", "--max_gaps", "0",
          "--max_fg_miss", "1", "--max_bg_miss", "1", "--max_miss", "2"]


@pytest.mark.parametrize("name,aln,thr,dropped", [
    ("gap", GAP_ALN, GAP_T, 1),      # position 1 is discarded by miss_pair
    ("miss", MISS_ALN, MISS_T, 0),   # every position is discarded by miss_pair
])
def test_kernel_matches_scalar(tmp_path, name, aln, thr, dropped):
    with_mp = _discover(tmp_path, aln, thr, miss_pair=True)
    without = _discover(tmp_path, aln, thr, miss_pair=False)
    assert with_mp[0] == with_mp[1], f"{name}: kernel != scalar with miss_pair"
    assert without[0] == without[1], f"{name}: kernel != scalar without miss_pair"
    # non-vacuous: miss_pair changes the scalar result, and only where expected
    assert without[0] != with_mp[0]
    assert dropped not in _positions(with_mp[0]) and dropped in _positions(without[0])
