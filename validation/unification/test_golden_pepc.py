"""b_0 through the perm-replay kernel vs the frozen CAAStools discovery of the Tier 1 PEPC run.

The golden files (golden/pepc_c4_complete) are the observed `discovery.tab` and `background.output`
of the genotypic PEPC run (100 hypotheses, 4 pairs, patterns 1,2,3, miss_pair, caap_mode,
no gaps or missing allowed, up to 2 conserved pairs). The kernel replays the same 100 hypotheses
as the labelings `b_0~H<m>` and must reproduce those rows and the tested-position list exactly.
"""
import csv
import gzip
import subprocess
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CT = ROOT / "subworkflows/CT/local/ct"
B0_SCRIPT = ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"
GOLD = HERE / "golden/pepc_c4_complete"

KEY = ["gene", "caap_group", "hyp", "position", "caas", "amino_encoded"]


def _run(cmd, cwd):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr


def test_b0_kernel_matches_frozen_discovery(tmp_path):
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (cfg / "contrast_hypotheses_pairs.tsv").write_text((GOLD / "contrast_hypotheses_pairs.tsv").read_text())
    # one traitfile per hypothesis (species, 1 = fg / 0 = bg, pair id), as the pipeline writes them
    hyps = {}
    for r in csv.DictReader(open(cfg / "contrast_hypotheses_pairs.tsv"), delimiter="\t"):
        hyps.setdefault(r["hypothesis_id"], []).append(r)
    for h, rows in hyps.items():
        (cfg / f"traitfile_{h}.tab").write_text(
            "".join(f"{r['species1']}\t1\t{r['pair']}\n{r['species2']}\t0\t{r['pair']}\n" for r in rows))
    assert len(hyps) == 100

    _run(["python3", B0_SCRIPT, "--config", cfg, "--fop", "--labelings-out", "b0.tab", "--pairs-out", "/dev/null"], tmp_path)
    _run([CT, "perm-replay", "-a", GOLD / "PEPC.fasta", "-t", cfg, "-s", "b0.tab", "--fmt", "fasta",
          "--patterns", "1,2,3", "--miss_pair", "--caap_mode", "--max_conserved", "2",
          "--max_fg_gaps", "0", "--max_bg_gaps", "0", "--max_gaps", "0",
          "--max_fg_miss", "0", "--max_bg_miss", "0", "--max_miss", "0",
          "--export_perm_discovery", "k.disc", "--export_b0_background", "k.bg"], tmp_path)

    gold = pd.read_csv(GOLD / "discovery.tab.gz", sep="\t")
    gold["hyp"] = gold["trait"].str.extract(r"(H\d+)")[0]
    kern = pd.read_csv(tmp_path / "k.disc", sep="\t")
    kern["hyp"] = kern["cycle"].str.extract(r"(H\d+)")[0]
    rows = lambda d: set(map(tuple, d[KEY].astype(str).itertuples(index=False, name=None)))
    assert len(gold) > 8000  # guard against an empty golden
    assert rows(kern) == rows(gold)
    # kernel positions are the same 0-based ids as the scalar path
    assert len(kern) == len(gold)

    gene, tested = gzip.open(GOLD / "background.output.gz", "rt").read().rstrip("\n").split("\t")
    k_gene, k_tested = (tmp_path / "k.bg").read_text().rstrip("\n").split("\t")
    assert (k_gene, set(k_tested.split(","))) == (gene, set(tested.split(",")))
