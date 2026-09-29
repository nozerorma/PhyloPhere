"""asr_path_score keeps every bit from the disambiguation master to filtered_discovery.tsv."""
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parents[2] / "subworkflows/CT_POSTPROC/local/src"


def test_scores_survive_prepare_and_gene_filter_unchanged(tmp_path):
    rng = random.Random(9)
    vals = [rng.random() for _ in range(3000)]
    master = tmp_path / "master.csv"
    rows = ["gene,msa_pos,caap_group,side,asr_path_score,domain_1_posterior"]
    rows += [f"g{i % 40},{i},US,top,{repr(v)},0.9" for i, v in enumerate(vals)]
    master.write_text("\n".join(rows) + "\n")

    prep = tmp_path / "prepared.tsv"
    subprocess.run([sys.executable, str(SRC / "prepare_postproc_input.py"), "--input", str(master),
                    "--output", str(prep), "--removed-output", str(tmp_path / "removed.tsv")],
                   check=True, capture_output=True, cwd=SRC)
    lens = tmp_path / "len.tsv"
    pd.DataFrame({"gene": [f"g{i}" for i in range(40)], "length": 1000}).to_csv(lens, sep="\t", index=False)
    out = tmp_path / "filtered.tsv"
    subprocess.run([sys.executable, str(SRC / "filter_caas_genes.py"), "-i", str(prep), "-l", str(lens),
                    "-m", "none", "-o", str(out)], check=True, capture_output=True, cwd=SRC)

    got = pd.read_csv(out, sep="\t", float_precision="round_trip").sort_values("Position")["asr_path_score"].tolist()
    assert got == vals
