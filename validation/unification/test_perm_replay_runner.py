"""run_ct_perm_replay_batch.sh: two-column manifest, no discovery filter, no FOP option.

One test runs the batch runner with a stand-in `ct` that records its arguments; the other runs it with the real `ct`
on the PEPC alignment and compares the exported perm-discovery rows with a direct `ct perm-replay` call.
"""
import csv
import os
import stat
import subprocess
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
CT = ROOT / "subworkflows/CT/local/ct"
RUNNER = ROOT / "subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh"
B0_SCRIPT = ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"
GOLD = HERE / "golden/pepc_c4_complete"
ARGS = ["--patterns", "1,2,3", "--miss_pair", "--caap_mode", "--max_conserved", "2", "--max_fg_gaps", "0",
        "--max_bg_gaps", "0", "--max_gaps", "0", "--max_fg_miss", "0", "--max_bg_miss", "0", "--max_miss", "0"]


def _runner(cwd, ct, manifest, cfg, labelings, *extra):
    (cwd / "m.tsv").write_text(manifest)
    (cwd / "args.txt").write_text(" ".join(ARGS))
    return subprocess.run(["bash", str(RUNNER), "--batch-id", "b1", "--manifest", "m.tsv", "--caas-config", str(cfg),
                           "--resampled-path", str(labelings), "--workers", "1", "--ali-format", "fasta", "--ct-bin", str(ct),
                           "--export-perm-discovery", "1", "--extra-args-file", "args.txt", *extra],
                          cwd=cwd, capture_output=True, text=True)


def test_the_runner_passes_no_discovery_filter_and_no_fop_option(tmp_path):
    fake = tmp_path / "ct"
    fake.write_text('#!/usr/bin/env bash\necho "$@" >> calls.log\n')
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    (tmp_path / "alignments").mkdir()
    r = _runner(tmp_path, fake, "G1\tG1.fa\nG2\tG2.fa\n", "cfg", "lab.tab")
    assert r.returncode == 0, r.stdout + r.stderr
    calls = (tmp_path / "calls.log").read_text().splitlines()
    assert len(calls) == 2
    for c in calls:
        assert "--discovery" not in c and "--fop" not in c.split()
        assert "--export_perm_discovery" in c
    assert _runner(tmp_path, fake, "G1\tG1.fa\n", "cfg", "lab.tab", "--fop", "1").returncode != 0  # the option is gone


def _pepc_inputs(tmp_path):
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (cfg / "contrast_hypotheses_pairs.tsv").write_text((GOLD / "contrast_hypotheses_pairs.tsv").read_text())
    hyps = {}
    for r in csv.DictReader(open(cfg / "contrast_hypotheses_pairs.tsv"), delimiter="\t"):
        hyps.setdefault(r["hypothesis_id"], []).append(r)
    for h, rows in hyps.items():
        (cfg / f"traitfile_{h}.tab").write_text(
            "".join(f"{r['species1']}\t1\t{r['pair']}\n{r['species2']}\t0\t{r['pair']}\n" for r in rows))
    p = subprocess.run(["python3", str(B0_SCRIPT), "--config", str(cfg), "--fop", "--labelings-out", "b0.tab",
                        "--pairs-out", "/dev/null"], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    (tmp_path / "alignments").mkdir()
    (tmp_path / "alignments/PEPC.fa").write_text((GOLD / "PEPC.fasta").read_text())
    return cfg, tmp_path / "b0.tab"


def test_the_runner_with_the_real_ct_exports_what_a_direct_call_exports(tmp_path):
    cfg, lab = _pepc_inputs(tmp_path)
    r = _runner(tmp_path, CT, "PEPC\tPEPC.fa\n", cfg, lab)
    assert r.returncode == 0, r.stdout + r.stderr
    direct = subprocess.run([str(CT), "perm-replay", "-a", "alignments/PEPC.fa", "-t", str(cfg), "-s", str(lab), "-o", "d.out",
                             "--fmt", "fasta", *ARGS, "--export_perm_discovery", "d.disc"],
                            cwd=tmp_path, capture_output=True, text=True)
    assert direct.returncode == 0, direct.stdout + direct.stderr
    a = pd.read_csv(tmp_path / "PEPC.perm_replay.discovery.output", sep="\t")
    b = pd.read_csv(tmp_path / "d.disc", sep="\t")
    assert len(b) > 8000 and a.equals(b)
