"""review_gate_2_2.sh: prints the three comparison commands, refuses the cluster login node, and is valid bash."""
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "review_gate_2_2.sh"


def _run(*args, **env):
    e = {k: v for k, v in os.environ.items() if k not in ("SLURM_JOB_ID", "ALLOW_LOGIN_NODE", "REVIEW_HOSTNAME")}
    e.update(env)
    return subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True, env=e)


def test_the_script_is_valid_bash():
    assert subprocess.run(["bash", "-n", str(SCRIPT)]).returncode == 0


def test_dry_run_prints_the_three_commands_with_the_run_directories():
    p = _run("--dry-run", "/runs/base", "/runs/new")
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert len(lines) == 3
    assert "compare_b0.py --run /runs/new --out /runs/new/gate_2_2/compare_b0.json" in lines[0]
    assert "compare_null.py --a /runs/base --b /runs/new --extra scoring/position_scores.tsv --extra scoring/gene_scores.tsv" in lines[1]
    assert "compare_contract.py --a /runs/base --b /runs/new --report /runs/new/gate_2_2/compare_contract.json" in lines[2]


def test_a_report_directory_can_be_given_and_python_chosen():
    p = _run("--dry-run", "/b", "/n", "/reports", PYTHON="/opt/py")
    assert p.returncode == 0 and all(l.startswith("/opt/py ") for l in p.stdout.splitlines()) and "/reports/compare_b0.json" in p.stdout


def test_the_login_node_is_refused_unless_a_job_or_the_override_is_present():
    refused = _run("--dry-run", "/b", "/n", REVIEW_HOSTNAME="correfoc-01.s.upf.edu")
    assert refused.returncode == 2 and "login node" in refused.stderr
    assert _run("--dry-run", "/b", "/n", REVIEW_HOSTNAME="correfoc-01", SLURM_JOB_ID="123").returncode == 0
    assert _run("--dry-run", "/b", "/n", REVIEW_HOSTNAME="correfoc-01", ALLOW_LOGIN_NODE="1").returncode == 0
    assert _run("--dry-run", "/b", "/n", REVIEW_HOSTNAME="cr-07-01").returncode == 0  # a compute node


def test_a_wrong_number_of_arguments_prints_the_usage():
    p = _run("/only/one")
    assert p.returncode == 2 and "review_gate_2_2.sh" in p.stderr
