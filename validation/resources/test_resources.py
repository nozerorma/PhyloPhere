"""Resource configuration: one source of per-process values, effective ceilings, no hardcoded partition.

Run: python3 -m pytest -q validation/resources/
The Nextflow tests are skipped when `nextflow` is not on PATH. They only use the `local`
profile, so running them on a cluster never submits a Slurm job.
"""
import glob
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gui.models.project import ProjectConfig  # noqa: E402
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.resource_defaults import DEFAULTS_FILE, load_defaults  # noqa: E402

CONF = ROOT / "conf" / "resources.config"
SELECTOR = re.compile(r"with(?:Name|Label):?\s*'?([A-Za-z0-9_]+)'?\s*\{")


# ── one source of per-process values ─────────────────────────────────────────

def test_the_conf_is_the_only_file_with_per_process_values():
    assert glob.glob(str(ROOT / "conf" / "resources.config.*")) == []


def test_defaults_parser_reads_the_conf_and_skips_alternation_selectors():
    rows = load_defaults()
    selectors = {r.selector for r in rows}
    assert DEFAULTS_FILE == CONF and len(rows) > 50
    assert {"process_reporting", "SCORING_REPORT"} <= selectors
    assert all(r.cpus == "" or r.cpus.isdigit() for r in rows)
    assert not any("|" in r.selector for r in rows)
    # a memory closure keeps its retry scaling, or a retry would not get more memory
    assert any("task.attempt" in r.memory for r in rows)


def test_saved_templates_carry_no_override_that_targets_a_missing_process():
    """An override for a selector nothing defines applies to nothing. Templates should hold none, or only real ones."""
    known = set(SELECTOR.findall(CONF.read_text()))
    known |= set(re.findall(r"process\s+([A-Z][A-Z0-9_]+)\s*\{", "\n".join(p.read_text() for p in ROOT.rglob("*.nf") if "validation" not in p.parts)))
    for t in glob.glob(str(ROOT / "gui" / "templates" / "*.json")):
        rows = json.load(open(t)).get("resources", {}).get("process_overrides", [])
        missing = [r["selector"] for r in rows if r["selector"] not in known]
        assert not missing, (t, missing)


def test_saved_templates_do_not_freeze_a_copy_of_the_defaults():
    for t in glob.glob(str(ROOT / "gui" / "templates" / "*.json")):
        rows = json.load(open(t)).get("resources", {}).get("process_overrides", [])
        assert len(rows) < 10, f"{t} carries {len(rows)} override rows; they pin values the conf should own"


# ── ceilings ─────────────────────────────────────────────────────────────────

def test_profiles_enforce_the_ceilings_and_the_dead_helper_is_gone():
    cfg = (ROOT / "nextflow.config").read_text()
    assert "def check_max" not in cfg
    assert cfg.count("resourceLimits") == 2  # slurm and local profiles


needs_nextflow = pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")

_MAIN = """
process BIG { cpus 16; memory 32.GB; time 100.h
  output: stdout
  script: "echo BIG cpus=${task.cpus} mem=${task.memory} time=${task.time}" }
process SMALL { cpus 2; memory 2.GB; time 1.h
  output: stdout
  script: "echo SMALL cpus=${task.cpus} mem=${task.memory} time=${task.time}" }
workflow { BIG().view(); SMALL().view() }
"""


def _run_local(tmp_path, *params):
    (tmp_path / "main.nf").write_text(_MAIN)
    r = subprocess.run(["nextflow", "run", "main.nf", "-c", str(ROOT / "nextflow.config"), "-profile", "local", *params],
                       cwd=tmp_path, capture_output=True, text=True, timeout=240)
    return {m.group(1): m.group(2) for m in re.finditer(r"^(BIG|SMALL) (.*)$", r.stdout, re.M)}, r


@needs_nextflow
def test_ceilings_clamp_requests_and_leave_small_ones_alone(tmp_path):
    got, r = _run_local(tmp_path, "--max_cpus", "4", "--max_memory", "8.GB", "--max_time", "90.m")
    assert got.get("BIG") == "cpus=4 mem=8 GB time=1h 30m", r.stderr[-400:]
    assert got.get("SMALL") == "cpus=2 mem=2 GB time=1h"


@needs_nextflow
def test_ceilings_do_not_touch_requests_below_them(tmp_path):
    got, r = _run_local(tmp_path)
    assert got.get("BIG") == "cpus=16 mem=32 GB time=4d 4h", r.stderr[-400:]


# ── SLURM partition ──────────────────────────────────────────────────────────

def _project(partition, batched):
    p = ProjectConfig()
    p.runtime.runtime_type = "slurm"
    p.runtime.sbatch_partition = partition
    p.runtime.batched = batched
    return p


def test_default_partition_is_empty_and_never_haswell():
    assert ProjectConfig().runtime.sbatch_partition == ""
    assert "haswell" not in (ROOT / "gui" / "generation" / "templates" / "run_single.sh.j2").read_text()


@pytest.mark.parametrize("partition", ["", "high-cpu"])
def test_single_runner_takes_the_queue_from_the_gui_partition(partition):
    text = render_single(_project(partition, batched=False))
    assert f'"slurm_queue": "${{SLURM_QUEUE:-{partition}}}"' in text
    assert "haswell" not in text


def test_batch_wrapper_omits_the_partition_line_when_empty_and_keeps_the_next_directive():
    empty = render_batch(_project("", batched=True))
    assert "#SBATCH --partition" not in empty and "#SBATCH -t 144:00:00" in empty
    assert "#SBATCH --partition=high-cpu" in render_batch(_project("high-cpu", batched=True))
