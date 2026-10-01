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


# ── processes sized from a full-genome run (11.8k genes, 12 hypotheses x 1000 cycles) ──
#
# RESAMPLE peaked at 3.0-3.5 GB, the perm-replay step of CAAS_CORE_BATCHED had six of 804 batches taking 42-179 min
# at 8 cpus (the four longest averaged 5.5 busy cores; 8.9-14.6 GB peak depending on whether Slurm
# or Nextflow's trace measures it), and SCORING_COMPUTE peaked at 22.5-23.5 GB. A request below
# these makes an attempt fail or the slowest batch time out.

_SIZED = """
process RESAMPLE { label 'process_resample'
  script:
  \"\"\"
  echo "RESAMPLE ${task.attempt} ${task.cpus} ${task.memory.toGiga()} ${task.time.toMinutes()}" >> ${params.out}
  if [ ${task.attempt} -lt 2 ]; then exit 137; fi
  \"\"\" }
process CAAS_CORE_BATCHED { label 'process_resample'
  script:
  \"\"\"
  echo "CAAS_CORE_BATCHED ${task.attempt} ${task.cpus} ${task.memory.toGiga()} ${task.time.toMinutes()}" >> ${params.out}
  if [ ${task.attempt} -lt 2 ]; then exit 137; fi
  \"\"\" }
process SCORING_COMPUTE { label 'error_retry'
  script:
  \"\"\"
  echo "SCORING_COMPUTE ${task.attempt} ${task.cpus} ${task.memory.toGiga()} ${task.time.toMinutes()}" >> ${params.out}
  if [ ${task.attempt} -lt 2 ]; then exit 137; fi
  \"\"\" }
workflow { RESAMPLE(); CAAS_CORE_BATCHED(); SCORING_COMPUTE() }
"""


@pytest.fixture(scope="module")
def sized(tmp_path_factory):
    if shutil.which("nextflow") is None:
        pytest.skip("nextflow not on PATH")
    d = tmp_path_factory.mktemp("sized")
    (d / "main.nf").write_text(_SIZED)
    (d / "big.config").write_text("executor { cpus = 64; memory = 512.GB }\n")
    out = d / "out.txt"
    r = subprocess.run(["nextflow", "run", "main.nf", "-c", str(CONF), "-c", "big.config", "--out", str(out)],
                       cwd=d, capture_output=True, text=True, timeout=300)
    seen = {}
    for line in out.read_text().splitlines() if out.exists() else []:
        name, attempt, cpus, mem, minutes = line.split()
        seen[(name, int(attempt))] = dict(cpus=int(cpus), mem=int(mem), minutes=int(minutes))
    assert seen, r.stderr[-500:]
    return seen


def test_resample_keeps_its_large_request_and_retries_with_more(sized):
    first, second = sized[("RESAMPLE", 1)], sized[("RESAMPLE", 2)]
    assert (first["cpus"], first["mem"], first["minutes"]) == (12, 24, 12 * 60)
    assert second["mem"] == 2 * first["mem"] and second["minutes"] == 2 * first["minutes"]


def test_core_batches_have_the_cpus_and_time_the_slowest_batch_needs(sized):
    first = sized[("CAAS_CORE_BATCHED", 1)]
    assert first["cpus"] >= 4              # 4 shortens the toy stage; the slowest full-scale batches averaged 5.5 of 8 cores
    assert first["minutes"] >= 12 * 60     # 179 min at 8 cpus, and more at 4
    assert first["mem"] >= 16              # up to 14.6 GB by the trace, 8.9 GB by Slurm
    assert sized[("CAAS_CORE_BATCHED", 2)]["mem"] > first["mem"]


def test_scoring_compute_starts_above_its_measured_peak_and_still_escalates(sized):
    first, second = sized[("SCORING_COMPUTE", 1)], sized[("SCORING_COMPUTE", 2)]
    assert first["mem"] >= 32              # 22.5 GB by Slurm, 23.5 GB by the trace, plus margin
    assert second["mem"] > first["mem"]
