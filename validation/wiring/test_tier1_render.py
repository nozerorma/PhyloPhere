"""render_tier1_scripts.py: the run scripts of a Tier 1 template go to a fresh directory and never to the template's own."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
RENDER = ROOT / "validation/tier1/scripts/render_tier1_scripts.py"
TEMPLATE = ROOT / "gui/templates/tier1_pepc_c4.json"


def _render(out, n="7"):
    return subprocess.run([sys.executable, str(RENDER), "--template", str(TEMPLATE), "--outdir", str(out), "--evidence-top-n", n],
                          capture_output=True, text=True)


def test_the_scripts_land_in_the_output_directory_with_fresh_work_results_and_asr_cache(tmp_path):
    p = _render(tmp_path / "run")
    assert p.returncode == 0, p.stderr
    batch = (tmp_path / "run/run_tier1_pepc_local_complete.sh").read_text()
    single = (tmp_path / "run/run_tier1_pepc_single_complete.sh").read_text()
    assert f'export ASR_CACHE_DIR="{tmp_path}/run/asr_cache"' in batch
    assert f"{tmp_path}/run/results" in batch and f"{tmp_path}/run/work" in batch
    assert 'export CAAS_EVIDENCE_TOP_N="7"' in batch and '"caas_evidence_top_n"' in single
    assert "validation/tier1/output/pepc/" not in batch


def test_both_traits_of_the_template_are_dispatched(tmp_path):
    assert _render(tmp_path / "run").returncode == 0
    batch = (tmp_path / "run/run_tier1_pepc_local_complete.sh").read_text()
    assert "c4_phenotypic" in batch and "for task_id in $(seq 1 2)" in batch


def test_an_invalid_template_is_refused(tmp_path):
    import json
    d = json.loads(TEMPLATE.read_text())
    d["runtime"]["toy_mode"], d["runtime"]["toy_n"] = True, ""
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(d))
    p = subprocess.run([sys.executable, str(RENDER), "--template", str(bad), "--outdir", str(tmp_path / "o")], capture_output=True, text=True)
    assert p.returncode != 0 and "toy sample size" in p.stderr + p.stdout
    assert not (tmp_path / "o").exists()
