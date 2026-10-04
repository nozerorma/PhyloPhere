"""`ct perm-replay` takes one labelings file and writes exports; the scalar discovery, the counts file, the
discovery-position filter, the directory mode and the progress log are gone.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
CT = ROOT / "subworkflows/CT/local/ct"
sys.path.insert(0, str(HERE))
from test_perm_replay_runner import ARGS, _pepc_inputs  # noqa: E402


def _ct(cwd, *args):
    return subprocess.run([str(CT), "perm-replay", *args], cwd=cwd, capture_output=True, text=True)


def _base(cfg, lab):
    return ["-a", "alignments/PEPC.fa", "-t", str(cfg), "-s", str(lab), "--fmt", "fasta", *ARGS]


def test_a_run_without_an_export_is_refused(tmp_path):
    cfg, lab = _pepc_inputs(tmp_path)
    p = _ct(tmp_path, *_base(cfg, lab))
    assert p.returncode != 0 and "No export requested" in p.stdout + p.stderr


def test_a_directory_of_labelings_is_refused(tmp_path):
    cfg, lab = _pepc_inputs(tmp_path)
    (tmp_path / "labelings").mkdir()
    p = _ct(tmp_path, *_base(cfg, tmp_path / "labelings"), "--export_b0_background", "bg")
    assert p.returncode != 0 and "one labelings file" in p.stdout + p.stderr


@pytest.mark.parametrize("option", [["-o", "x.out"], ["--output", "x.out"], ["--fop"], ["--discovery", "d.tab"], ["--progress_log", "p.log"]])
def test_the_removed_options_are_rejected(tmp_path, option):
    cfg, lab = _pepc_inputs(tmp_path)
    p = _ct(tmp_path, *_base(cfg, lab), "--export_b0_background", "bg", *option)
    assert p.returncode != 0 and "no such option" in p.stderr


def test_the_scalar_discovery_is_not_a_tool_of_ct(tmp_path):
    p = subprocess.run([str(CT), "discovery", "-a", "x"], cwd=tmp_path, capture_output=True, text=True)
    assert "no tool named discovery" in p.stdout + p.stderr


def test_the_background_alone_is_enough_and_equals_the_one_written_with_the_other_exports(tmp_path):
    cfg, lab = _pepc_inputs(tmp_path)
    only = _ct(tmp_path, *_base(cfg, lab), "--export_b0_background", "only.bg")
    both = _ct(tmp_path, *_base(cfg, lab), "--export_b0_background", "both.bg", "--export_b0_discovery", "both.discovery")
    assert only.returncode == 0 and both.returncode == 0, only.stderr + both.stderr
    assert (tmp_path / "only.bg").read_text() == (tmp_path / "both.bg").read_text() and (tmp_path / "only.bg").read_text().startswith("PEPC\t")


def test_the_module_keeps_the_helpers_the_kernel_uses_and_none_of_the_removed_ones():
    sys.path.insert(0, str(ROOT / "subworkflows/CT/local"))
    import importlib
    pr = importlib.import_module("modules.perm_replay")
    assert pr._fop_base_cycle("b_12~H3") == "b_12" and pr._FOP_H_SUFFIX.search("b_1~H2")
    for gone in ("parse_discovery_positions", "pval", "collapse_fop_hits_by_base", "log_progress", "format_time", "calculate_eta"):
        assert not hasattr(pr, gone), gone
    assert "fop_mode" not in pr.run_perm_replay_on_alignment.__code__.co_varnames
    with pytest.raises(ImportError):
        importlib.import_module("modules.disco")
    with pytest.raises(ValueError, match="at least one export"):
        pr.run_perm_replay_on_alignment("cfg", object(), object(), "NO", "NO", "NO", "NO", "NO", "NO", "1,2,3")
