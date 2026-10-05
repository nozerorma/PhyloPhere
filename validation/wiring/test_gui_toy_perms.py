"""The number of permulation cycles of a toy run is a field of the Runtime tab, and the draw budget follows it."""
import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.generation.validate import validate  # noqa: E402
from gui.i18n import TRANSLATIONS  # noqa: E402
from gui.models.runtime import RuntimeConfig  # noqa: E402
from gui.models.serialization import from_dict  # noqa: E402
from gui.project_io import load_project  # noqa: E402

TEMPLATE = ROOT / "gui/templates/cancer_no_prune_multi.json"
LABEL = "Toy permulation cycles"


def _project(toy_mode=True, perms=None):
    proj = load_project(TEMPLATE)
    proj.runtime.toy_mode = toy_mode
    if perms is not None:
        proj.runtime.toy_perms = perms
    return proj


def _toy_environment(script):
    """CAAS_FULL_PERMS and MAX_TRIES as the generated script's run-config block leaves them (the block is executed)."""
    lines = script.splitlines()
    start = next(i for i, l in enumerate(lines) if re.match(r'if \[ "(true|false)" = true \]; then', l))
    end = next(i for i in range(start, len(lines)) if lines[i] == "fi")
    out = subprocess.run(["bash", "-c", "set -u\n" + "\n".join(lines[start:end + 1]) + '\necho "$CAAS_FULL_PERMS $MAX_TRIES"'],
                         capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def test_the_default_is_the_former_fixed_toy_run():
    assert RuntimeConfig().toy_perms == "100"
    assert _toy_environment(render_batch(_project())) == (100, 20000)


@pytest.mark.parametrize("perms", [1, 37, 100, 1000, 2500])
def test_the_cycles_and_the_draw_budget_move_together(perms):
    cycles, tries = _toy_environment(render_batch(_project(perms=str(perms))))
    assert cycles == perms and tries == 200 * perms


def test_a_full_run_does_not_read_the_toy_field():
    proj = _project(toy_mode=False, perms="5")
    proj.modules.caas.caas_full_perms, proj.modules.caas.max_tries = "777", "888"
    assert _toy_environment(render_batch(proj)) == (777, 888)


def test_the_script_has_no_fixed_toy_cycle_count_left():
    text = (ROOT / "gui/generation/templates/sbatch_array.sh.j2").read_text()
    assert 'CAAS_FULL_PERMS="100"' not in text and 'MAX_TRIES="20000"' not in text


@pytest.mark.parametrize("bad", ["", "abc", "-3", "1.5", " "])
def test_a_toy_run_needs_a_non_negative_integer_number_of_cycles(bad):
    errors = validate(_project(perms=bad))
    assert any("toy" in e.lower() and "cycles" in e.lower() for e in errors), errors
    assert not any("cycles" in e.lower() for e in validate(_project(toy_mode=False, perms=bad)))


@pytest.mark.parametrize("ok", ["0", "1", "1000"])
def test_a_valid_number_of_cycles_adds_no_error(ok):
    assert not any("cycles" in e.lower() for e in validate(_project(perms=ok)))


def test_no_cycle_keeps_a_positive_draw_budget():
    assert _toy_environment(render_batch(_project(perms="0"))) == (0, 200)


def test_a_project_saved_before_the_field_loads_with_the_default():
    d = json.loads(TEMPLATE.read_text())
    assert "toy_perms" not in d["runtime"]
    assert from_dict(copy.deepcopy(d)).runtime.toy_perms == "100"


def test_the_runtime_tab_has_the_field_next_to_the_toy_sample_size_and_the_label_is_translated():
    src = (ROOT / "gui/widgets/tabs/runtime_tab.py").read_text()
    # the field is wired both ways: shown from the model, written back to it, and its label retranslated
    assert "QLineEdit(self._config.toy_perms)" in src and "self.toy_perms.textChanged.connect(self._on_toy_perms_changed)" in src
    assert re.search(r"def _on_toy_perms_changed\(self, value: str\) -> None:\n\s+self\._config\.toy_perms = value\n\s+self\.changed\.emit\(\)", src)
    assert f'"{LABEL}' in src and 'self.toy_perms_label.setText(tr("Toy permulation cycles' in src
    assert src.index("self.toy_n = QLineEdit") < src.index("self.toy_perms = QLineEdit") < src.index("def _build_execution_group")
    label = next(k for k in TRANSLATIONS if k.startswith(LABEL))
    assert set(TRANSLATIONS[label]) >= {"en", "es", "ca", "fr", "it", "de"}


@pytest.mark.parametrize("bad", ["", " ", "abc", "-5", "2.5", "0"])
def test_a_toy_run_needs_a_positive_integer_sample_size(bad):
    """The pipeline reads a blank or zero toy_n as 50 alignments, so the validator refuses what would silently shrink a run."""
    proj = _project()
    proj.runtime.toy_n = bad
    errors = validate(proj)
    assert any("toy" in e.lower() and "sample" in e.lower() for e in errors), errors
    proj.runtime.toy_mode = False
    assert not any("sample" in e.lower() for e in validate(proj))


def test_every_saved_template_with_toy_mode_names_its_sample_size():
    import glob
    for f in glob.glob(str(ROOT / "gui/templates/*.json")):
        proj = load_project(Path(f))
        if proj.runtime.toy_mode:
            assert not any("sample" in e.lower() for e in validate(proj)), f
