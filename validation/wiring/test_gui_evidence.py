"""The number of positions to explain after SCORING (caas_evidence_top_n) reaches every surface of the pipeline:
the model, the Scoring tab, both generated scripts, the validator, the help and the README."""
import copy
import json
import os
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.generation.validate import validate  # noqa: E402
from gui.i18n import TRANSLATIONS  # noqa: E402
from gui.models.modules import ScoringConfig  # noqa: E402
from gui.models.serialization import from_dict  # noqa: E402
from gui.project_io import load_project  # noqa: E402

TEMPLATE = ROOT / "gui/templates/cancer_no_prune_multi.json"
LABEL = "Evidence of the N best positions"


def _project(n=None, **changes):
    proj = load_project(TEMPLATE)
    if n is not None:
        proj.modules.scoring.caas_evidence_top_n = n
    for dotted, value in changes.items():
        obj, attr = dotted.rsplit("__", 1)
        target = proj
        for part in obj.split("__"):
            target = getattr(target, part)
        setattr(target, attr, value)
    return proj


def _evidence_errors(proj):
    return [e for e in validate(proj) if "evidence" in e.lower()]


def test_the_default_is_off_and_the_templates_validate_with_it():
    assert ScoringConfig().caas_evidence_top_n == "0"
    assert not _evidence_errors(_project())


def test_both_scripts_carry_the_value_to_the_parameter():
    assert 'export CAAS_EVIDENCE_TOP_N="10"' in render_batch(_project("10"))
    assert '"caas_evidence_top_n": "${CAAS_EVIDENCE_TOP_N:-0}"' in render_single(_project("10"))


def test_the_single_script_reads_the_variable_the_batch_script_exports():
    single = render_single(_project())
    batch = render_batch(_project("7"))
    var = re.search(r'"caas_evidence_top_n": "\$\{(\w+):-', single).group(1)
    assert f'export {var}="7"' in batch


@pytest.mark.parametrize("bad", ["", "abc", "-1", "2.5", " "])
def test_the_number_of_positions_must_be_a_non_negative_integer(bad):
    assert any("integer" in e for e in _evidence_errors(_project(bad)))


@pytest.mark.parametrize("ok", ["0", "1", "10"])
def test_a_valid_number_adds_no_error_on_a_live_run(ok):
    assert not _evidence_errors(_project(ok))


def test_evidence_needs_scoring():
    errors = _evidence_errors(_project("5", modules__scoring__enabled=False))
    assert any("Scoring" in e for e in errors), errors
    assert not _evidence_errors(_project("0", modules__scoring__enabled=False))


def test_evidence_needs_the_observed_discovery_of_the_run_or_a_reused_one():
    live_off = _project("5", modules__caas__enabled=False)
    assert any("discovery" in e for e in _evidence_errors(live_off))
    reused = _project("5", modules__caas__enabled=False, precomputed__use_discovery=True, precomputed__base_path="/x")
    assert not any("discovery" in e for e in _evidence_errors(reused))
    no_scoring_input = _project("5", modules__disambiguation__enabled=False)
    assert any("discovery" in e for e in _evidence_errors(no_scoring_input))


def test_a_project_saved_before_the_field_loads_with_the_default():
    d = json.loads(TEMPLATE.read_text())
    assert "caas_evidence_top_n" not in d["modules"]["scoring"]
    assert from_dict(copy.deepcopy(d)).modules.scoring.caas_evidence_top_n == "0"


def test_the_scoring_tab_shows_the_field_and_its_label_is_translated():
    src = (ROOT / "gui/widgets/tabs/scoring_tab.py").read_text()
    assert re.search(r'FieldSpec\(name="caas_evidence_top_n", label="Evidence of the N best positions', src)
    label = next(k for k in TRANSLATIONS if k.startswith(LABEL))
    assert set(TRANSLATIONS[label]) >= {"en", "es", "ca", "fr", "it", "de"}


def test_the_parameter_is_documented_where_the_others_are():
    assert re.search(r"^\s*caas_evidence_top_n\s*=\s*0\b", (ROOT / "conf/scoring.config").read_text(), re.M)
    assert re.search(r"^--caas_evidence_top_n\s", (ROOT / "workflows/help.nf").read_text(), re.M)
    assert "`caas_evidence_top_n`" in (ROOT / "README.md").read_text()
