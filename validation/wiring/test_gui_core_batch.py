"""ct_core_batch_size is the one batch-size parameter of the permulation null: GUI model, tab, templates, generators, migration."""
import copy
import json
import logging
import os
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.models.modules import CaasConfig  # noqa: E402
from gui.models.serialization import from_dict  # noqa: E402
from gui.project_io import load_project  # noqa: E402
from test_gui_map_dir import _conf_default, _tab_field_names  # noqa: E402

OLD = ("ct_perm_replay_batch_size", "ct_disambig_perms_batch_size")
TEMPLATES = sorted((ROOT / "gui/templates").glob("*.json")) + [ROOT / "validation/tier1/templates/tier1_pepc_c4.json"]
REMOVED_SELECTORS = {"process_perm_replay", "process_perm_replay_batched", "CAAS_PERMS_DISAMBIGUATE", "CAAS_PERMS_DISAMBIGUATE_BATCHED",
                     "CAAS_PERMS_AGGREGATE", "CAAS_PERMS_MERGE_DETAIL", "CAAS_PERMS_REBUILD"}


def test_the_model_has_one_batch_size_for_the_null_with_the_config_default():
    fields = set(CaasConfig.__dataclass_fields__)
    assert "ct_core_batch_size" in fields and not fields & set(OLD)
    assert CaasConfig().ct_core_batch_size == _conf_default("ct.config", "ct_core_batch_size")


def test_the_caas_tab_shows_the_single_batch_size():
    names = _tab_field_names("caas_tab.py", "essential_fields") + _tab_field_names("caas_tab.py", "advanced_fields")
    assert "ct_core_batch_size" in names and not set(names) & set(OLD)


@pytest.mark.parametrize("path", TEMPLATES, ids=lambda p: p.name)
def test_every_template_has_the_single_batch_size_and_no_selector_of_a_removed_process(path):
    data = json.loads(path.read_text())
    caas = data["modules"]["caas"]
    assert caas["ct_core_batch_size"] == "20" and not set(caas) & set(OLD)
    selectors = {row["selector"] for row in data["resources"].get("process_overrides", [])}
    assert not selectors & REMOVED_SELECTORS
    assert load_project(path).modules.caas.ct_core_batch_size == "20"


def test_the_generators_write_and_export_the_single_batch_size():
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    proj.modules.caas.ct_core_batch_size = "7"
    single, batch = render_single(proj), render_batch(proj)
    default = re.search(r'"ct_core_batch_size": "\$\{CT_CORE_BATCH_SIZE:-([^}]*)\}"', single).group(1)
    assert default == CaasConfig().ct_core_batch_size
    assert 'export CT_CORE_BATCH_SIZE="7"' in batch
    for text in (single, batch):
        assert not any(k in text or k.upper() in text for k in OLD)


def test_a_project_saved_with_the_old_batch_sizes_loads_with_the_default_and_a_warning(caplog):
    d = json.loads((ROOT / "gui/templates/cancer_no_prune_multi.json").read_text())
    d["modules"]["caas"].pop("ct_core_batch_size", None)
    d["modules"]["caas"].update({"ct_perm_replay_batch_size": "50", "ct_disambig_perms_batch_size": "5"})
    with caplog.at_level(logging.WARNING):
        p = from_dict(copy.deepcopy(d))
    assert p.modules.caas.ct_core_batch_size == CaasConfig().ct_core_batch_size
    assert not any(hasattr(p.modules.caas, k) for k in OLD)
    warned = " ".join(r.getMessage() for r in caplog.records)
    assert all(k in warned for k in OLD) and "ct_core_batch_size" in warned


def test_a_project_with_the_new_batch_size_loads_without_a_warning(caplog):
    d = json.loads((ROOT / "gui/templates/cancer_no_prune_multi.json").read_text())
    d["modules"]["caas"]["ct_core_batch_size"] = "33"
    with caplog.at_level(logging.WARNING):
        p = from_dict(copy.deepcopy(d))
    assert p.modules.caas.ct_core_batch_size == "33" and not caplog.records


# ── b_0 is always replayed: the switch that enabled it is gone ───────────────

def test_the_b0_diagnostic_switch_is_gone_from_model_tab_templates_and_generators():
    assert "caas_b0_diagnostic" not in CaasConfig.__dataclass_fields__
    names = _tab_field_names("caas_tab.py", "essential_fields") + _tab_field_names("caas_tab.py", "advanced_fields")
    assert "caas_b0_diagnostic" not in names
    for path in TEMPLATES:
        assert "caas_b0_diagnostic" not in json.loads(path.read_text())["modules"]["caas"], path.name
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    assert "b0_diagnostic" not in render_single(proj) + render_batch(proj)
    assert not re.search(r"^\s*caas_b0_diagnostic\s*=", (ROOT / "conf/ct.config").read_text(), re.M)


def test_a_project_saved_with_the_b0_diagnostic_switch_loads_with_a_warning(caplog):
    d = json.loads((ROOT / "gui/templates/cancer_no_prune_multi.json").read_text())
    d["modules"]["caas"]["caas_b0_diagnostic"] = True
    with caplog.at_level(logging.WARNING):
        p = from_dict(copy.deepcopy(d))
    assert not hasattr(p.modules.caas, "caas_b0_diagnostic")
    assert "caas_b0_diagnostic" in " ".join(r.getMessage() for r in caplog.records)


# ── parameters nothing reads: the per-gene discovery batch, the observed disambiguation batch, the ASR mode ──

RETIRED_PARAMS = {"caas": ("ct_discovery_batch_size",), "disambiguation": ("ct_disambig_asr_mode", "ct_disambig_batch_size")}
_ALL_RETIRED = [name for names in RETIRED_PARAMS.values() for name in names]


def test_the_retired_parameters_are_gone_from_the_model_the_tabs_and_the_config():
    from gui.models.modules import DisambiguationConfig
    for model, module in ((CaasConfig, "caas"), (DisambiguationConfig, "disambiguation")):
        assert not set(model.__dataclass_fields__) & set(RETIRED_PARAMS[module])
    names = (_tab_field_names("caas_tab.py", "essential_fields") + _tab_field_names("caas_tab.py", "advanced_fields")
             + _tab_field_names("disambiguation_tab.py", "essential_fields") + _tab_field_names("disambiguation_tab.py", "advanced_fields"))
    assert not set(names) & set(_ALL_RETIRED)
    conf = " ".join((ROOT / "conf" / f).read_text() for f in ("ct.config", "ct_disambiguation.config"))
    assert not any(re.search(rf"^\s*{name}\s*=", conf, re.M) for name in _ALL_RETIRED)


@pytest.mark.parametrize("path", TEMPLATES, ids=lambda p: p.name)
def test_no_template_carries_a_retired_parameter(path):
    modules = json.loads(path.read_text())["modules"]
    for module, names in RETIRED_PARAMS.items():
        assert not set(modules[module]) & set(names), (path.name, module)


def test_the_generators_write_none_of_the_retired_parameters():
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    text = render_single(proj) + render_batch(proj)
    for name in _ALL_RETIRED:
        assert not re.search(rf"\b{name}\b", text) and not re.search(rf"\b{name.upper()}\b", text), name
    assert "ct_disambig_asr_cache_dir" in text and "ct_disambig_asr_model" in text  # the ones that stay


def test_a_project_saved_with_the_retired_parameters_loads_with_one_warning_per_module(caplog):
    d = json.loads((ROOT / "gui/templates/cancer_no_prune_multi.json").read_text())
    d["modules"]["caas"]["ct_discovery_batch_size"] = "100"
    d["modules"]["disambiguation"].update({"ct_disambig_asr_mode": "compute", "ct_disambig_batch_size": "20"})
    with caplog.at_level(logging.WARNING):
        p = from_dict(copy.deepcopy(d))
    assert not any(hasattr(p.modules.caas, n) for n in RETIRED_PARAMS["caas"])
    assert not any(hasattr(p.modules.disambiguation, n) for n in RETIRED_PARAMS["disambiguation"])
    warned = " ".join(r.getMessage() for r in caplog.records)
    assert all(n in warned for n in _ALL_RETIRED)


def test_the_asr_cache_directory_is_always_required_when_disambiguation_runs():
    from gui.generation.validate import validate
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    proj.modules.disambiguation.ct_disambig_asr_cache_dir = ""
    assert any("ASR cache directory is required" in e for e in validate(proj))
    proj.modules.disambiguation.ct_disambig_asr_cache_dir = "/some/cache"
    assert not any("ASR cache directory is required" in e for e in validate(proj))
