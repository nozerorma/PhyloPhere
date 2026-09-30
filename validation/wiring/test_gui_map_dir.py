"""The per-gene MAP directory is a post-processing parameter (caas_map_dir) that VEP reuses: GUI model, templates, generators."""
import copy
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.generation.validate import validate  # noqa: E402
from gui.models.serialization import from_dict, migrate  # noqa: E402
from gui.project_io import load_project  # noqa: E402

TEMPLATES = sorted((ROOT / "gui/templates").glob("*.json")) + [ROOT / "validation/tier1/templates/tier1_pepc_c4.json"]


def _project_dict(name="cancer_no_prune_multi.json"):
    return json.loads((ROOT / "gui/templates" / name).read_text())


def test_an_older_project_keeps_its_map_directory_under_the_new_name():
    d = _project_dict()
    d["modules"]["disambiguation"].pop("caas_map_dir", None)
    d["modules"]["vep"]["vep_map_dir"] = "/old/maps"
    p = from_dict(copy.deepcopy(d))
    assert p.modules.disambiguation.caas_map_dir == "/old/maps"
    assert not hasattr(p.modules.vep, "vep_map_dir")


def test_migration_does_not_overwrite_a_map_directory_already_set_in_post_processing():
    d = _project_dict()
    d["modules"]["disambiguation"]["caas_map_dir"] = "/new/maps"
    d["modules"]["vep"]["vep_map_dir"] = "/old/maps"
    assert migrate(copy.deepcopy(d))["modules"]["disambiguation"]["caas_map_dir"] == "/new/maps"


@pytest.mark.parametrize("path", TEMPLATES, ids=lambda p: p.name)
def test_every_template_carries_the_map_directory_in_post_processing_only(path):
    data = json.loads(path.read_text())
    assert "caas_map_dir" in data["modules"]["disambiguation"]
    assert "vep_map_dir" not in data["modules"]["vep"]
    assert isinstance(load_project(path).modules.disambiguation.caas_map_dir, str)


def test_the_generators_pass_the_post_processing_directory_to_nextflow_and_to_the_batch_environment():
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    proj.modules.disambiguation.caas_map_dir = "/some/maps"
    single, batch = render_single(proj), render_batch(proj)
    assert '"caas_map_dir": "${MAP_DIR:-}"' in single and "vep_map_dir" not in single
    assert 'export MAP_DIR="/some/maps"' in batch


def test_vep_requires_the_post_processing_map_directory():
    proj = load_project(ROOT / "gui/templates/cancer_no_prune_multi.json")
    proj.modules.vep.enabled = True
    proj.modules.disambiguation.caas_map_dir = ""
    assert any("per-gene MAP directory" in e for e in validate(proj))
    proj.modules.disambiguation.caas_map_dir = "/some/maps"
    assert not any("per-gene MAP directory" in e for e in validate(proj))


def _conf_default(conf, name):
    import re
    m = re.search(rf'^\s*{name}\s*=\s*"?([^"\s/]+)"?', (ROOT / "conf" / conf).read_text(), re.M)
    assert m, (conf, name)
    return m.group(1)


def test_the_single_run_script_falls_back_to_the_config_and_model_defaults():
    import re
    from gui.models.modules import CaasConfig, DisambiguationConfig
    script = (ROOT / "gui/generation/templates/run_single.sh.j2").read_text()
    strategy = re.search(r'"perm_strategy": "\$\{PERM_STRATEGY:-([^}]*)\}"', script).group(1)
    posterior = re.search(r'"ct_disambig_posterior_threshold": "\$\{POSTERIOR_THRESHOLD:-([^}]*)\}"', script).group(1)
    assert strategy == CaasConfig().perm_strategy == _conf_default("ct.config", "perm_strategy")
    assert posterior == DisambiguationConfig().ct_disambig_posterior_threshold == _conf_default(
        "ct_disambiguation.config", "ct_disambig_posterior_threshold")
