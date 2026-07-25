#!/usr/bin/env python3
# test_models_roundtrip.py — dataclass <-> dict <-> JSON round-trip for ProjectConfig.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
import json

# ── Local ─────────────────────────────────────────────────────────────────────
from conftest import golden_project
from gui.models.project import ProjectConfig
from gui.models.runtime import PhenotypeRow
from gui.models.serialization import from_dict, to_dict


def _roundtrip(project: ProjectConfig) -> ProjectConfig:
    return from_dict(json.loads(json.dumps(to_dict(project))))


def test_default_project_roundtrips():
    assert _roundtrip(ProjectConfig()) == ProjectConfig()


def test_golden_project_roundtrips():
    project = golden_project()
    assert _roundtrip(project) == project


def test_project_with_disabled_module_roundtrips():
    project = golden_project()
    project.modules.rer.enabled = False
    project.modules.enrichment.posenrich_enabled = False
    assert _roundtrip(project) == project


def test_phenotype_rows_roundtrip():
    project = ProjectConfig()
    project.runtime.phenotype_rows = [
        PhenotypeRow(trait_class=1, trait="a", secondary="b"),
        PhenotypeRow(trait_class=2, trait="c", discrete_method="decile"),
    ]
    restored = _roundtrip(project)
    assert restored.runtime.phenotype_rows == project.runtime.phenotype_rows


def test_json_output_is_human_diffable_dict():
    data = to_dict(golden_project())
    # Top-level keys follow field-declaration order (schema_version, flavor, ...).
    assert list(data.keys())[:2] == ["schema_version", "flavor"]
    assert data["flavor"] == "primates"
    assert isinstance(data["runtime"]["phenotype_rows"], list)


def test_unsupported_schema_version_raises():
    import pytest

    from gui.models.serialization import from_dict

    data = to_dict(ProjectConfig())
    data["schema_version"] = 999
    with pytest.raises(ValueError):
        from_dict(data)
