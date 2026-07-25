#!/usr/bin/env python3
# vep_tab.py — VEP module tab (variant effect / pathogenicity annotation).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import VepConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="VEP",
    blurb=(
        "Annotates CAAS hits with PrimateAI-3D pathogenicity scores and COSMIC "
        "cancer-mutation overlap."
    ),
    disclaimer="Scoring needs a vep_caas_input fallback when this is off.",
    essential_fields=(
        FieldSpec(name="vep_primateai_db", label="PrimateAI-3D database", kind="path_file"),
        FieldSpec(name="vep_map_dir", label="Per-gene MAP directory", kind="path_dir"),
        FieldSpec(name="scoring_vep_cosmic", label="COSMIC database", kind="path_file"),
    ),
    fallback_fields=(
        FieldSpec(name="vep_caas_input", label="Precomputed CAAS input", kind="path_file"),
    ),
)


class VepTab(ModuleTabWidget):
    def __init__(self, config: VepConfig, parent=None):
        super().__init__(SPEC, config, parent)
