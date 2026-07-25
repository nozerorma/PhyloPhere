#!/usr/bin/env python3
# accumulation_tab.py — Accumulation module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import AccumulationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Accumulation",
    blurb=(
        "Tests whether CAAS hits accumulate on specific genes more than expected under "
        "a randomized background, using per-site Valdar-entropy conservation weighting."
    ),
    disclaimer="Enrichment's accumulation gene lists need this module's output.",
    essential_fields=(
        FieldSpec(name="accumulation_n_randomizations", label="Randomizations"),
        FieldSpec(name="accumulation_entropy_dir", label="Entropy directory", kind="path_dir"),
    ),
    fallback_fields=(
        FieldSpec(name="accumulation_caas_input", label="CAAS input file", kind="path_file"),
        FieldSpec(
            name="accumulation_background_input", label="Background input", kind="path_file"
        ),
    ),
)


class AccumulationTab(ModuleTabWidget):
    def __init__(self, config: AccumulationConfig, parent=None):
        super().__init__(SPEC, config, parent)
