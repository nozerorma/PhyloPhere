#!/usr/bin/env python3
# accumulation_tab.py — Accumulation module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import AccumulationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Accumulation",
    blurb=(
        "Tests whether CAAS hits accumulate on specific genes more than expected under "
        "a randomized background, using per-site Valdar-entropy conservation weighting."
    ),
    disclaimer=(
        "Enrichment's accumulation gene lists need this module's output. Check 'Use "
        "precomputed Accumulation output' on the Precomputed Run tab instead."
    ),
    essential_fields=(
        Section("Randomization and burden parameters"),
        FieldSpec(
            name="accumulation_randomization_type",
            label="Randomization type",
            kind="choice",
            choices=("cons_decile", "naive"),
        ),
        FieldSpec(name="accumulation_n_randomizations", label="Randomizations"),
        FieldSpec(name="accumulation_fdr", label="FDR threshold"),
    ),
    advanced_fields=(
        Section("Entropy and standalone inputs"),
        FieldSpec(name="accumulation_entropy_dir", label="Entropy directory", kind="path_dir"),
    ),
)


class AccumulationTab(ModuleTabWidget):
    def __init__(self, config: AccumulationConfig, parent=None):
        super().__init__(SPEC, config, parent)
