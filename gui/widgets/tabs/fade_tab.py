#!/usr/bin/env python3
# fade_tab.py — FADE module tab (directional selection, HyPhy).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""Off by default (RUN_FADE=false in the reference scripts)."""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import FadeConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="FADE",
    blurb=(
        "Directional selection analysis (HyPhy FADE) across the phylogeny, run on "
        "either every gene or just the CAAS hit set (fade_mode)."
    ),
    disclaimer=(
        "Scoring needs a scoring_fade_summary_top/bottom fallback per phenotype row "
        "(Runtime tab) when this is off, plus optionally scoring_fade_site_top/bottom "
        "(Scoring tab) for site-level fallback."
    ),
    essential_fields=(
        FieldSpec(name="fade_mode", label="Mode", kind="choice", choices=("all", "gene_set")),
    ),
)


class FadeTab(ModuleTabWidget):
    def __init__(self, config: FadeConfig, parent=None):
        super().__init__(SPEC, config, parent)
