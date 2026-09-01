#!/usr/bin/env python3
# scoring_tab.py — Scoring module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import ScoringConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Scoring",
    blurb=(
        "Combines CAAS, RERconverge, FADE, Accumulation, and VEP outputs into a "
        "composite per-gene / per-position score."
    ),
    disclaimer=(
        "Enrichment needs Scoring's gene lists/background when this is off. All of "
        "Scoring's standalone-input fallbacks (including the FADE site-level "
        "fallback) live on the Precomputed Run tab."
    ),
    essential_fields=(
        FieldSpec(name="scoring_stress", label="Run stress-enrichment analysis", kind="bool"),
        FieldSpec(name="scoring_window_size_bp", label="Genomic window size (bp)"),
        FieldSpec(name="gene_ensembl_file", label="Gene-Ensembl mapping file", kind="path_file"),
    ),
    advanced_fields=(
        FieldSpec(name="scoring_stress_top_n", label="Stress-enrichment top-N"),
        FieldSpec(
            name="scoring_stress_rank_metric",
            label="Stress-enrichment rank metric",
            kind="choice",
            choices=("spearman", "pearson"),
        ),
        FieldSpec(name="scoring_position_top_pct", label="Top position percentile"),
        FieldSpec(name="scoring_gene_top_pct", label="Top gene percentile"),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
