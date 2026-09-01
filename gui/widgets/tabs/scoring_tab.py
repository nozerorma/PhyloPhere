#!/usr/bin/env python3
# scoring_tab.py — Scoring module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import ScoringConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

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
        Section("Ranking cutoffs and composite weights"),
        FieldSpec(name="scoring_gene_top_pct", label="Top gene percentile"),
        FieldSpec(name="scoring_position_top_pct", label="Top position percentile"),
        FieldSpec(name="scoring_weight_caas", label="CAAS score weight"),
        FieldSpec(name="scoring_weight_rer", label="RERconverge weight"),
        FieldSpec(name="scoring_weight_fade", label="FADE weight"),
        FieldSpec(
            name="scoring_rer_direction",
            label="RER direction filter",
            kind="choice",
            choices=("both", "accelerated", "decelerated"),
        ),
        FieldSpec(name="gene_ensembl_file", label="Gene-Ensembl mapping file", kind="path_file"),
    ),
    advanced_fields=(
        Section("Downstream characterization and active modules"),
        FieldSpec(name="scoring_ami", label="Enable AMI active module report", kind="bool"),
        FieldSpec(name="scoring_string", label="Enable STRING DB integration", kind="bool"),
        FieldSpec(name="scoring_compare_fdr", label="Cross-tool comparison FDR cutoff"),
        FieldSpec(name="scoring_compare_top_n", label="Cross-tool comparison top-N genes"),
        Section("Robustness and stress testing"),
        FieldSpec(name="scoring_stress", label="Run stress-enrichment analysis", kind="bool"),
        FieldSpec(name="scoring_stress_top_n", label="Stress-enrichment top-N"),
        FieldSpec(
            name="scoring_stress_rank_metric",
            label="Stress-enrichment rank metric",
            kind="choice",
            choices=("spearman", "pearson"),
        ),
        FieldSpec(name="scoring_window_size_bp", label="Genomic window size (bp)"),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
