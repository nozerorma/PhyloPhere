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
        FieldSpec(name="scoring_gene_top_pct", label="Top gene percentile", importance="default"),
        FieldSpec(name="scoring_position_top_pct", label="Top position percentile", importance="default"),
        FieldSpec(name="scoring_weight_caas", label="CAAS score weight", importance="default"),
        FieldSpec(name="scoring_weight_rer", label="RERconverge weight", importance="default"),
        FieldSpec(name="scoring_weight_fade", label="FADE weight", importance="default"),
        FieldSpec(
            name="scoring_rer_direction",
            label="RER direction filter",
            kind="choice",
            choices=("both", "accelerated", "decelerated"),
            importance="default",
        ),
        FieldSpec(
            name="gene_ensembl_file",
            label="Gene-Ensembl mapping file",
            kind="path_file",
            importance="optional",
        ),
        FieldSpec(
            name="auto_generate_ensembl",
            label="Auto-generate Ensembl mapping via BioMart if unset",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(
            name="ensembl_dataset",
            label="Ensembl BioMart dataset (optional)",
            importance="optional",
        ),
    ),
    advanced_fields=(
        Section("Downstream characterization and active modules"),
        FieldSpec(name="scoring_ami", label="Enable AMI active module report", kind="bool", importance="optional"),
        FieldSpec(name="scoring_string", label="Enable STRING DB integration", kind="bool", importance="optional"),
        FieldSpec(name="scoring_compare_fdr", label="Cross-tool comparison FDR cutoff", importance="default"),
        FieldSpec(name="scoring_compare_top_n", label="Cross-tool comparison top-N genes", importance="optional"),
        Section("Robustness and stress testing"),
        FieldSpec(name="scoring_stress", label="Run stress-enrichment analysis", kind="bool", importance="optional"),
        FieldSpec(name="scoring_stress_top_n", label="Stress-enrichment top-N", importance="optional"),
        FieldSpec(
            name="scoring_stress_rank_metric",
            label="Stress-enrichment rank metric",
            kind="choice",
            choices=("spearman", "pearson"),
            importance="default",
        ),
        FieldSpec(name="scoring_window_size_bp", label="Genomic window size (bp)", importance="optional"),
        FieldSpec(
            name="scoring_pos_perm_p_thr",
            label="Position-level permulation p (BH) threshold",
            importance="default",
        ),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
