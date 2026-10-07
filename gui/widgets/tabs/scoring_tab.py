#!/usr/bin/env python3
# scoring_tab.py — Scoring module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
ScoringTab: field specification of the Scoring module (conf/scoring.config).

Declares the ModuleTabSpec and binds it to ScoringConfig through ModuleTabWidget.

Imported by: gui/widgets/main_window.py
"""

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
        Section("Ranking cutoffs"),
        FieldSpec(name="scoring_gene_top_pct", label="Top gene percentile", importance="default"),
        FieldSpec(name="scoring_position_top_pct", label="Top position percentile", importance="default"),
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
        FieldSpec(
            name="scoring_hypotheses_pairs",
            label="Contrast hypotheses pairs file for scoring (optional)",
            kind="path_file",
            importance="optional",
        ),
        Section("Evidence of the best positions"),
        FieldSpec(name="caas_evidence_top_n", label="Evidence of the N best positions (0 = off)", importance="optional"),
    ),
    advanced_fields=(
        Section("Genomic windows and position-level significance"),
        FieldSpec(
            name="caas_score_aggregation",
            label="Scheme aggregation of the position score (cumulative: sum over the 5 schemes, missing = 0)",
            kind="choice",
            choices=("cumulative", "mean"),
            importance="default",
        ),
        FieldSpec(name="scoring_window_size_bp", label="Genomic window size (bp)", importance="optional"),
        FieldSpec(
            name="scoring_p_emp_thr",
            label="Position-level permulation p.adj_bh / p.adj_sam threshold",
            importance="default",
        ),
        FieldSpec(name="scoring_gene_perm_pooled", label="Pooled-null gene permulation p (n-stratified, opt-in)", kind="bool", importance="optional"),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
