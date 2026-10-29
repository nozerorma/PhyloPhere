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
        "Combines CAAS, RERconverge, FADE, and VEP outputs into a "
        "composite per-gene / per-position score."
    ),
    disclaimer=(
        "Enrichment needs Scoring's gene lists/background when this is off. All of "
        "Scoring's standalone-input fallbacks (including the FADE site-level "
        "fallback) live on the Precomputed Run tab."
    ),
    essential_fields=(
        Section("Empirical evidence (CAAS and RER)"),
        FieldSpec(name="evidence_p_threshold", label="Nominal empirical p-value for activity", importance="default", help="Nominal p-value cutoff that selects active CAAS genes and positions and RER genes; FCS keeps the full ranking."),
        FieldSpec(name="evidence_fdr_threshold", label="BH threshold for highlighting", importance="default", help="BH-adjusted cutoff that highlights the subset in report tables and AMI networks; does not add a module run."),
        FieldSpec(
            name="gene_ensembl_file", help="Gene annotation TSV (gene, coordinates, length, human_protein_id); blank is filled from BioMart only when auto-generation is on.",
            label="Gene-Ensembl mapping file",
            kind="path_file",
            importance="optional",
        ),
        FieldSpec(
            name="auto_generate_ensembl", help="Build gene_ensembl_file from BioMart when it is blank; needs network access.",
            label="Auto-generate Ensembl mapping via BioMart if unset",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(
            name="ensembl_dataset", help="BioMart dataset for auto-generation; blank derives it from the reference species (Homo_sapiens gives hsapiens_gene_ensembl).",
            label="Ensembl BioMart dataset (optional)",
            importance="optional",
        ),
        FieldSpec(
            name="scoring_hypotheses_pairs", help="contrast_hypotheses_pairs.tsv of per-pair PSS weights; blank uses this run's contrast selection (standalone runs auto-discover it).",
            label="Contrast hypotheses pairs file for scoring (optional)",
            kind="path_file",
            importance="optional",
        ),
        Section("Evidence of the best positions"),
        FieldSpec(name="caas_evidence_top_n", help="Number of best positions in the evidence table written to scoring/evidence/; 0 = off.", label="Evidence of the N best positions (0 = off)", importance="optional"),
    ),
    advanced_fields=(
        Section("Genomic windows and position-level significance"),
        FieldSpec(name="scoring_window_size_bp", help="Genomic window size in bp used by the overlap reports; default 1000000.", label="Genomic window size (bp)", importance="optional"),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
