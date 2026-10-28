#!/usr/bin/env python3
# enrichment_tab.py — Enrichment module tab (bundles FCS, STRING/DOMINO, COMPARE).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
EnrichmentTab: field specification of the Enrichment module (conf/enrichment.config).

One tab covers gene-set enrichment (FCS, STRING/DOMINO, COMPARE). The position-wise
POSENRICH parameters live in posenrich_tab.py (conf/posenrich.config); the PP2A position sets
(posenrich_custom_file) stay here.

scoring_string (with its alias scoring_ami) is the single toggle for the
DOMINO-based active-module identification (AMI) run, which workflows/enrichment.nf
enables when either parameter is set; its report is 13.AMI_analysis.Rmd.

Imported by: gui/widgets/main_window.py
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import EnrichmentConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Enrichment",
    blurb=(
        "Gene-set enrichment (FCS ranked-Wilcoxon + optional STRING/DOMINO + cross-"
        "module COMPARE). Position-wise enrichment has its own POSENRICH tab."
    ),
    disclaimer="This is the final module in the pipeline — nothing downstream depends on it.",
    essential_fields=(
        Section("FCS & Pathway Resources"),
        FieldSpec(name="fcs_enabled", label="Run FCS gene-set enrichment", kind="bool", importance="optional"),
        FieldSpec(
            name="auto_fetch_gmt",
            label="Also download the current GO and WikiPathways gene sets",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(
            name="gmt_dir",
            label="GMT pathway directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "Custom GMT gene-set files. Leave blank to use the curated default "
                "pathway sets in subworkflows/ENRICHMENT/dat/ (or auto-fetch if "
                "enabled)."
            ),
        ),
        Section("Position annotations — PP2A"),
        FieldSpec(
            name="posenrich_custom_file",
            label="PP2A / custom position sets",
            kind="path_file",
            importance="optional",
            help=(
                "Validated position sets, including PP2A phosphosites. Select the "
                "pp2a_custom.tsv for the alignment reference used by this run "
                "(PPGA or marine). The neoplasia and marine templates supply the "
                "corresponding cluster path. Reload the template to load it into "
                "an existing form."
            ),
        ),
        Section("STRING DB Resources"),
        FieldSpec(
            name="string_db_dir",
            label="STRING database directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "Pre-downloaded STRING files, checked before downloading. Leave "
                "blank to auto-download+cache instead (see STRING cache directory, "
                "Advanced) — no need to pre-populate this by hand."
            ),
        ),
        FieldSpec(
            name="string_species",
            label="STRING species (NCBI taxid; blank = reference species taxid)",
            importance="default",
        ),
    ),
    advanced_fields=(
        Section("Caching"),
        FieldSpec(
            name="string_cache_dir",
            label="STRING cache directory",
            importance="optional",
            help=(
                "Where auto-downloaded STRING links/info files are cached across "
                "runs. Defaults to ~/.cache/phylophere/string — override only to "
                "share a cache location across users/clusters."
            ),
        ),
        Section("FCS (Wilcoxon-AUC) gene-set enrichment parameters"),
        FieldSpec(name="fcs_min_genes", label="FCS minimum genes per set", importance="default"),
        FieldSpec(name="fcs_max_genes", label="FCS maximum genes per set (0 = uncapped)", importance="default"),
        FieldSpec(name="fcs_fdr", label="FCS FDR threshold", importance="default"),
        FieldSpec(name="fcs_pperm_thr", label="FCS Redlich permulation p threshold", importance="default"),
        FieldSpec(name="fcs_analytic_levels", label="FCS analytic rule (Redlich): BH levels, comma-separated", importance="default"),
        FieldSpec(name="fcs_cluster_threshold", label="FCS term blocks: minimum Jaccard inside a block", importance="default"),
        FieldSpec(name="fcs_top_n", label="FCS top-N leading edge", importance="optional"),
        FieldSpec(name="fcs_batch_size", label="FCS GMTs per task", importance="optional"),
        FieldSpec(name="fcs_block_batch_genes", label="GENE_P genes per task", importance="optional"),
        Section("DOMINO active module identification thresholds"),
        FieldSpec(
            name="scoring_ami",
            label="Run centralized DOMINO AMI + STRING PPI + COMPARE",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(name="domino_network_score_thr", label="DOMINO network score threshold", importance="default"),
        FieldSpec(name="domino_slice_thr", label="DOMINO slice threshold", importance="default"),
        FieldSpec(name="domino_module_thr", label="DOMINO module significance threshold", importance="default"),
        FieldSpec(name="ami_max_blocks", label="AMI term blocks: blocks shown per module", importance="default"),
        FieldSpec(name="ami_block_threshold", label="AMI term blocks: minimum Jaccard of module genes", importance="default"),
        FieldSpec(
            name="publish_domino_intermediates",
            label="Publish DOMINO network.sif/modules/edge scores (debug)",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(name="scoring_compare_fdr", label="COMPARE report FDR threshold", importance="default"),
        FieldSpec(name="scoring_compare_top_n", label="COMPARE report top-N", importance="optional"),
    ),
)


class EnrichmentTab(ModuleTabWidget):
    def __init__(self, config: EnrichmentConfig, parent=None):
        super().__init__(SPEC, config, parent)
