#!/usr/bin/env python3
# enrichment_tab.py — Enrichment module tab (bundles FCS, STRING/DOMINO, COMPARE).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
EnrichmentTab: field specification of the Enrichment module (conf/enrichment.config).

One tab covers gene-set enrichment (FCS, STRING/DOMINO, COMPARE). The position-wise
POSENRICH parameters live in posenrich_tab.py (conf/posenrich.config).

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
        FieldSpec(
            name="fcs_enabled",
            label="Run FCS gene-set enrichment",
            kind="bool",
            importance="optional",
            help="Runs FCS ranked-Wilcoxon gene-set enrichment on the ranked gene list; default on.",
        ),
        FieldSpec(
            name="auto_fetch_gmt",
            label="Also download the current GO and WikiPathways gene sets",
            kind="bool",
            importance="optional",
            help="Downloads current GO and WikiPathways gene sets and adds them to the versioned ones; applies only when the GMT directory is blank; default off.",
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
        Section("STRING DB Resources"),
        FieldSpec(
            name="string_db_dir",
            label="STRING database directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "Pre-downloaded STRING files, checked before any download; blank downloads and caches them."
            ),
        ),
        FieldSpec(
            name="string_species",
            label="STRING species (NCBI taxid; blank = reference species taxid)",
            importance="default",
            help="NCBI taxid of the STRING species; blank uses the reference species taxid (default 9606, human).",
        ),
    ),
    advanced_fields=(
        Section("Caching"),
        FieldSpec(
            name="string_cache_dir",
            label="STRING cache directory",
            importance="optional",
            help=(
                "Cache of downloaded STRING files; blank uses ~/.cache/phylophere/string. Set it to share a cache across users."
            ),
        ),
        Section("FCS (Wilcoxon-AUC) gene-set enrichment parameters"),
        FieldSpec(
            name="fcs_min_genes",
            label="FCS minimum genes per set",
            importance="default",
            help="Minimum number of genes in a gene set for FCS testing; default 10.",
        ),
        FieldSpec(
            name="fcs_max_genes",
            label="FCS maximum genes per set (0 = uncapped)",
            importance="default",
            help="Maximum number of genes in a gene set for FCS testing; 0 = no upper limit; default 0.",
        ),
        FieldSpec(
            name="fcs_fdr",
            label="FCS FDR threshold",
            importance="default",
            help="Benjamini-Hochberg threshold on the permutation p of each set (q.perm), per ranking and database; default 0.15.",
        ),
        FieldSpec(
            name="fcs_pperm_thr",
            label="FCS Redlich permulation p threshold",
            importance="default",
            help="Raw permutation p threshold of the Redlich analytic rule, independent of the q.perm significance; default 0.025.",
        ),
        FieldSpec(
            name="fcs_analytic_levels",
            label="FCS analytic rule (Redlich): BH levels, comma-separated",
            importance="default",
            help="Benjamini-Hochberg levels of the analytic rule of Redlich et al., comma-separated; reported as redlich_level; default 0.05,0.1,0.15.",
        ),
        FieldSpec(
            name="fcs_cluster_threshold",
            label="FCS term blocks: minimum Jaccard inside a block",
            importance="default",
            help="Minimum Jaccard overlap of the full gene sets needed to place two terms in one block (complete linkage); default 0.2.",
        ),
        FieldSpec(
            name="fcs_top_n",
            label="FCS top-N leading edge",
            importance="optional",
            help="Number of top gene sets shown in the report tables; default 20.",
        ),
        FieldSpec(
            name="fcs_batch_size",
            label="FCS GMTs per task",
            importance="optional",
            help="GMT files per FCS compute task; batching does not change the adjusted p-values; default 4.",
        ),
        FieldSpec(
            name="fcs_block_batch_genes",
            label="GENE_P genes per task",
            importance="optional",
            help="Genes per GENE_P task; memory grows with genes per task times total cycles; default 1000.",
        ),
        Section("DOMINO active module identification thresholds"),
        FieldSpec(
            name="scoring_ami",
            label="Run centralized DOMINO AMI + STRING PPI + COMPARE",
            kind="bool",
            importance="optional",
            help="Runs the DOMINO active-module (AMI) analysis with STRING PPI and the COMPARE report; scoring_string is a legacy alias; default on.",
        ),
        FieldSpec(
            name="domino_network_score_thr",
            label="DOMINO network score threshold",
            importance="default",
            help="Minimum STRING combined score (0-1000) for DOMINO network edges; default 700.",
        ),
        FieldSpec(
            name="domino_slice_thr",
            label="DOMINO slice threshold",
            importance="default",
            help="DOMINO network slice relevance threshold (0-1); default 0.3.",
        ),
        FieldSpec(
            name="domino_module_thr",
            label="DOMINO module significance threshold",
            importance="default",
            help="DOMINO module significance threshold; default 0.05.",
        ),
        FieldSpec(
            name="ami_max_blocks",
            label="AMI term blocks: blocks shown per module",
            importance="default",
            help="Term blocks shown in the label and summary of each AMI module; the tables keep all blocks; default 8.",
        ),
        FieldSpec(
            name="ami_block_threshold",
            label="AMI term blocks: minimum Jaccard of module genes",
            importance="default",
            help="Minimum Jaccard overlap of module genes between two terms for them to share an AMI block (complete linkage); default 0.2.",
        ),
        FieldSpec(
            name="publish_domino_intermediates",
            label="Publish DOMINO network.sif/modules/edge scores (debug)",
            kind="bool",
            importance="optional",
            help="Publishes DOMINO network.sif, module and edge-score files for debugging; default off.",
        ),
        FieldSpec(
            name="scoring_compare_fdr",
            label="COMPARE report FDR threshold",
            importance="default",
            help="FDR cutoff of the cross-tool COMPARE report; default 0.15.",
        ),
        FieldSpec(
            name="scoring_compare_top_n",
            label="COMPARE report top-N",
            importance="optional",
            help="Number of top genes highlighted in the COMPARE report; default 20.",
        ),
    ),
)


class EnrichmentTab(ModuleTabWidget):
    def __init__(self, config: EnrichmentConfig, parent=None):
        super().__init__(SPEC, config, parent)
