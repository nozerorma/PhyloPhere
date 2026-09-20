#!/usr/bin/env python3
# enrichment_tab.py — Enrichment module tab (bundles FCS, STRING/DOMINO, COMPARE, POSENRICH).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
POSENRICH is bundled here rather than its own tab, matching the reference scripts'
RUN_ENRICHMENT / RUN_POSENRICH pairing and conf/enrichment.config, which covers
both FCS/STRING gene-set enrichment and position-wise enrichment in one file.

scoring_string (alias scoring_ami) is the single AMI toggle: the centralized
DOMINO-based AMI run + cross-module COMPARE report, driven from
workflows/enrichment.nf, covering CAAS/FADE/RER in one unified
13.AMI_analysis.Rmd. RER/FADE/Accumulation's own gene lists are always computed
automatically whenever those tools run — no separate --ami flag; the old
standalone per-module AMI reports (one HTML per tool) were retired since they
never produced usable output.
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import EnrichmentConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Enrichment",
    blurb=(
        "Gene-set enrichment (FCS ranked-Wilcoxon + optional STRING/DOMINO + cross-"
        "module COMPARE) and, if POSENRICH is on, position-wise Fisher-exact "
        "enrichment against domain/UCR/FUBAR annotations."
    ),
    disclaimer="This is the final module in the pipeline — nothing downstream depends on it.",
    essential_fields=(
        Section("STRING DB & GMT Resources"),
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
            name="gmt_dir",
            label="GMT pathway directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "Custom GMT gene-set files. Leave blank to auto-fetch the default "
                "GO Biological Process/Molecular Function + Reactome + WikiPathways "
                "set (fetches current copies; falls back to vendored assets/gmt/ "
                "if offline)."
            ),
        ),
        # Determines which species' STRING/GO background is queried (default 9606 =
        # human) — changing it changes what the enrichment background *means*
        # scientifically, not just where files are cached.
        FieldSpec(name="string_species", label="STRING species (NCBI taxid)", importance="default"),
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
        FieldSpec(name="fcs_pperm_thr", label="FCS permulation p threshold", importance="default"),
        FieldSpec(name="fcs_top_n", label="FCS top-N leading edge", importance="optional"),
        FieldSpec(name="fcs_batch_size", label="FCS GMTs per task", importance="optional"),
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
        FieldSpec(
            name="publish_domino_intermediates",
            label="Publish DOMINO network.sif/modules/edge scores (debug)",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(name="scoring_compare_fdr", label="COMPARE report FDR threshold", importance="default"),
        FieldSpec(name="scoring_compare_top_n", label="COMPARE report top-N", importance="optional"),
        # Borderline default/optional: gates whether the concordance null chunk is
        # computed at all, but disabling it doesn't change any existing result's
        # validity — it just skips an extra corroborating statistical test, same
        # spirit as scoring_stress below.
        FieldSpec(
            name="comparison_perm_null",
            label="COMPARE: CAAS x RER concordance null",
            kind="bool",
            importance="optional",
        ),
        FieldSpec(
            name="comparison_perm_stat",
            label="COMPARE concordance statistic (spearman | topk_overlap)",
            kind="choice",
            choices=("spearman", "topk_overlap"),
            importance="default",
        ),
        FieldSpec(name="comparison_perm_topk", label="COMPARE concordance top-k fraction", importance="default"),
        Section("POSENRICH position-wise enrichment parameters"),
        FieldSpec(name="posenrich_enabled", label="Run POSENRICH", kind="bool", importance="optional"),
        FieldSpec(name="posenrich_min_size", label="POSENRICH min set size", importance="default"),
        FieldSpec(name="posenrich_max_size", label="POSENRICH max set size (0 = uncapped)", importance="default"),
        FieldSpec(name="posenrich_padj_thr", label="POSENRICH adjusted p threshold", importance="default"),
        FieldSpec(name="posenrich_p_perm_thr", label="POSENRICH CAAS-null p.perm threshold", importance="default"),
        FieldSpec(
            name="posenrich_batch_size",
            label="POSENRICH GMTs per task (1 = no batching)",
            importance="optional",
        ),
        FieldSpec(name="domain_variability_file", label="Domain variability file", kind="path_file", importance="optional"),
        FieldSpec(name="ucr_positions_file", label="UCR positions file", kind="path_file", importance="optional"),
        # validate.py's Enrichment section requires this whenever POSENRICH is
        # enabled (require(enrichment.fubar_sites_file, ...)) — it's the one
        # POSENRICH input with no auto-generation fallback.
        FieldSpec(
            name="fubar_sites_file",
            label="FUBAR sites file",
            kind="path_file",
            importance="required",
            help="Cannot be generated in-house: HyPhy's per-site FUBAR fit needs "
                 "the full phylogeny + codon alignment + MCMC/VB inference, not "
                 "just the alignment plus a public DB. See "
                 "github.com/nozerorma/ortholog_characterizator. Required for "
                 "POSENRICH.",
        ),
        FieldSpec(
            name="egg_members_file",
            label="eggNOG members file (optional)",
            kind="path_file",
            importance="optional",
            help=(
                "eggNOG5 Primates orthogroup members. Leave blank to auto-fetch "
                "(falls back to the vendored human-subset copy in assets/eggnog/ "
                "if offline)."
            ),
        ),
        FieldSpec(
            name="egg_annotations_file",
            label="eggNOG annotations file (optional)",
            kind="path_file",
            importance="optional",
            help=(
                "eggNOG5 Primates orthogroup annotations, paired with the members "
                "file above. Leave blank to auto-fetch alongside it."
            ),
        ),
    ),
)


class EnrichmentTab(ModuleTabWidget):
    def __init__(self, config: EnrichmentConfig, parent=None):
        super().__init__(SPEC, config, parent)
