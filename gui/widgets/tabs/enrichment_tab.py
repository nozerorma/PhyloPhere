#!/usr/bin/env python3
# enrichment_tab.py — Enrichment module tab (bundles FCS, STRING/DOMINO, COMPARE, POSENRICH).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
POSENRICH is bundled here rather than its own tab, matching the reference scripts'
RUN_ENRICHMENT / RUN_POSENRICH pairing and conf/enrichment.config, which covers
both FCS/STRING gene-set enrichment and position-wise enrichment in one file.

Two independent STRING-related toggles, both real and both worth keeping distinct
rather than merging into one checkbox:
  - scoring_string: the centralized DOMINO-based PPI run + cross-module COMPARE
    report, driven from workflows/enrichment.nf.
  - string: each of RER/FADE/Accumulation can *also* run its own STRING report
    directly from its own workflow file (workflows/rerconverge.nf, fade.nf,
    ct_accumulation.nf) — publishing under that module's own output tree instead
    of string/.
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import EnrichmentConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Enrichment",
    blurb=(
        "Gene-set enrichment (FCS ranked-Wilcoxon + optional STRING/DOMINO + cross-"
        "module COMPARE) and, if POSENRICH is on, position-wise Fisher-exact "
        "enrichment against domain/UCR/FUBAR annotations."
    ),
    disclaimer="This is the final module in the pipeline — nothing downstream depends on it.",
    essential_fields=(
        FieldSpec(name="gmt_dir", label="GMT pathway directory", kind="path_dir"),
        FieldSpec(name="string_db_dir", label="STRING database directory", kind="path_dir"),
        FieldSpec(name="string", label="Run per-module STRING (RER/FADE/Accumulation)", kind="bool"),
        FieldSpec(name="scoring_ami", label="Run centralized DOMINO AMI + STRING PPI + COMPARE", kind="bool"),
        FieldSpec(name="posenrich_enabled", label="Run POSENRICH", kind="bool"),
        FieldSpec(name="egg_members_file", label="eggNOG members file", kind="path_file"),
        FieldSpec(name="egg_annotations_file", label="eggNOG annotations file", kind="path_file"),
        FieldSpec(name="cosmic_db", label="COSMIC database (also used by VEP)", kind="path_file"),
        FieldSpec(
            name="domain_variability_file", label="Domain variability file", kind="path_file"
        ),
        FieldSpec(name="ucr_positions_file", label="UCR positions file", kind="path_file"),
        FieldSpec(name="fubar_sites_file", label="FUBAR sites file", kind="path_file"),
    ),
    advanced_fields=(
        FieldSpec(name="fcs_min_genes", label="FCS minimum genes per set"),
        FieldSpec(name="fcs_fdr", label="FCS FDR threshold"),
        FieldSpec(name="fcs_pperm_thr", label="FCS permulation p threshold"),
        FieldSpec(name="fcs_top_n", label="FCS top-N leading edge"),
        FieldSpec(name="rer_permulation_enrichment", label="RER permulation-excess null", kind="bool"),
        # NOTE: CAAS's own permulation-excess null is controlled from the CAAS tab
        # (caas_permulation_enrichment) — one param, gates both CT's own
        # CAAS_PERMULATION subworkflow and this report's null, so it has one home.

        FieldSpec(name="string_species", label="STRING species (NCBI taxid)"),
        FieldSpec(name="string_required_score", label="STRING required combined score"),

        FieldSpec(name="domino_network_score_thr", label="DOMINO network score threshold"),
        FieldSpec(name="domino_slice_thr", label="DOMINO slice threshold"),
        FieldSpec(name="domino_module_thr", label="DOMINO module significance threshold"),
        FieldSpec(name="scoring_compare_fdr", label="COMPARE report FDR threshold"),
        FieldSpec(name="scoring_compare_top_n", label="COMPARE report top-N"),

        FieldSpec(name="posenrich_min_size", label="POSENRICH min set size"),
        FieldSpec(name="posenrich_max_size", label="POSENRICH max set size (0 = uncapped)"),
        FieldSpec(name="posenrich_padj_thr", label="POSENRICH adjusted p threshold"),
        FieldSpec(name="posenrich_fold_thr", label="POSENRICH fold-enrichment threshold"),
    ),
)


class EnrichmentTab(ModuleTabWidget):
    def __init__(self, config: EnrichmentConfig, parent=None):
        super().__init__(SPEC, config, parent)
