#!/usr/bin/env python3
# posenrich_tab.py — POSENRICH module tab (position-wise analysis; Accumulation is its level 1).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
PosenrichTab: field specification of POSENRICH (conf/posenrich.config).

The tab edits the posenrich section of the project (PosenrichConfig). Its enable switch is
--posenrich, the gate of Accumulation (level 1: gene versus its share of the position universe)
and of the position registry with its components (which need SCORING).

Imported by: gui/widgets/main_window.py
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import PosenrichConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="POSENRICH",
    blurb=(
        "Position-level analysis of the CAAS detections: gene accumulation, concentration in domains, "
        "regions and site lists within genes, and overlap or proximity to sites. Needs Scoring."
    ),
    disclaimer="no POSENRICH component runs.",
    essential_fields=(
        FieldSpec(name="accumulation_fdr", label="FDR threshold of the report tables", importance="default"),
        FieldSpec(name="posenrich_domains", label="Pfam domains", kind="bool", importance="optional"),
        FieldSpec(name="posenrich_ucr", label="Conserved regions (UCR)", kind="bool", importance="optional"),
        FieldSpec(name="fubar_sites_file", label="FUBAR sites file (optional)", kind="path_file", importance="optional"),
        FieldSpec(
            name="posenrich_custom_file",
            label="Custom position sets",
            kind="path_file",
            importance="optional",
            help="TSV without header: gene, zero-based CAAS position, term, optional description.",
        ),
    ),
    advanced_fields=(
        Section("Universe"),
        FieldSpec(name="accumulation_caas_input", label="Standalone filtered_discovery.tsv", kind="path_file", importance="optional",
                  help="Blank: the Post-processing output of the run."),
        FieldSpec(name="accumulation_background_input", label="Standalone cleaned background", kind="path_file", importance="optional",
                  help="Blank: the cleaned background of the run."),
        Section("Annotation sources"),
        FieldSpec(name="domain_variability_file", label="Domain variability file", kind="path_file", importance="optional"),
        FieldSpec(name="ucr_positions_file", label="UCR positions file", kind="path_file", importance="optional"),
        FieldSpec(name="domain_ref_species", label="Domain reference species", placeholder="Homo_sapiens", importance="optional"),
        FieldSpec(name="pfam_cache_dir", label="Pfam cache directory", kind="path_dir", importance="optional"),
        Section("Level 2 and proximity"),
        FieldSpec(name="posenrich_collections", label="Collections", placeholder="pfam_domains,ucr_core", importance="optional",
                  help="Blank: the collections informative inside a protein."),
        FieldSpec(name="posenrich_min_genes", label="Minimum genes per term", placeholder="2", importance="optional"),
        FieldSpec(name="posenrich_holdout", label="Hold-out null cycles (validation, 0 = off)", placeholder="0", importance="optional"),
        FieldSpec(name="posenrich_radios", label="Proximity radii (residues)", placeholder="0,5,10", importance="optional"),
        Section("Memory"),
        FieldSpec(name="posenrich_within_gene_memory", label="Level 2", placeholder="16 GB", importance="optional"),
        FieldSpec(name="position_registry_memory", label="Registry", placeholder="64 GB", importance="optional"),
        FieldSpec(name="position_module_memory", label="Components", placeholder="32 GB", importance="optional"),
    ),
)


class PosenrichTab(ModuleTabWidget):
    def __init__(self, config: PosenrichConfig, parent=None):
        super().__init__(SPEC, config, parent)
