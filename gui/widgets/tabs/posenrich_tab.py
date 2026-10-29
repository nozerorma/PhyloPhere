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
        FieldSpec(name="accumulation_fdr", help="Benjamini-Hochberg FDR cutoff for the report tables of levels 1 and 2 (Accumulation); default 0.1.", label="FDR threshold of the report tables", importance="default"),
        FieldSpec(name="posenrich_domains", help="Include Pfam domain annotations; generated from productive alignments when domain_variability_file is blank.", label="Pfam domains", kind="bool", importance="optional"),
        FieldSpec(name="posenrich_ucr", help="Include UCR annotations; generated from productive alignments when ucr_positions_file is blank.", label="Conserved regions (UCR)", kind="bool", importance="optional"),
        FieldSpec(name="fubar_sites_file", help="Per-site FUBAR table from ortholog_characterizator; this pipeline does not generate it.", label="FUBAR sites file (optional)", kind="path_file", importance="optional"),
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
                  help="Standalone observed positives of the universe U (filtered_discovery.tsv); blank uses this run's Post-processing output."),
        FieldSpec(name="accumulation_background_input", label="Standalone cleaned background", kind="path_file", importance="optional",
                  help="Standalone gene list of U (cleaned background); blank uses this run's cleaned background."),
        Section("Annotation sources"),
        FieldSpec(name="domain_variability_file", help="Precomputed Pfam domain variability table; blank generates it, and a supplied table needs an explicit coordinate_space.", label="Domain variability file", kind="path_file", importance="optional"),
        FieldSpec(name="ucr_positions_file", help="Precomputed UCR positions table; blank generates sliding UCR windows, and a supplied table needs an explicit coordinate_space.", label="UCR positions file", kind="path_file", importance="optional"),
        FieldSpec(name="domain_ref_species", help="Reference species for the domain-variability hmmscan; default is the run's reference species, else Homo_sapiens.", label="Domain reference species", placeholder="Homo_sapiens", importance="optional"),
        FieldSpec(name="pfam_cache_dir", help="Directory holding the Pfam-A.hmm and clans files, downloaded once and reused; default ~/.cache/phylophere/pfam.", label="Pfam cache directory", kind="path_dir", importance="optional"),
        Section("Level 2 and proximity"),
        FieldSpec(name="posenrich_collections", label="Collections", placeholder="pfam_domains,ucr_core", importance="optional",
                  help="Comma-separated level 2 collection ids; blank uses the collections informative inside a protein."),
        FieldSpec(name="posenrich_min_genes", help="Minimum genes of U a term needs to be tested in level 2; default 2, so single-gene terms are never tested.", label="Minimum genes per term", placeholder="2", importance="optional"),
        FieldSpec(name="posenrich_holdout", help="Validation mode: this many null cycles stand in for the observed run and only holdout_fpr.tsv is written; 0 = off.", label="Hold-out null cycles (validation, 0 = off)", placeholder="0", importance="optional"),
        FieldSpec(name="posenrich_radios", help="Radii of the protein proximity window, in reference residues, comma-separated; default 0,5,10.", label="Proximity radii (residues)", placeholder="0,5,10", importance="optional"),
        Section("Memory"),
        FieldSpec(name="posenrich_within_gene_memory", help="Memory for level 2, which reads one block of the null at a time; blank uses 16 GB.", label="Level 2", placeholder="16 GB", importance="optional"),
        FieldSpec(name="position_registry_memory", help="Memory for the full normalized position registry; blank sizes it from the background.", label="Registry", placeholder="64 GB", importance="optional"),
        FieldSpec(name="position_module_memory", help="Memory for the full positional tables; blank sizes it from the registry.", label="Components", placeholder="32 GB", importance="optional"),
    ),
)


class PosenrichTab(ModuleTabWidget):
    def __init__(self, config: PosenrichConfig, parent=None):
        super().__init__(SPEC, config, parent)
