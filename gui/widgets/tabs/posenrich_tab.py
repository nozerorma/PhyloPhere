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
        "Position-wise analysis of the CAAS positions, after Scoring and before Enrichment. Level 1, Accumulation, "
        "tests gene by gene whether a gene holds a larger share of the detections than it does in the null "
        "cycles (Benjamini-Hochberg over the genes with observed CAAS); its output is published under "
        "posenrich/accumulation. Level 2 asks whether a term (domain, clan, conserved region, site list) "
        "concentrates the detections inside each gene more than the rest of the same gene; its statistic is "
        "calibrated with the null cycles. The registry, characterization, overlap and proximity components and "
        "report 14 follow; they need Scoring.\n\n"
        "ℹ️ Universe U: the positions CAAStools tested (in genes that survive post-processing) that score "
        "in at least one cycle of the filtered permulation null Scoring also uses, plus the observed CAAS "
        "positions. Levels 1 and 2 read the same universe table and the cycles of the null. The null floor of p "
        "is 1 / (cycles + 1): read the results as a ranking with a transparent error rate, not as a list of "
        "significant genes.\n\n"
        "⚠️ Level 2 has not been validated at scale (10,000 null cycles, several term sets): treat its p and q "
        "values as provisional. Overlaps and proximity are descriptive, without p-values."
    ),
    disclaimer="no POSENRICH component runs (registry, characterization, level 2, proximity, Accumulation).",
    essential_fields=(
        Section("Universe and Accumulation (level 1)"),
        FieldSpec(name="accumulation_fdr", label="FDR threshold of the report tables", importance="default"),
        FieldSpec(
            name="accumulation_caas_input",
            label="Standalone filtered_discovery.tsv (optional)",
            kind="path_file",
            importance="optional",
            help="Blank: the Post-processing output of the run (live or precomputed).",
        ),
        FieldSpec(
            name="accumulation_background_input",
            label="Standalone cleaned background (optional)",
            kind="path_file",
            importance="optional",
            help="Blank: the cleaned background of the run (live or precomputed). A directory holding cleaned_background_main.txt also works.",
        ),
        Section("Position references"),
        FieldSpec(name="posenrich_domains", label="Include Pfam annotations (blank file = automatic)", kind="bool", importance="optional"),
        FieldSpec(name="domain_variability_file", label="Domain variability file (optional; blank = automatic)", kind="path_file", importance="optional"),
        FieldSpec(
            name="domain_ref_species",
            label="Domain reference species (blank = reference species name)",
            placeholder="Homo_sapiens",
            importance="optional",
            help=(
                "Reference species name in alignment headers used to extract "
                "sequences for Pfam hmmscan (default: reference species name). Specify "
                "this when running non-human/non-primate alignments."
            ),
        ),
        FieldSpec(name="pfam_cache_dir", label="Pfam cache directory (blank = ~/.cache/phylophere/pfam)", kind="path_dir", importance="optional"),
        FieldSpec(name="posenrich_ucr", label="Include UCR annotations (blank file = automatic)", kind="bool", importance="optional"),
        FieldSpec(name="ucr_positions_file", label="UCR positions file (optional; blank = automatic sliding UCR)", kind="path_file", importance="optional"),
        FieldSpec(
            name="fubar_sites_file",
            label="FUBAR sites file (optional)",
            kind="path_file",
            importance="optional",
            help="Cannot be generated in-house: HyPhy's per-site FUBAR fit needs "
                 "the full phylogeny + codon alignment + MCMC/VB inference, not "
                 "just the alignment plus a public DB. See "
                 "github.com/nozerorma/ortholog_characterizator.",
        ),
        Section("Custom position sets"),
        FieldSpec(
            name="posenrich_custom_file",
            label="Custom position sets (CAAS coordinates)",
            kind="path_file",
            importance="optional",
            help=(
                "Any list of positions to characterize and test as terms: a TSV without header with gene, "
                "zero-based CAAS position, term and an optional description. The positions must already be "
                "validated in the CAAS coordinates of the alignment reference of this run; rows outside the MAP "
                "are reported as pending, never as negatives. Level 2 tests every term with members and "
                "non-members in a gene; terms with few members inside the universe U are reported as not "
                "identifiable."
            ),
        ),
        Section("Level 2: concentration of the signal inside genes"),
        FieldSpec(
            name="posenrich_collections",
            label="Collections (blank = those informative inside a protein)",
            placeholder="pfam_domains,ucr_core",
            importance="optional",
            help="Comma-separated collection ids tested by level 2. Blank uses the collections the registry marks as "
                 "informative inside a protein.",
        ),
        FieldSpec(
            name="posenrich_min_genes",
            label="Minimum genes per term",
            placeholder="2",
            importance="optional",
            help="A term is tested only when it has members in at least this many genes of the universe. Fixed by "
                 "the annotation, not by the observed positives; a term in one gene is a test inside that gene.",
        ),
        FieldSpec(
            name="posenrich_holdout",
            label="Hold-out null cycles (validation, 0 = off)",
            placeholder="0",
            importance="optional",
            help="Validation mode: this many null cycles play the observed one and only the false-positive table is written.",
        ),
        FieldSpec(name="posenrich_within_gene_memory", label="Level 2 memory (blank = 16 GB)", placeholder="16 GB", importance="optional"),
        Section("Position registry: characterization, proximity"),
        FieldSpec(name="posenrich_radios", label="Proximity radii (residues, descriptive)", placeholder="0,5,10", importance="optional"),
        FieldSpec(name="position_registry_memory", label="Registry memory (blank = sized from the background)", placeholder="64 GB", importance="optional"),
        FieldSpec(name="position_module_memory", label="Component memory (blank = sized from the registry)", placeholder="32 GB", importance="optional"),
    ),
)


class PosenrichTab(ModuleTabWidget):
    def __init__(self, config: PosenrichConfig, parent=None):
        super().__init__(SPEC, config, parent)
