#!/usr/bin/env python3
# posenrich_tab.py — POSENRICH module tab (position-wise analysis; Accumulation is its level 1).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
PosenrichTab: field specification of POSENRICH (conf/posenrich.config).

The tab edits the posenrich section of the project (PosenrichConfig). Its enable switch is
--posenrich, the gate of Accumulation (level 1: gene versus its share of the position universe)
and of the position registry with its components (which need SCORING). The PP2A position sets
(posenrich_custom_file) stay on the Enrichment tab.

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
        "tests gene by gene whether CAAS positions accumulate on a gene more than its share of the position "
        "universe U predicts (exact hypergeometric test per gene and direction, Benjamini-Hochberg over the genes "
        "with observed CAAS); its output is published under posenrich/accumulation. Level 2 asks whether a term "
        "(domain, clan, conserved region, site list) concentrates the detections inside each gene, beyond what the "
        "gene already accumulates; its statistic is calibrated with the null cycles. The registry, "
        "characterization and proximity components and report 14 follow; they need Scoring.\n\n"
        "ℹ️ Universe U: the positions CAAStools tested (in genes that survive post-processing) that score "
        "in at least one cycle of the filtered permulation null Scoring also uses, plus the observed CAAS "
        "positions. Levels 1 and 2 read the same universe table. Level 1 treats the positions of U as "
        "equivalent, which understates the variance of genes with many positions in U: read its results as a "
        "ranking, not as a list of significant genes.\n\n"
        "⚠️ Level 2 has not been validated at scale (10,000 null cycles, several term sets): treat its p and q "
        "values as provisional. The proximity control is not validated either."
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
            name="posenrich_min_nonzero_cycles",
            label="Minimum null cycles with signal per gene",
            placeholder="30",
            importance="optional",
            help="A gene enters the statistic of a term only when its residual is nonzero in at least this many null "
                 "cycles (stability of its standard deviation). A heuristic rule, fixed before the observed results.",
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
