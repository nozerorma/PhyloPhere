#!/usr/bin/env python3
# enrichment_tab.py — Enrichment module tab (bundles POSENRICH position-wise enrichment).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
POSENRICH is bundled here rather than its own tab, matching the reference scripts'
RUN_ENRICHMENT / RUN_POSENRICH pairing and conf/enrichment.config, which covers
both FCS/STRING gene-set enrichment and position-wise enrichment in one file.
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import EnrichmentConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Enrichment",
    blurb=(
        "Gene-set enrichment (FCS ranked-Wilcoxon + optional STRING network) and, if "
        "POSENRICH is on, position-wise Fisher-exact enrichment against domain/UCR/FUBAR "
        "annotations."
    ),
    disclaimer="This is the final module in the pipeline — nothing downstream depends on it.",
    essential_fields=(
        FieldSpec(name="gmt_dir", label="GMT pathway directory", kind="path_dir"),
        FieldSpec(name="string_db_dir", label="STRING database directory", kind="path_dir"),
        FieldSpec(name="posenrich_enabled", label="Run POSENRICH", kind="bool"),
        FieldSpec(name="egg_members_file", label="eggNOG members file", kind="path_file"),
        FieldSpec(name="egg_annotations_file", label="eggNOG annotations file", kind="path_file"),
        FieldSpec(name="cosmic_db", label="COSMIC database", kind="path_file"),
        FieldSpec(
            name="domain_variability_file", label="Domain variability file", kind="path_file"
        ),
        FieldSpec(name="ucr_positions_file", label="UCR positions file", kind="path_file"),
        FieldSpec(name="fubar_sites_file", label="FUBAR sites file", kind="path_file"),
    ),
    fallback_fields=(
        FieldSpec(
            name="enrichment_background_input", label="Background input", kind="path_file"
        ),
        FieldSpec(
            name="enrichment_gene_lists_input", label="Gene lists input", kind="path_file"
        ),
        FieldSpec(
            name="accumulation_enrichment_gene_lists_input",
            label="Accumulation gene lists input",
            kind="path_file",
        ),
        FieldSpec(
            name="posenrich_background_file", label="POSENRICH background file", kind="path_file"
        ),
    ),
)


class EnrichmentTab(ModuleTabWidget):
    def __init__(self, config: EnrichmentConfig, parent=None):
        super().__init__(SPEC, config, parent)
