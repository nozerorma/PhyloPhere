#!/usr/bin/env python3
# vep_tab.py — VEP module tab (variant effect / pathogenicity annotation).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
VepTab: field specification of the VEP module (variant effect annotation; conf/vep.config).

cosmic_db (conf/vep.config) is the COSMIC database path. It is one pipeline
parameter with two readers: workflows/vep.nf builds its COSMIC channel from it, and
workflows/enrichment.nf hands it to POSENRICH. Its only GUI field is on this tab,
so it is editable only while the VEP module is enabled. scoring_vep_cosmic
(Scoring's standalone fallback for a precomputed COSMIC scores TSV,
conf/scoring.config) is a different parameter and lives on the Precomputed Run tab.

Imported by: gui/widgets/main_window.py
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import VepConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="VEP",
    blurb=(
        "Annotates CAAS hits with PrimateAI-3D pathogenicity scores and COSMIC "
        "cancer-mutation overlap."
    ),
    disclaimer=(
        "Scoring can fall back to this module's PrimateAI-3D/COSMIC score outputs "
        "when it's off — check 'Use precomputed VEP output' on the Precomputed Run tab."
    ),
    essential_fields=(
        Section("Variant effect prediction databases"),
        FieldSpec(
            name="vep_primateai_db",
            label="PrimateAI-3D database (optional)",
            kind="path_file",
            help="PrimateAI-3D pathogenicity table you supply; it is Illumina-licensed, so it is neither bundled nor downloaded.",
            importance="optional",
        ),
        FieldSpec(
            name="cosmic_db",
            label="COSMIC database (optional)",
            kind="path_file",
            help="COSMIC mutation database you supply; it is neither bundled nor downloaded, as with PrimateAI-3D.",
            importance="optional",
        ),
        # gui/generation/validate.py requires disambiguation.caas_map_dir whenever
        # vep is enabled: the per-gene MAP directory is a post-processing
        # parameter that VEP reuses.
        Section("Per-gene MAP directory: set under Disambiguation, post-processing (required by VEP)"),
        Section("Ensembl VEP (independent of the databases above)"),
        FieldSpec(
            name="vep_ensembl",
            label="Enable Ensembl VEP consequence annotation",
            kind="bool",
            help="Run the official Ensembl VEP CLI for consequence prediction, independent of the two databases above.",
            importance="optional",
        ),
    ),
    advanced_fields=(
        FieldSpec(
            name="vep_cache_dir",
            label="Ensembl VEP cache directory",
            kind="path_dir",
            help="Ensembl VEP offline cache; blank uses ~/.cache/phylophere/vep/<species>_<assembly>, filled on first use. Set only to share a cache.",
            importance="default",
        ),
        # Rated "default" although borderline: validate.py does not check it and it
        # has a working default ("homo_sapiens"/"GRCh38"), but a mismatch with the
        # actual reference data gives wrong annotations without an error.
        FieldSpec(name="vep_species", help="Species passed to Ensembl VEP; default homo_sapiens.", label="Ensembl VEP species", importance="default"),
        FieldSpec(name="vep_assembly", help="Assembly passed to Ensembl VEP; must match your reference data, since a mismatch gives wrong annotations silently. Default GRCh38.", label="Ensembl VEP assembly", importance="default"),
        FieldSpec(
            name="vep_mode",
            label="Ensembl VEP source",
            kind="choice",
            choices=("auto", "database", "offline"),
            help="auto: Ensembl database when vep_db_host is reachable, else offline cache; database: no download, slow for many variants; offline: local cache only.",
            importance="default",
        ),
        FieldSpec(name="vep_db_host", help="Ensembl MySQL host queried in database mode; default ensembldb.ensembl.org.", label="Ensembl VEP database host", importance="optional"),
        FieldSpec(name="vep_db_port", help="Port of the Ensembl MySQL host; default 3306.", label="Ensembl VEP database port", importance="optional"),
    ),
)


class VepTab(ModuleTabWidget):
    def __init__(self, config: VepConfig, parent=None):
        super().__init__(SPEC, config, parent)
