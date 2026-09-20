#!/usr/bin/env python3
# vep_tab.py — VEP module tab (variant effect / pathogenicity annotation).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
cosmic_db (conf/vep.config) is VEP's own COSMIC database path — workflows/vep.nf
reads params.cosmic_db directly to build its cosmic_db_file channel. An earlier
version of this tab instead exposed scoring_vep_cosmic here (Scoring's *own*
standalone fallback for a precomputed COSMIC *scores* TSV, conf/scoring.config),
which VEP never reads — so filling in "COSMIC database" on this tab silently did
nothing for VEP's actual COSMIC annotation. scoring_vep_cosmic now lives on the
Precomputed Run tab where it belongs, alongside the rest of Scoring's fallbacks.
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
            help="Illumina-licensed; redistribution isn't permitted, so this stays "
                 "a file you supply yourself — there's no vendoring or auto-download "
                 "path for it.",
            importance="optional",
        ),
        FieldSpec(
            name="cosmic_db",
            label="COSMIC database (optional)",
            kind="path_file",
            help="License terms for redistribution haven't been confirmed yet, so "
                 "this stays external for now (same treatment as PrimateAI-3D) "
                 "pending that check.",
            importance="optional",
        ),
        # validate.py: require(vep.vep_map_dir, ...) unconditionally whenever
        # vep.enabled — matches this field's own help text ("required whenever
        # any annotation source on this tab is used").
        FieldSpec(
            name="vep_map_dir",
            label="Per-gene MAP directory",
            kind="path_dir",
            help="Cannot be generated in-house: maps each alignment codon column "
                 "to its real hg38 genomic/protein coordinate, which needs the full "
                 "alignment-to-protein pipeline, not just the alignment plus a "
                 "public DB. See github.com/nozerorma/ortholog_characterizator. "
                 "Required whenever any annotation source on this tab is used.",
            importance="required",
        ),
        Section("Ensembl VEP (independent of the databases above)"),
        FieldSpec(
            name="vep_ensembl",
            label="Enable Ensembl VEP consequence annotation",
            kind="bool",
            help="Runs the official Ensembl VEP CLI for consequence prediction — "
                 "works even when neither PrimateAI-3D nor COSMIC is supplied. "
                 "Needs a pre-downloaded offline cache (below).",
            importance="optional",
        ),
    ),
    advanced_fields=(
        # validate.py: require(vep.vep_cache_dir, ...) conditionally, only when
        # vep.vep_ensembl is checked — required within that branch, not globally.
        FieldSpec(
            name="vep_cache_dir",
            label="Ensembl VEP cache directory",
            kind="path_dir",
            help="Local offline VEP cache. Multi-GB, not auto-downloaded — "
                 "populate once with: vep_install -a cf -s <species> -y <assembly> "
                 "-c <this dir> --NO_HTSLIB, then reuse across runs. Required when "
                 "'Enable Ensembl VEP consequence annotation' is checked.",
            importance="required",
        ),
        # Borderline default/optional: not in validate.py and has a working
        # default ("homo_sapiens"/"GRCh38"), but silently mismatching your actual
        # reference data produces wrong annotations rather than an obvious error —
        # treated as "default" so changing it prompts a second look.
        FieldSpec(name="vep_species", label="Ensembl VEP species", importance="default"),
        FieldSpec(name="vep_assembly", label="Ensembl VEP assembly", importance="default"),
    ),
)


class VepTab(ModuleTabWidget):
    def __init__(self, config: VepConfig, parent=None):
        super().__init__(SPEC, config, parent)
