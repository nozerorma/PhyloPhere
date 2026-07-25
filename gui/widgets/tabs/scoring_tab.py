#!/usr/bin/env python3
# scoring_tab.py — Scoring module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import ScoringConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Scoring",
    blurb=(
        "Combines CAAS, RERconverge, FADE, Accumulation, and VEP outputs into a "
        "composite per-gene / per-position score."
    ),
    disclaimer="Enrichment needs Scoring's gene lists/background when this is off.",
    essential_fields=(
        FieldSpec(name="scoring_stress", label="Run stress-enrichment analysis", kind="bool"),
        FieldSpec(name="scoring_window_size_bp", label="Genomic window size (bp)"),
        FieldSpec(name="gene_ensembl_file", label="Gene-Ensembl mapping file", kind="path_file"),
        FieldSpec(
            name="scoring_fade_site_top",
            label="FADE site-level fallback — top (if FADE off)",
            kind="path_file",
        ),
        FieldSpec(
            name="scoring_fade_site_bottom",
            label="FADE site-level fallback — bottom (if FADE off)",
            kind="path_file",
        ),
    ),
    fallback_fields=(
        FieldSpec(name="scoring_postproc_input", label="Post-processing input", kind="path_file"),
        FieldSpec(name="scoring_accum_dir", label="Accumulation directory", kind="path_dir"),
        FieldSpec(name="scoring_vep_primateai", label="VEP PrimateAI scores", kind="path_file"),
        FieldSpec(name="scoring_background_input", label="Background input", kind="path_file"),
        FieldSpec(name="caas_perms_file", label="CAAS permulation RDS", kind="path_file"),
    ),
)


class ScoringTab(ModuleTabWidget):
    def __init__(self, config: ScoringConfig, parent=None):
        super().__init__(SPEC, config, parent)
