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
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

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
        FieldSpec(name="vep_primateai_db", label="PrimateAI-3D database", kind="path_file"),
        FieldSpec(name="vep_map_dir", label="Per-gene MAP directory", kind="path_dir"),
        FieldSpec(name="cosmic_db", label="COSMIC database", kind="path_file"),
    ),
)


class VepTab(ModuleTabWidget):
    def __init__(self, config: VepConfig, parent=None):
        super().__init__(SPEC, config, parent)
