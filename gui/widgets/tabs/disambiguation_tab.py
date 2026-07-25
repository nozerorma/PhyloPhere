#!/usr/bin/env python3
# disambiguation_tab.py — Disambiguation module tab (bundles Post-processing sub-section).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Post-processing (conf/ct_postproc.config) stays bundled here rather than its own
tab: the reference scripts only ever expose it as a bare --ct_postproc boolean
nested inside RUN_DISAMBIGUATION's block, and matching that keeps the GUI aligned
with the scripts it replaces (see implementation plan §5, "Decision").
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import DisambiguationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="Disambiguation",
    blurb=(
        "Resolves CAAS convergence direction via ancestral state reconstruction (ASR), "
        "then (if Post-processing is on) filters/clusters the disambiguated hits."
    ),
    disclaimer=(
        "Accumulation and Scoring need this module's output. Supply "
        "signification_from/disambiguation_input/disambiguation_dir/background_input below."
    ),
    essential_fields=(
        FieldSpec(
            name="ct_disambig_asr_mode",
            label="ASR mode",
            kind="choice",
            choices=("precomputed",),
        ),
        FieldSpec(name="ct_disambig_asr_cache_dir", label="ASR cache directory", kind="path_dir"),
        FieldSpec(name="ct_postproc", label="Run Post-processing (--ct_postproc)", kind="bool"),
        FieldSpec(
            name="caas_postproc_mode",
            label="Post-processing mode",
            kind="choice",
            choices=("filter", "exploratory"),
        ),
        FieldSpec(name="filter_minlen", label="Cluster min length"),
        FieldSpec(name="filter_maxcaas", label="Cluster max CAAS value"),
        FieldSpec(
            name="gene_filter_mode",
            label="Gene filter mode",
            kind="choice",
            choices=("none", "extreme", "dubious", "both"),
        ),
    ),
    fallback_fields=(
        FieldSpec(name="signification_from", label="Signification output", kind="path_dir"),
        FieldSpec(name="disambiguation_input", label="Disambiguation input file", kind="path_file"),
        FieldSpec(name="disambiguation_dir", label="Disambiguation directory", kind="path_dir"),
        FieldSpec(name="background_input", label="Background input", kind="path_file"),
    ),
)


class DisambiguationTab(ModuleTabWidget):
    def __init__(self, config: DisambiguationConfig, parent=None):
        super().__init__(SPEC, config, parent)
