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

ASR mode: disambiguation_main.py's own argparse only accepts "compute" or
"precomputed" (subworkflows/CT_DISAMBIGUATION/local/disambiguation_main.py) — an
earlier version of this tab only offered "precomputed" as a choice, silently ruling
out live ASR computation.
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
        "signification_from/disambiguation_input/disambiguation_dir/background_input "
        "on the Precomputed Run tab instead."
    ),
    essential_fields=(
        FieldSpec(
            name="ct_disambig_asr_mode",
            label="ASR mode",
            kind="choice",
            choices=("precomputed", "compute"),
        ),
        FieldSpec(name="ct_disambig_asr_cache_dir", label="ASR cache directory", kind="path_dir"),
        FieldSpec(name="asr_robustness", label="Run ASR Robustness diagnostics report", kind="bool"),
        FieldSpec(name="ct_postproc", label="Run Post-processing (--ct_postproc)", kind="bool"),
        FieldSpec(
            name="caas_postproc_mode",
            label="Post-processing mode",
            kind="choice",
            choices=("filter", "exploratory"),
        ),
        FieldSpec(name="filter_minlen", label="Cluster min length (filter mode)"),
        FieldSpec(name="filter_maxcaas", label="Cluster max CAAS value (filter mode)"),
        FieldSpec(
            name="gene_filter_mode",
            label="Gene filter mode",
            kind="choice",
            choices=("none", "extreme", "dubious", "both"),
        ),
    ),
    advanced_fields=(
        FieldSpec(name="ct_disambig_asr_model", label="ASR substitution model"),
        FieldSpec(
            name="ct_disambig_convergence_mode",
            label="Convergence mode",
            kind="choice",
            choices=("focal_clade", "mrca"),
        ),
        FieldSpec(name="ct_disambig_posterior_threshold", label="Posterior probability threshold"),
        FieldSpec(name="ct_disambig_max_tasks_per_child", label="Max tasks per worker child"),
        FieldSpec(name="minlen_values", label="Cluster min length sweep (exploratory mode)"),
        FieldSpec(name="maxcaas_values", label="Cluster max CAAS sweep (exploratory mode)"),
        FieldSpec(name="extreme_threshold", label="Extreme-gene quantile threshold"),
        FieldSpec(name="iqr_multiplier", label="Dubious-gene IQR multiplier"),
    ),
)


class DisambiguationTab(ModuleTabWidget):
    def __init__(self, config: DisambiguationConfig, parent=None):
        super().__init__(SPEC, config, parent)
