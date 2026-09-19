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

--ct_postproc has no separate on/off checkbox: it always runs whenever this tab's
own enable toggle is on (gui/generation/context.py's ct_postproc_enabled is derived
straight from disambig.enabled, not a sibling field) — Post-processing was never
meaningfully independent of Disambiguation in practice, so the extra checkbox only
added a state (Disambiguation on, Post-processing off) nobody used on purpose.

ASR mode: disambiguation_main.py's own argparse only accepts "compute" or
"precomputed" (subworkflows/CT_DISAMBIGUATION/local/disambiguation_main.py) — an
earlier version of this tab only offered "precomputed" as a choice, silently ruling
out live ASR computation.
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import DisambiguationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Disambiguation",
    blurb=(
        "Resolves CAAS convergence direction via ancestral state reconstruction (ASR), "
        "and manages CAAS cluster & gene-level Post-processing.\n\n"
        "ℹ️ Post-Processing Execution Modes:\n"
        "• Exploratory Sweep (run_postproc_exploratory): Performs parameter grid search across "
        "minlen_values × maxcaas_values. Writes output to ${OUTDIR}_exploratory/ and automatically skips downstream tasks.\n"
        "• Production Filtering (run_postproc_filter): Executes single filtering pass with "
        "filter_minlen & filter_maxcaas. Writes output to ${OUTDIR}_final/ and runs downstream modules.\n"
        "• Both Selected: Produces separate script pairs (sbatch_exploratory.sh and sbatch_filtering.sh) for parallel or sequential execution.\n\n"
        "ℹ️ Gene Filter Mode (Advanced):\n"
        "• none: no gene-level filtering.\n"
        "• extreme: drops genes whose CAAS count exceeds the Extreme-gene quantile threshold.\n"
        "• dubious: drops genes flagged by cluster-based IQR screening (Dubious-gene IQR multiplier; needs a cluster file).\n"
        "• both: applies extreme and dubious together."
    ),
    disclaimer=(
        "Accumulation and Scoring need this module's output. Check 'Use precomputed "
        "Disambiguation output' on the Precomputed Run tab instead."
    ),
    essential_fields=(
        Section("Disambiguation parameters (conf/ct_disambiguation.config)"),
        FieldSpec(
            name="ct_disambig_asr_mode",
            label="ASR mode",
            kind="choice",
            choices=("precomputed", "compute"),
        ),
        FieldSpec(name="ct_disambig_asr_model", label="ASR substitution model"),
        FieldSpec(
            name="ct_disambig_asr_cache_dir",
            label="ASR cache directory",
            kind="path_dir",
            help=(
                "Required in both ASR modes — this doubles as the ASR result cache, "
                "not just a place to point at precomputed states. Defaults to "
                "caches/.asr_cache under the run directory (created automatically); "
                "override only to reuse a cache from a previous run or share one "
                "across runs."
            ),
        ),
        Section("Post-processing mode (conf/ct_postproc.config)"),
        FieldSpec(
            name="run_postproc_exploratory",
            label="Run Exploratory Post-Processing Sweep",
            kind="bool",
            help=(
                "Grid-searches Cluster min length x Cluster max CAAS value across the "
                "sweep ranges below (Advanced) instead of a single pass. Writes to "
                "${OUTDIR}_exploratory/ and skips downstream modules — use this to "
                "find good filter thresholds before committing to a production run."
            ),
        ),
        FieldSpec(
            name="run_postproc_filter",
            label="Run Filtering Production Post-Processing",
            kind="bool",
            help=(
                "Single filtering pass using the fixed Cluster min length / Cluster "
                "max CAAS value thresholds below (Advanced). Writes to ${OUTDIR}_final/ "
                "and runs downstream modules — this is the normal production path."
            ),
        ),
    ),
    advanced_fields=(
        FieldSpec(
            name="ct_disambig_posterior_threshold",
            label="Posterior probability threshold",
            help=(
                "Minimum ASR posterior probability for an ancestral state call to be "
                "trusted; calls below this are treated as ambiguous. Statistical "
                "parameter — changing it affects which positions get a directional "
                "call, not just performance."
            ),
        ),
        Section("Post-processing filter thresholds (conf/ct_postproc.config)"),
        FieldSpec(
            name="filter_minlen",
            label="Cluster min length (filter mode)",
            help="Minimum CAAS cluster length to keep. Used by Production Filtering; "
                 "ignored when Exploratory Sweep is selected instead.",
        ),
        FieldSpec(
            name="filter_maxcaas",
            label="Cluster max CAAS value (filter mode)",
            help="Maximum per-cluster CAAS fraction to keep (0-1). Used by Production "
                 "Filtering; ignored when Exploratory Sweep is selected instead.",
        ),
        FieldSpec(
            name="gene_filter_mode",
            label="Gene filter mode",
            kind="choice",
            choices=("none", "extreme", "dubious", "both"),
            help=(
                "'extreme' drops genes whose CAAS count is a statistical outlier "
                "(Extreme-gene quantile threshold, below). 'dubious' drops genes "
                "flagged by cluster-based IQR screening (Dubious-gene IQR multiplier, "
                "below; needs a cluster file). 'both' applies both filters; 'none' "
                "disables gene-level filtering entirely."
            ),
        ),
        Section("Sensitivity and robustness testing"),
        FieldSpec(name="asr_robustness", label="Run ASR Robustness diagnostics report", kind="bool"),
        Section("Performance and batching"),
        FieldSpec(name="ct_disambig_max_tasks_per_child", label="Max tasks per worker child"),
        FieldSpec(name="ct_disambig_batch_size", label="Disambiguation genes per batch"),
        Section("Exploratory parameter sweep values (conf/ct_postproc.config)"),
        FieldSpec(
            name="minlen_values",
            label="Cluster min length sweep (exploratory mode)",
            help="Comma-separated Cluster min length values to grid-search. Used only "
                 "when Exploratory Sweep is selected.",
        ),
        FieldSpec(
            name="maxcaas_values",
            label="Cluster max CAAS sweep (exploratory mode)",
            help="Comma-separated Cluster max CAAS value values to grid-search. Used "
                 "only when Exploratory Sweep is selected.",
        ),
        Section("Gene-level outlier thresholds (conf/ct_postproc.config)"),
        FieldSpec(
            name="extreme_threshold",
            label="Extreme-gene quantile threshold",
            help="Quantile above which a gene's CAAS count is considered an outlier. "
                 "Only used when Gene filter mode is 'extreme' or 'both'.",
        ),
        FieldSpec(
            name="iqr_multiplier",
            label="Dubious-gene IQR multiplier",
            help="IQR multiplier for the cluster-based dubious-gene screen. Only used "
                 "when Gene filter mode is 'dubious' or 'both' (needs a cluster file).",
        ),
    ),
)


class DisambiguationTab(ModuleTabWidget):
    def __init__(self, config: DisambiguationConfig, parent=None):
        super().__init__(SPEC, config, parent)

    def retranslate(self, lang: str = "en") -> None:
        super().retranslate(lang)
        from gui.i18n import tr
        main_blurb = tr(
            "Resolves CAAS convergence direction via ancestral state reconstruction (ASR), "
            "and manages CAAS cluster & gene-level Post-processing.",
            lang
        )
        guidance = tr("Postproc Guidance", lang)
        self.blurb_label.setText(f"{main_blurb}\n\n{guidance}")
