#!/usr/bin/env python3
# disambiguation_tab.py — Disambiguation module tab (bundles Post-processing sub-section).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
DisambiguationTab: field specification of the Disambiguation module, which also
holds the Post-processing parameters (conf/ct_disambiguation.config and
conf/ct_postproc.config).

Post-processing has no tab or enable checkbox of its own: gui/generation/context.py
derives ct_postproc_enabled from this tab's enable toggle (and from the Precomputed
Run tab not supplying a post-processing result).

ASR: a gene's reconstruction is read from the cache directory when it is there and
computed with PAML (and written to the cache) when it is not, so there is no mode to
choose.

Imported by: gui/widgets/main_window.py
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
            name="ct_disambig_asr_model",
            label="ASR substitution model",
            kind="choice",
            # Empirical amino-acid matrices run through PAML codeml; the full set of
            # models is MODEL_SPECS in
            # subworkflows/CT_DISAMBIGUATION/local/src/asr/reconstruct.py.
            choices=("lg", "wag", "jtt", "dayhoff"),
            importance="default",
        ),
        FieldSpec(
            name="ct_disambig_asr_cache_dir",
            label="ASR cache directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "This is the ASR result cache: a gene found here is read, a gene "
                "missing from it is computed and written to it. Leaving it blank auto-"
                "generates a working default (repo_dir/caches/.asr_cache, created "
                "automatically; see gui/generation/templates/run_single.sh.j2's bash "
                "fallback). Override only to reuse a cache from a previous run or "
                "share one across runs — cache location has no effect on statistical "
                "validity."
            ),
        ),
        FieldSpec(
            name="ct_disambig_hypotheses_pairs",
            label="Contrast hypotheses pairs file for the observed scoring (optional)",
            kind="path_file",
            importance="optional",
        ),
        Section("Post-processing mode (conf/ct_postproc.config)"),
        FieldSpec(
            name="run_postproc_exploratory",
            label="Run Exploratory Post-Processing Sweep",
            kind="bool",
            importance="optional",
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
            importance="optional",
            help=(
                "Single filtering pass using the fixed Cluster min length / Cluster "
                "max CAAS value thresholds below (Advanced). Writes to ${OUTDIR}_final/ "
                "and runs downstream modules — this is the normal production path."
            ),
        ),
        FieldSpec(
            name="caas_map_dir",
            label="Per-gene MAP directory (optional)",
            kind="path_dir",
            help="Directory of the trimmer's per-gene MAP tables. When set, cluster "
                 "trains measure their span in untrimmed alignment columns (columns the "
                 "trimmer removed cannot hold a CAAS), in the observed filter and in the "
                 "permulation null alike; a gene without a MAP file keeps trimmed "
                 "coordinates. VEP maps positions with the same tables and requires it. "
                 "Cannot be generated in-house: see "
                 "github.com/nozerorma/ortholog_characterizator.",
            importance="optional",
        ),
    ),
    advanced_fields=(
        FieldSpec(
            name="ct_disambig_posterior_threshold",
            label="Posterior probability threshold",
            importance="default",
            help=(
                "Residues with an ASR posterior below this value are dropped from the "
                "recorded distribution of a node (the most probable residue is always "
                "kept). The dropped mass is charged as worst case in the position "
                "score, so a higher value can only lower scores. Statistical "
                "parameter — changing it changes the scores, not just performance."
            ),
        ),
        Section("Post-processing filter thresholds (conf/ct_postproc.config)"),
        FieldSpec(
            name="filter_minlen",
            label="Cluster min length (filter mode)",
            importance="default",
            help="Minimum CAAS cluster length to keep. Used by Production Filtering; "
                 "ignored when Exploratory Sweep is selected instead.",
        ),
        FieldSpec(
            name="filter_maxcaas",
            label="Cluster max CAAS value (filter mode)",
            importance="default",
            help="Maximum per-cluster CAAS fraction to keep (0-1). Used by Production "
                 "Filtering; ignored when Exploratory Sweep is selected instead.",
        ),
        FieldSpec(
            name="gene_filter_mode",
            label="Gene filter mode",
            kind="choice",
            choices=("none", "extreme", "dubious", "both"),
            importance="default",
            help=(
                "'extreme' drops genes whose CAAS count is a statistical outlier "
                "(Extreme-gene quantile threshold, below). 'dubious' drops genes "
                "flagged by cluster-based IQR screening (Dubious-gene IQR multiplier, "
                "below; needs a cluster file). 'both' applies both filters; 'none' "
                "disables gene-level filtering entirely."
            ),
        ),
        FieldSpec(
            name="remove_caas_clusters",
            label="Remove CAAS spatial clusters",
            kind="bool",
            importance="default",
            help=(
                "Discard individual positions flagged as spatial CAAS clusters from the "
                "discovery dataset. Decoupled from gene-level filtering: when enabled, "
                "spatial clusters are removed even if Gene filter mode is set to 'none'."
            ),
        ),
        Section("ASR diagnostics"),
        FieldSpec(
            name="asr_diagnostics",
            label="Run ASR diagnostics report",
            kind="bool",
            importance="optional",
        ),
        Section("Performance and batching"),
        FieldSpec(name="ct_disambig_max_tasks_per_child", label="Max tasks per worker child", importance="optional"),
        Section("Exploratory parameter sweep values (conf/ct_postproc.config)"),
        # Rated "optional" although borderline: these values shape only the
        # diagnostic sweep grid of Exploratory mode, not the production filter
        # thresholds above, so a poor sweep range wastes a diagnostic run without
        # affecting a final result.
        FieldSpec(
            name="minlen_values",
            label="Cluster min length sweep (exploratory mode)",
            importance="optional",
            help="Comma-separated Cluster min length values to grid-search. Used only "
                 "when Exploratory Sweep is selected.",
        ),
        FieldSpec(
            name="maxcaas_values",
            label="Cluster max CAAS sweep (exploratory mode)",
            importance="optional",
            help="Comma-separated Cluster max CAAS value values to grid-search. Used "
                 "only when Exploratory Sweep is selected.",
        ),
        Section("Gene-level outlier thresholds (conf/ct_postproc.config)"),
        FieldSpec(
            name="extreme_threshold",
            label="Extreme-gene quantile threshold",
            importance="default",
            help="Quantile above which a gene's CAAS count is considered an outlier. "
                 "Only used when Gene filter mode is 'extreme' or 'both'.",
        ),
        FieldSpec(
            name="iqr_multiplier",
            label="Dubious-gene IQR multiplier",
            importance="default",
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
