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
        "filter_minlen & filter_maxcaas. Writes output to ${OUTDIR}_complete/ and runs downstream modules.\n"
        "• Both Selected: Run the exploratory scripts first; the complete scripts reuse their output.\n"
        "• Strict Analysis: Reuses the normal run's strict tables for post-processing and downstream modules.\n\n"
        "ℹ️ Gene filter (Advanced):\n"
        "• Extreme: a gene whose CAAS count exceeds the far fence of its gene-length stratum (Q3 + k_far x scale), any gene. "
        "Uses gene_filter_strata and gene_filter_k_far.\n"
        "• Dubious: a gene whose count exceeds the near fence of its stratum (Q3 + k_near x scale) AND that carries a cluster train. "
        "Uses gene_filter_strata, gene_filter_k_near and the cluster parameters (filter_minlen, filter_maxcaas, filter_min_positions).\n"
        "• Workflow: the exploratory sweep reports what each rule would remove; in the complete run, decide from that report "
        "to remove Extreme genes, Dubious genes, both or none. The same switches apply to the observed table and to the permulation null."
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
            help="Empirical amino-acid substitution matrix for PAML codeml ancestral reconstruction: lg (default), wag, jtt or dayhoff.",
        ),
        FieldSpec(
            name="ct_disambig_asr_cache_dir",
            label="ASR cache directory (optional)",
            kind="path_dir",
            importance="optional",
            help=(
                "ASR result cache: genes found there are read, others computed and stored. Blank uses repo_dir/caches/.asr_cache; reuse a cache to skip reconstruction."
            ),
        ),
        FieldSpec(
            name="ct_disambig_hypotheses_pairs",
            label="Contrast hypotheses pairs file for the observed scoring (optional)",
            kind="path_file",
            importance="optional",
            help="contrast_hypotheses_pairs.tsv weighting the observed scoring; blank uses this run's contrast selection, or auto-discovery in standalone runs.",
        ),
        Section("Post-processing mode (conf/ct_postproc.config)"),
        FieldSpec(
            name="run_postproc_exploratory",
            label="Run Exploratory Post-Processing Sweep",
            kind="bool",
            importance="optional",
            help=(
                "Sweeps cluster length and CAAS density thresholds (Advanced) to choose filters; writes to ${OUTDIR}_exploratory/ and skips downstream modules."
            ),
        ),
        FieldSpec(
            name="run_postproc_filter",
            label="Run Filtering Production Post-Processing",
            kind="bool",
            importance="optional",
            help=(
                "Single filtering pass with the fixed cluster thresholds below; writes to ${OUTDIR}_complete/ and runs downstream modules (normal production path)."
            ),
        ),
        FieldSpec(
            name="run_strict_analysis",
            label="Generate the strict analysis scripts",
            kind="bool",
            importance="optional",
            help=(
                "Also generates scripts that repeat the run from post-processing with only hypotheses without conserved pairs; run the normal scripts first. Needs a min divergent fraction below 1."
            ),
        ),
        FieldSpec(
            name="caas_map_dir",
            label="Per-gene MAP directory (optional)",
            kind="path_dir",
            help="Per-gene MAP tables from the Ortholog Characterizator; cluster spans use untrimmed columns, and VEP needs them.",
            importance="optional",
        ),
    ),
    advanced_fields=(
        FieldSpec(
            name="ct_disambig_posterior_threshold",
            label="Posterior probability threshold",
            importance="default",
            help=(
                "Residues with an ASR posterior below this are dropped from a node; the dropped mass counts as worst case, so higher values only lower scores."
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
            name="filter_min_positions",
            label="Minimum candidates per interval (filter mode)",
            importance="optional",
            help="Positive integer. Blank keeps the legacy gene-count gate of the minimum cluster length.",
        ),
        FieldSpec(
            name="filter_maxcaas",
            label="Cluster max CAAS value (filter mode)",
            importance="default",
            help="Maximum per-cluster CAAS fraction to keep (0-1). Used by Production "
                 "Filtering; ignored when Exploratory Sweep is selected instead.",
        ),
        FieldSpec(
            name="remove_extreme_genes",
            label="Remove Extreme genes",
            kind="bool",
            importance="default",
            help=(
                "Remove genes whose CAAS count exceeds the far fence of their gene-length stratum (Q3 + far multiplier x scale), "
                "whether or not they carry a cluster train. Uses Gene-length strata and Far-fence IQR multiplier. "
                "Applied to the observed table and to every cycle of the permulation null. The exploratory run reports "
                "what this rule would remove; decide here for the complete run."
            ),
        ),
        FieldSpec(
            name="remove_dubious_genes",
            label="Remove Dubious genes",
            kind="bool",
            importance="default",
            help=(
                "Remove genes whose CAAS count exceeds the near fence of their stratum (Q3 + near multiplier x scale) AND that "
                "carry a cluster train. Uses Gene-length strata, Near-fence IQR multiplier and the cluster parameters "
                "(minimum length, maximum CAAS value, minimum candidates per interval). Applied to the observed table "
                "and to every cycle of the permulation null."
            ),
        ),
        FieldSpec(
            name="remove_caas_clusters",
            label="Remove CAAS spatial clusters",
            kind="bool",
            importance="default",
            help=(
                "Discard positions in spatial CAAS clusters from the discovery, independent of the gene removal switches."
            ),
        ),
        Section("ASR diagnostics"),
        FieldSpec(
            name="asr_diagnostics",
            label="Run ASR diagnostics report",
            kind="bool",
            importance="optional",
            help="Renders the ASR diagnostics report (9.ASR_diagnostics.html); default on.",
        ),
        Section("Performance and batching"),
        FieldSpec(
            name="ct_disambig_max_tasks_per_child",
            label="Max tasks per worker child",
            importance="optional",
            help="Number of tasks a worker process runs before it is recycled; default 50.",
        ),
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
        Section("Gene-level outlier fences (conf/ct_postproc.config)"),
        FieldSpec(
            name="gene_filter_strata",
            label="Gene-length strata",
            importance="default",
            help="Number of equal-size gene-length strata in which the outlier fences are computed, so that a gene "
                 "is compared with genes of its own length. Used by both Remove Extreme genes and Remove Dubious genes.",
        ),
        FieldSpec(
            name="gene_filter_k_far",
            label="Far-fence IQR multiplier",
            importance="default",
            help="A gene whose CAAS count exceeds Q3 + k x scale of its length stratum is Extreme (scale: the larger of the "
                 "stratum IQR and the IQR of all genes). Only used by Remove Extreme genes.",
        ),
        FieldSpec(
            name="gene_filter_k_near",
            label="Near-fence IQR multiplier",
            importance="default",
            help="A gene above this lower fence is Dubious only when it also carries a cluster train "
                 "(needs a cluster file). Only used by Remove Dubious genes.",
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
        gene_guidance = tr("Gene filter guidance", lang)
        self.blurb_label.setText(f"{main_blurb}\n\n{guidance}\n\n{gene_guidance}")
