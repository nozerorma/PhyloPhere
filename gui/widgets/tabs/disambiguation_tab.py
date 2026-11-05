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
        "ℹ️ Scripts to generate:\n"
        "• Exploratory (diagnostic): works on the real labeling only and never builds the permulation null. It sweeps "
        "minlen_values × maxcaas_values, reports what the gene filter would remove, skips the downstream modules and writes to "
        "<results dir>/<phenotype>_exploratory/.\n"
        "• Production (complete analysis): one script that runs CAAS, Disambiguation, post-processing with the parameters set here, "
        "the permulation null built with them and the enabled downstream modules; it writes to <results dir>/<phenotype>_complete/ and does not "
        "reuse the exploratory output.\n"
        "• Strict analysis (hypotheses without conserved pairs): the exploratory script also runs a strict exploratory pass with "
        "its own report, after the normal one and in the same task, and the production script also runs the strict pass after "
        "the complete one.\n\n"
        "ℹ️ Gene filter (Advanced):\n"
        "• Extreme: a gene whose CAAS count exceeds the far fence of its gene-length stratum (Q3 + k_far x scale), any gene. "
        "Uses gene_filter_strata and gene_filter_k_far.\n"
        "• Dubious: a gene whose count exceeds the near fence of its stratum (Q3 + k_near x scale) AND that carries a cluster train. "
        "Uses gene_filter_strata, gene_filter_k_near and the cluster parameters (filter_minlen, filter_maxcaas, filter_min_positions).\n"
        "• Workflow: read the exploratory report (it shows both marks whatever the switches say), then set the cluster values and "
        "the switches here and run the production script: the null is built there with the same values as the observed table. "
        "Cluster parameters are shared by the normal and strict analyses, because they act in the replay; the strict gene filter "
        "can have its own settings (Advanced)."
    ),
    disclaimer=(
        "Accumulation and Scoring need this module's output. Check 'Use precomputed "
        "Disambiguation output' on the Precomputed Run tab instead."
    ),
    essential_fields=(
        Section("Disambiguation: ancestral reconstruction (conf/ct_disambiguation.config)"),
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
        Section("Scripts to generate"),
        FieldSpec(
            name="run_postproc_exploratory",
            label="Generate exploratory script (diagnostic)",
            kind="bool",
            importance="optional",
            help=(
                "Real labeling only, never the permulation null: sweeps cluster length and CAAS density (Advanced) and reports what "
                "the gene filter would remove. Skips the downstream modules; writes to <results dir>/<phenotype>_exploratory/."
            ),
        ),
        FieldSpec(
            name="run_postproc_filter",
            label="Generate production script (complete analysis)",
            kind="bool",
            importance="optional",
            help=(
                "The complete analysis in one script: CAAS, Disambiguation, post-processing with the parameters below, the "
                "permulation null built with them and the enabled downstream modules. Writes to <results dir>/<phenotype>_complete/."
            ),
        ),
        FieldSpec(
            name="run_strict_analysis",
            label="Run strict analysis",
            kind="bool",
            importance="optional",
            help=(
                "Over the hypotheses without conserved pairs: the exploratory script also runs a strict exploratory pass (its own "
                "report) after the normal one, and the production script also runs the strict pass after the complete one. "
                "Writes the strict tables and null of the normal run. Needs a min divergent fraction below 1. Unchecked: every "
                "strict option below is disabled."
            ),
        ),
        Section("Cluster filter (conf/ct_postproc.config)"),
        FieldSpec(
            name="filter_minlen",
            label="Cluster min length",
            importance="default",
            help="Minimum interval span, in alignment columns, of a CAAS cluster. Used by the production script and as the selected "
                 "pair of the exploratory sweep. It acts in the replay, so the normal and the strict analyses share it.",
        ),
        FieldSpec(
            name="filter_min_positions",
            label="Minimum candidates per interval",
            importance="optional",
            help="Positive integer. Blank keeps the legacy gene-count gate of the minimum cluster length.",
        ),
        FieldSpec(
            name="filter_maxcaas",
            label="Cluster max CAAS value",
            importance="default",
            help="Maximum CAAS density of a cluster window (0-1). Used by the production script and as the selected pair of the "
                 "exploratory sweep; shared by the normal and the strict analyses.",
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
        FieldSpec(
            name="caas_map_dir",
            label="Per-gene MAP directory (optional)",
            kind="path_dir",
            help="Per-gene MAP tables from the Ortholog Characterizator; cluster spans and gene lengths use untrimmed columns, and VEP needs them.",
            importance="optional",
        ),
        Section("Gene filter: Extreme and Dubious genes (conf/ct_postproc.config)"),
        FieldSpec(
            name="remove_extreme_genes",
            label="Remove Extreme genes",
            kind="bool",
            importance="default",
            help=(
                "Remove genes whose CAAS count exceeds the far fence of their gene-length stratum (Q3 + far multiplier x scale), "
                "whether or not they carry a cluster train. Uses Gene-length strata and Far-fence IQR multiplier. The exploratory "
                "report shows what this rule would remove; the production script applies it to the observed table and to every "
                "cycle of the permulation null."
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
    ),
    advanced_fields=(
        Section("Disambiguation: fine tuning (conf/ct_disambiguation.config)"),
        FieldSpec(
            name="ct_disambig_posterior_threshold",
            label="Posterior probability threshold",
            importance="default",
            help=(
                "Residues with an ASR posterior below this are dropped from a node; the dropped mass counts as worst case, so higher values only lower scores."
            ),
        ),
        FieldSpec(
            name="asr_diagnostics",
            label="Run ASR diagnostics report",
            kind="bool",
            importance="optional",
            help="Renders the ASR diagnostics report (9.ASR_diagnostics.html); default on.",
        ),
        FieldSpec(
            name="ct_disambig_max_tasks_per_child",
            label="Max tasks per worker child",
            importance="optional",
            help="Number of tasks a worker process runs before it is recycled; default 50.",
        ),
        Section("Cluster filter: exploratory sweep values (conf/ct_postproc.config)"),
        # Rated "optional" although borderline: these values shape only the
        # diagnostic sweep grid of the exploratory script, not the production filter
        # thresholds, so a poor sweep range wastes a diagnostic run without
        # affecting a final result.
        FieldSpec(
            name="minlen_values",
            label="Cluster min length sweep (exploratory mode)",
            importance="optional",
            help="Comma-separated Cluster min length values to grid-search. Used only by the exploratory script.",
        ),
        FieldSpec(
            name="maxcaas_values",
            label="Cluster max CAAS sweep (exploratory mode)",
            importance="optional",
            help="Comma-separated Cluster max CAAS value values to grid-search. Used only by the exploratory script.",
        ),
        Section("Gene filter: length strata and fences (conf/ct_postproc.config)"),
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
        Section("Strict analysis: gene filter (conf/ct_postproc.config)"),
        FieldSpec(
            name="strict_same_gene_filter",
            label="Strict: same gene filter as the normal run",
            kind="bool",
            importance="optional",
            enabled_when="run_strict_analysis",
            help=(
                "Checked: the strict analysis removes genes with the settings above, and the fields below are disabled. Unchecked: "
                "it uses the fields below. Disabled while 'Run strict analysis' is unchecked. The strict gene marks come from the strict counts, so they can differ from the normal "
                "ones. The strict null takes these values when the production script builds it. Cluster parameters have no strict "
                "counterpart: they act in the replay, which produces the normal and the strict tables together."
            ),
        ),
        FieldSpec(
            name="strict_remove_extreme_genes",
            label="Strict: remove Extreme genes",
            kind="bool",
            importance="optional",
            enabled_when="run_strict_analysis",
            disabled_when="strict_same_gene_filter",
            help="Remove genes above the far fence of their length stratum in the strict analysis (strict counts).",
        ),
        FieldSpec(
            name="strict_remove_dubious_genes",
            label="Strict: remove Dubious genes",
            kind="bool",
            importance="optional",
            enabled_when="run_strict_analysis",
            disabled_when="strict_same_gene_filter",
            help="Remove genes above the near fence of their stratum that carry a cluster train in the strict analysis (strict counts).",
        ),
        FieldSpec(
            name="strict_gene_filter_strata",
            label="Strict: gene-length strata",
            importance="optional",
            enabled_when="run_strict_analysis",
            disabled_when="strict_same_gene_filter",
            help="Number of equal-size gene-length strata of the strict analysis.",
        ),
        FieldSpec(
            name="strict_gene_filter_k_far",
            label="Strict: far-fence IQR multiplier",
            importance="optional",
            enabled_when="run_strict_analysis",
            disabled_when="strict_same_gene_filter",
            help="Far-fence multiplier of the strict analysis (Extreme above Q3 + k x scale).",
        ),
        FieldSpec(
            name="strict_gene_filter_k_near",
            label="Strict: near-fence IQR multiplier",
            importance="optional",
            enabled_when="run_strict_analysis",
            disabled_when="strict_same_gene_filter",
            help="Near-fence multiplier of the strict analysis (Dubious above it with a cluster train).",
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
