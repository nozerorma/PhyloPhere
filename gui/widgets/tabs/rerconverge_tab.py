#!/usr/bin/env python3
# rerconverge_tab.py — RERconverge module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Off by default (RUN_RER=false in the reference scripts). The 4 rer_tool sub-steps
are individual checkboxes jointly building --rer_tool, same pattern as CAAS's
ct_tool discovery/resample/bootstrap checkboxes.

No "gene_set" mode here: RERconverge operates genome-wide across all
gene trees. Subworkflows take full gene trees as standard input.
main.nf always passes Channel.empty() for both gene-set channels, and RER_MAIN's
own take: comment says so explicitly ("Kept for backwards-compatible signature;
currently unused."). Exposing a mode picker here would silently do nothing, so
this is left as a real pipeline gap rather than a GUI gap — see the spawned
follow-up task for wiring it, not a fix to make here.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import QCheckBox, QLabel

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import RerConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="RERconverge",
    blurb=(
        "Correlates relative evolutionary rates (RER) across gene trees with the "
        "phenotype, using permutation batches for significance."
    ),
    disclaimer=(
        "Scoring needs this module's output when it's off. Check 'Use precomputed "
        "RERconverge output' on the Precomputed Run tab — it's auto-derived per "
        "phenotype from one base path, no per-row entry needed."
    ),
    essential_fields=(
        Section("RER trait analysis"),
        # Routes which statistical model applies to the trait (continuous vs
        # binary RERconverge machinery) — model-selection parameter.
        FieldSpec(
            name="rer_trait_mode",
            label="Trait type routing",
            kind="choice",
            choices=("auto", "continuous", "binary"),
            importance="default",
        ),
        # Transform choice directly affects the RER-trait correlation's
        # statistical validity (comparable to a model-selection parameter).
        FieldSpec(
            name="rer_transform",
            label="Trait transform",
            kind="choice",
            choices=("ha_logit", "auto", "logit", "arcsin", "log10", "none"),
            importance="default",
        ),
        # p.perm (raw permutation p) vs p.adj (BH-adjusted) changes what
        # "significant" means downstream — an FDR/p-value-threshold-like choice.
        FieldSpec(
            name="rer_pval_column",
            label="Report p-value column",
            kind="choice",
            choices=("p.perm", "p.adj"),
            importance="default",
        ),
        # validate.py: require(rer.gene_trees, ...) whenever rer.enabled; also
        # user-supplied with no in-house generation path (promoted from advanced
        # to essential per commit history) — a strong "required" field.
        FieldSpec(
            name="gene_trees",
            label="Gene trees file",
            kind="path_file",
            help=(
                "Per-gene phylogenetic trees (one per alignment), required whenever "
                "RERconverge runs. No in-house generation path: building per-gene "
                "trees is a separate phylogenetics task outside this pipeline's scope "
                "— supply one built with your own tree-inference tool of choice."
            ),
            importance="required",
        ),
    ),
    advanced_fields=(
        Section("Gene tree and matrix output paths"),
        FieldSpec(name="trait_out", label="Trait output path (build_trait)", importance="optional"),
        FieldSpec(name="trees_out", label="Trees output path (build_tree)", importance="optional"),
        FieldSpec(name="matrix_out", label="Matrix output path (build_matrix)", importance="optional"),
        Section("Quality thresholds and winsorization"),
        FieldSpec(name="rer_minsp", label="Minimum species per gene", importance="default"),
        FieldSpec(name="winsorize_rer", label="Winsorize RER threshold", importance="default"),
        FieldSpec(name="winsorize_trait", label="Winsorize trait threshold", importance="default"),
        Section("Permutation null calibration"),
        FieldSpec(name="rer_perm_batches", label="Permutation batches", importance="default"),
        FieldSpec(name="rer_perms_per_batch", label="Permutations per batch", importance="default"),
        Section("Binary mode options"),
        FieldSpec(
            name="rer_binary_clade",
            label="Binary foreground clade",
            kind="choice",
            choices=("all", "ancestral", "terminal"),
            importance="default",
        ),
        FieldSpec(name="rer_min_pos", label="Min independent foreground lineages", importance="default"),
        Section("Reporting and visualization thresholds"),
        FieldSpec(name="rer_pval_threshold", label="Report p-value threshold", importance="default"),
        FieldSpec(name="rer_top_n_labels", label="Top-N labels on correlation plot", importance="optional"),
        FieldSpec(name="rer_universe_file", label="RER tested-gene universe file", kind="path_file", importance="optional"),
        FieldSpec(name="rer_gene_scores", label="Cross-module gene scores (fcs_stats.tsv)", kind="path_file", importance="optional"),
    ),
)


class RerConvergeTab(ModuleTabWidget):
    def __init__(self, config: RerConfig, parent=None):
        super().__init__(SPEC, config, parent)

        self.rer_tool_build_trait = QCheckBox("build_trait")
        self.rer_tool_build_trait.setChecked(config.rer_tool_build_trait)
        self.rer_tool_build_trait.toggled.connect(self._on_build_trait)

        self.rer_tool_build_tree = QCheckBox("build_tree")
        self.rer_tool_build_tree.setChecked(config.rer_tool_build_tree)
        self.rer_tool_build_tree.toggled.connect(self._on_build_tree)

        self.rer_tool_build_matrix = QCheckBox("build_matrix")
        self.rer_tool_build_matrix.setChecked(config.rer_tool_build_matrix)
        self.rer_tool_build_matrix.toggled.connect(self._on_build_matrix)

        self.rer_tool_continuous = QCheckBox("continuous")
        self.rer_tool_continuous.setChecked(config.rer_tool_continuous)
        self.rer_tool_continuous.toggled.connect(self._on_continuous)

        self._rer_tool_label = QLabel("RER tools (--rer_tool)")
        self._essential_form.insertRow(0, self._rer_tool_label, self.rer_tool_build_trait)
        self._essential_form.insertRow(1, "", self.rer_tool_build_tree)
        self._essential_form.insertRow(2, "", self.rer_tool_build_matrix)
        self._essential_form.insertRow(3, "", self.rer_tool_continuous)

    def retranslate(self, lang: str = "en") -> None:
        super().retranslate(lang)
        from gui.i18n import tr
        self._rer_tool_label.setText(tr("RER tools (--rer_tool)", lang))
        self.rer_tool_build_trait.setText(tr("build_trait", lang))
        self.rer_tool_build_tree.setText(tr("build_tree", lang))
        self.rer_tool_build_matrix.setText(tr("build_matrix", lang))
        self.rer_tool_continuous.setText(tr("continuous", lang))

    def _on_build_trait(self, value: bool) -> None:
        self._config.rer_tool_build_trait = value
        self.changed.emit()

    def _on_build_tree(self, value: bool) -> None:
        self._config.rer_tool_build_tree = value
        self.changed.emit()

    def _on_build_matrix(self, value: bool) -> None:
        self._config.rer_tool_build_matrix = value
        self.changed.emit()

    def _on_continuous(self, value: bool) -> None:
        self._config.rer_tool_continuous = value
        self.changed.emit()
