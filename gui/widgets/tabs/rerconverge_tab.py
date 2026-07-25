#!/usr/bin/env python3
# rerconverge_tab.py — RERconverge module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Off by default (RUN_RER=false in the reference scripts). The 4 rer_tool sub-steps
are individual checkboxes jointly building --rer_tool, same pattern as CAAS's
ct_tool discovery/resample/bootstrap checkboxes.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import QCheckBox

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import RerConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="RERconverge",
    blurb=(
        "Correlates relative evolutionary rates (RER) across gene trees with the "
        "phenotype, using permutation batches for significance."
    ),
    disclaimer="Scoring needs a scoring_rer_input fallback per phenotype row (Runtime tab) when this is off.",
    essential_fields=(
        FieldSpec(name="gene_trees", label="Gene trees file", kind="path_file"),
        FieldSpec(name="rer_perm_batches", label="Permutation batches"),
        FieldSpec(name="rer_perms_per_batch", label="Permutations per batch"),
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

        self._essential_form.insertRow(0, "RER tools (--rer_tool)", self.rer_tool_build_trait)
        self._essential_form.insertRow(1, "", self.rer_tool_build_tree)
        self._essential_form.insertRow(2, "", self.rer_tool_build_matrix)
        self._essential_form.insertRow(3, "", self.rer_tool_continuous)

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
