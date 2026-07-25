#!/usr/bin/env python3
# caas_tab.py — CAAS / CT module tab (contrast selection: discovery, resample, bootstrap).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
The discovery/resample/bootstrap sub-steps are separate boolean checkboxes on the
config (ct_tool_discovery/resample/bootstrap, see gui/models/modules.py) rather than
FieldSpec entries, since ModuleTabWidget's field kinds don't cover "3 checkboxes
that jointly build one comma-separated flag" — they're added directly in __init__.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import QCheckBox

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import CaasConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="CAAS / Contrast Selection",
    blurb=(
        "Runs CAAStools' discovery, resample, and bootstrap steps to find convergent "
        "amino-acid substitutions (CAAS) associated with the phenotype."
    ),
    disclaimer=(
        "Disambiguation and Accumulation need this module's output. Supply "
        "discovery_from/resample_from/bootstrap_from below to feed them precomputed results."
    ),
    essential_fields=(
        FieldSpec(name="caas_config_path", label="CAAS config file", kind="path_file"),
        FieldSpec(name="perms_cycles", label="Permutation cycles"),
        FieldSpec(name="caas_full_perms", label="Full-pool permulations"),
        FieldSpec(
            name="caas_permulation_enrichment", label="Permulation-excess for enrichment", kind="bool"
        ),
    ),
    fallback_fields=(
        FieldSpec(name="discovery_from", label="Discovery output (discovery_from)", kind="path_dir"),
        FieldSpec(name="resample_from", label="Resample output (resample_from)", kind="path_dir"),
        FieldSpec(name="bootstrap_from", label="Bootstrap output (bootstrap_from)", kind="path_dir"),
    ),
)


class CaasTab(ModuleTabWidget):
    def __init__(self, config: CaasConfig, parent=None):
        super().__init__(SPEC, config, parent)

        # ct_tool discovery/resample/bootstrap: 3 checkboxes jointly building --ct_tool.
        self.ct_tool_discovery = QCheckBox("discovery")
        self.ct_tool_discovery.setChecked(config.ct_tool_discovery)
        self.ct_tool_discovery.toggled.connect(self._on_ct_tool_discovery)

        self.ct_tool_resample = QCheckBox("resample")
        self.ct_tool_resample.setChecked(config.ct_tool_resample)
        self.ct_tool_resample.toggled.connect(self._on_ct_tool_resample)

        self.ct_tool_bootstrap = QCheckBox("bootstrap")
        self.ct_tool_bootstrap.setChecked(config.ct_tool_bootstrap)
        self.ct_tool_bootstrap.toggled.connect(self._on_ct_tool_bootstrap)

        self._essential_form.insertRow(0, "CT tools (--ct_tool)", self.ct_tool_discovery)
        self._essential_form.insertRow(1, "", self.ct_tool_resample)
        self._essential_form.insertRow(2, "", self.ct_tool_bootstrap)

    def _on_ct_tool_discovery(self, value: bool) -> None:
        self._config.ct_tool_discovery = value
        self.changed.emit()

    def _on_ct_tool_resample(self, value: bool) -> None:
        self._config.ct_tool_resample = value
        self.changed.emit()

    def _on_ct_tool_bootstrap(self, value: bool) -> None:
        self._config.ct_tool_bootstrap = value
        self.changed.emit()
