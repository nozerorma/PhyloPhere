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

--contrast_selection itself has no separate on/off checkbox: the reference scripts
always bundle it with CAAS (run_phenotype_single_primates.sh's RUN_CAAS block emits
it unconditionally), so gui/generation/templates emit it the same way rather than
exposing a redundant toggle here.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import QCheckBox, QLabel

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import CaasConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="CAAS / Contrast Selection",
    blurb=(
        "Runs CAAStools' discovery, resample, and bootstrap steps to find convergent "
        "amino-acid substitutions (CAAS) associated with the phenotype."
    ),
    disclaimer=(
        "Disambiguation and Accumulation need this module's output. Check Discovery/"
        "Resample/Bootstrap on the Precomputed Run tab to feed them precomputed "
        "results instead."
    ),
    essential_fields=(
        Section("Discovery params"),
        FieldSpec(name="patterns", label="CT patterns"),
        FieldSpec(name="min_divergent_fraction", label="Min divergent fraction"),
        FieldSpec(name="caap_mode", label="CAAP mode (properties-based)", kind="bool"),
        FieldSpec(name="multi_hypothesis", label="Multi-hypothesis mode", kind="bool"),
        FieldSpec(
            name="max_fop",
            label="Max FOP hypotheses (H1..Hn)",
            placeholder="alternative Dunn-independent hypotheses per contrast (observed + null); default 100",
        ),
        Section("Resample / Bootstrap params"),
        FieldSpec(name="chunk_size", label="Resampled groups per output file"),
        FieldSpec(name="include_b0", label="Include main hypothesis (b0)", kind="bool"),
        FieldSpec(name="resample_use_n", label="Use sample size counts (n/c)", kind="bool"),
        FieldSpec(name="perm_strategy", label="Permutation strategy", kind="choice", choices=("auto", "OU", "BM")),
        FieldSpec(
            name="perm_pool_size",
            label="Permulation pool size",
            placeholder="accepted permulations to harvest (position-level null)",
        ),
        FieldSpec(
            name="caas_full_perms",
            label="Full-pool permulations",
            placeholder="drawn from the pool for the CAAS FCS null",
        ),
        FieldSpec(
            name="max_tries",
            label="Max permulation tries",
            placeholder="draw budget; raised 50% up to twice if the pool falls short",
        ),
        FieldSpec(
            name="caas_perms_fop",
            label="Mirror FOP harvest in the CAAS permulation null",
            kind="bool",
        ),
    ),
    advanced_fields=(
        Section("Missingness parameters in discovery/bootstrap modes and otherwise"),
        FieldSpec(
            name="caas_config_path",
            label="CAAS config file (auto-derived if empty)",
            kind="path_file",
            placeholder="auto-derived from contrast_selection.nf when empty",
        ),
        FieldSpec(name="max_bg_gaps_fraction", label="Max background gaps fraction"),
        FieldSpec(name="max_fg_gaps_fraction", label="Max foreground gaps fraction"),
        FieldSpec(name="max_gaps_fraction", label="Max any-gaps fraction"),
        FieldSpec(name="max_bg_miss_fraction", label="Max background missing fraction"),
        FieldSpec(name="max_fg_miss_fraction", label="Max foreground missing fraction"),
        FieldSpec(name="max_miss_fraction", label="Max any-missing fraction"),
        FieldSpec(name="miss_pair", label="Enforce missing pairs", kind="bool"),
        Section("Batching logic (performance)"),
        FieldSpec(name="ct_discovery_batch_size", label="Discovery genes per task"),
        FieldSpec(name="ct_bootstrap_batch_size", label="Bootstrap genes per task"),
        Section("Publishing norms (debug)"),
        FieldSpec(name="publish_intermediates", label="Publish intermediate files", kind="bool"),
        FieldSpec(name="export_groups", label="Export groups (DEBUG)", kind="bool"),
        FieldSpec(name="export_perm_discovery", label="Export permuted discovery (DEBUG)", kind="bool"),
        Section("Contrast-selection tuning (conf/common.config)"),
        FieldSpec(name="pss_top_pct", label="PSS top percentile candidate gate"),
        FieldSpec(name="max_contrasts", label="Max contrasts (0 = dynamic)"),
        FieldSpec(name="min_contrasts", label="Minimum foreground contrasts"),
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

        self._ct_tool_label = QLabel("CT tools (--ct_tool)")
        self._essential_form.insertRow(0, self._ct_tool_label, self.ct_tool_discovery)
        self._essential_form.insertRow(1, "", self.ct_tool_resample)
        self._essential_form.insertRow(2, "", self.ct_tool_bootstrap)

    def retranslate(self, lang: str = "en") -> None:
        super().retranslate(lang)
        from gui.i18n import tr
        self._ct_tool_label.setText(tr("CT tools (--ct_tool)", lang))
        self.ct_tool_discovery.setText(tr("discovery", lang))
        self.ct_tool_resample.setText(tr("resample", lang))
        self.ct_tool_bootstrap.setText(tr("bootstrap", lang))

    def _on_ct_tool_discovery(self, value: bool) -> None:
        self._config.ct_tool_discovery = value
        self.changed.emit()

    def _on_ct_tool_resample(self, value: bool) -> None:
        self._config.ct_tool_resample = value
        self.changed.emit()

    def _on_ct_tool_bootstrap(self, value: bool) -> None:
        self._config.ct_tool_bootstrap = value
        self.changed.emit()
