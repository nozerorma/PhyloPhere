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
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="CAAS / Contrast Selection",
    blurb=(
        "Runs CAAStools' discovery, resample, and bootstrap steps to find convergent "
        "amino-acid substitutions (CAAS) associated with the phenotype."
    ),
    disclaimer=(
        "Disambiguation and Accumulation need this module's output. Supply "
        "discovery_from/resample_from/bootstrap_from on the Precomputed Run tab to "
        "feed them precomputed results instead."
    ),
    essential_fields=(
        FieldSpec(name="caas_config_path", label="CAAS config file", kind="path_file"),
        FieldSpec(name="perms_cycles", label="Permulation cycles"),
        FieldSpec(name="caas_full_perms", label="Full-pool permulations"),
        FieldSpec(
            name="caas_permulation_enrichment", label="Permulation-excess for enrichment", kind="bool"
        ),
    ),
    advanced_fields=(
        FieldSpec(name="patterns", label="CT patterns"),
        # Contrast-selection tuning (conf/common.config)
        FieldSpec(name="top_quantile", label="Top quantile (parameterized method)"),
        FieldSpec(name="bottom_quantile", label="Bottom quantile (parameterized method)"),
        FieldSpec(name="contrast_max_iter", label="Contrast max iterations"),
        FieldSpec(name="min_contrasts", label="Minimum foreground contrasts"),
        # Discovery/resample fine-tuning (conf/ct.config)
        FieldSpec(name="publish_intermediates", label="Publish intermediate files", kind="bool"),
        FieldSpec(name="ct_discovery_batch_size", label="Discovery genes per task"),
        FieldSpec(name="ct_bootstrap_batch_size", label="Bootstrap genes per task"),
        FieldSpec(name="min_divergent_fraction", label="Min divergent fraction"),
        FieldSpec(name="max_bg_gaps_fraction", label="Max background gaps fraction"),
        FieldSpec(name="max_fg_gaps_fraction", label="Max foreground gaps fraction"),
        FieldSpec(name="max_gaps_fraction", label="Max any-gaps fraction"),
        FieldSpec(name="max_bg_miss_fraction", label="Max background missing fraction"),
        FieldSpec(name="max_fg_miss_fraction", label="Max foreground missing fraction"),
        FieldSpec(name="max_miss_fraction", label="Max any-missing fraction"),
        FieldSpec(name="miss_pair", label="Enforce missing pairs", kind="bool"),
        FieldSpec(name="caap_mode", label="CAAP mode (properties-based)", kind="bool"),
        FieldSpec(name="fgsize", label="Foreground size (strategy random)"),
        FieldSpec(name="bgsize", label="Background size (strategy random)"),
        FieldSpec(name="perm_strategy", label="Permutation strategy", kind="choice", choices=("FGBG", "BM")),
        FieldSpec(name="traitvalues", label="Trait values file (strategy BM)", kind="path_file"),
        FieldSpec(name="chunk_size", label="Resampled groups per output file"),
        FieldSpec(name="include_b0", label="Include main hypothesis (b0)", kind="bool"),
        FieldSpec(name="alpha_threshold", label="Alpha threshold (significance)"),
        FieldSpec(name="export_groups", label="Export groups (DEBUG)", kind="bool"),
        FieldSpec(name="export_perm_discovery", label="Export permuted discovery (DEBUG)", kind="bool"),
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
