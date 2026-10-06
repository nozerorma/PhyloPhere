#!/usr/bin/env python3
# caas_tab.py — CAAS / CT module tab (contrast selection: discovery, resample).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
CaasTab: field specification of the CAAS / contrast-selection module (conf/ct.config).

The discovery and resample sub-steps are two checkboxes on the config
(ct_tool_discovery, ct_tool_resample; see gui/models/modules.py) that jointly build
the comma-separated --ct_tool value (gui/generation/context.py). ModuleTabWidget
has no field kind for that, so CaasTab adds them directly in __init__.

--contrast_selection has no checkbox of its own: run_single.sh.j2 turns it on
whenever CAAS or Disambiguation is enabled, because contrast selection is the only
producer of the foreground/background trait file that both need.

Imported by: gui/widgets/main_window.py
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
        "Runs CAAStools' discovery and resample steps to find convergent "
        "amino-acid substitutions (CAAS) associated with the phenotype."
    ),
    disclaimer=(
        "Disambiguation and Accumulation need this module's output. Check Discovery/"
        "Resample on the Precomputed Run tab to feed them precomputed "
        "results instead."
    ),
    essential_fields=(
        Section("Discovery params"),
        FieldSpec(
            name="patterns",
            label="CT patterns",
            kind="multichoice",
            choices=("1", "2", "3", "4"),
            importance="default",
        ),
        FieldSpec(
            name="min_divergent_fraction",
            label="Min divergent fraction",
            kind="choice",
            choices=("0.5", "0.75", "1.0"),
            editable=True,
            importance="default",
        ),
        FieldSpec(name="caap_mode", label="CAAP mode (properties-based)", kind="bool", importance="default"),
        FieldSpec(name="multi_hypothesis", label="Multi-hypothesis mode", kind="bool", importance="default"),
        FieldSpec(
            name="max_fop",
            label="Max FOP hypotheses (H1..Hn)",
            placeholder="alternative Dunn-independent hypotheses per contrast (observed + null); default 100",
            importance="default",
        ),
        Section("Resample / permulation params"),
        FieldSpec(name="chunk_size", label="Resampled groups per output file", importance="optional"),
        FieldSpec(name="resample_use_n", label="Use sample size counts (n/c)", kind="bool", importance="default"),
        FieldSpec(
            name="perm_strategy",
            label="Permutation strategy",
            kind="choice",
            choices=("auto", "OU", "BM"),
            importance="default",
        ),
        FieldSpec(
            name="caas_full_perms",
            label="Permulations",
            placeholder="accepted permulations to harvest AND replay for the CAAS FCS null",
            importance="default",
        ),
        # Rated "optional" although borderline: it is a draw-budget cap (raised 50%
        # up to twice if the pool falls short) and does not redefine the null, but
        # a very low value can truncate the permulation pool.
        FieldSpec(
            name="max_tries",
            label="Max permulation tries",
            placeholder="draw budget; raised 50% up to twice if the pool falls short",
            importance="optional",
        ),
    ),
    advanced_fields=(
        Section("Missingness parameters in discovery mode and otherwise"),
        FieldSpec(
            name="caas_config_path",
            label="CAAS config file (auto-derived if empty)",
            kind="path_file",
            placeholder="auto-derived from contrast_selection.nf when empty",
            importance="optional",
        ),
        FieldSpec(name="max_bg_gaps_fraction", label="Max background gaps fraction", importance="default"),
        FieldSpec(name="max_fg_gaps_fraction", label="Max foreground gaps fraction", importance="default"),
        FieldSpec(name="max_gaps_fraction", label="Max any-gaps fraction", importance="default"),
        FieldSpec(name="max_bg_miss_fraction", label="Max background missing fraction", importance="default"),
        FieldSpec(name="max_fg_miss_fraction", label="Max foreground missing fraction", importance="default"),
        FieldSpec(name="max_miss_fraction", label="Max any-missing fraction", importance="default"),
        FieldSpec(name="miss_pair", label="Enforce missing pairs", kind="bool", importance="default"),
        Section("Batching logic (performance)"),
        FieldSpec(
            name="ct_core_batch_size",
            label="Permulation-null genes per task (replay + ASR replay)",
            importance="optional",
        ),
        Section("Permulation null"),
        FieldSpec(name="caas_perms_postproc", label="Apply post-processing filters to the permulation null", kind="bool", importance="default"),
        FieldSpec(name="caas_permulation_enrichment", label="Use CAAS permulations for the FCS Wilcoxon null", kind="bool", importance="default"),
        Section("Publishing norms (debug)"),
        FieldSpec(name="publish_intermediates", label="Publish intermediate files", kind="bool", importance="optional"),
        Section("Contrast-selection tuning (conf/common.config)"),
        FieldSpec(name="pss_top_pct", label="PSS top percentile candidate gate", importance="default"),
        FieldSpec(name="max_contrasts", label="Max contrasts (0 = dynamic)", importance="default"),
        FieldSpec(name="min_contrasts", label="Minimum foreground contrasts", importance="default"),
    ),
)


class CaasTab(ModuleTabWidget):
    def __init__(self, config: CaasConfig, parent=None):
        super().__init__(SPEC, config, parent)

        # discovery and resample checkboxes jointly build --ct_tool; inserted at
        # the top of the essential fields.
        self.ct_tool_discovery = QCheckBox("discovery")
        self.ct_tool_discovery.setChecked(config.ct_tool_discovery)
        self.ct_tool_discovery.toggled.connect(self._on_ct_tool_discovery)

        self.ct_tool_resample = QCheckBox("resample")
        self.ct_tool_resample.setChecked(config.ct_tool_resample)
        self.ct_tool_resample.toggled.connect(self._on_ct_tool_resample)

        self._ct_tool_label = QLabel("CT tools (--ct_tool)")
        self._essential_form.insertRow(0, self._ct_tool_label, self.ct_tool_discovery)
        self._essential_form.insertRow(1, "", self.ct_tool_resample)

    def retranslate(self, lang: str = "en") -> None:
        super().retranslate(lang)
        from gui.i18n import tr
        self._ct_tool_label.setText(tr("CT tools (--ct_tool)", lang))
        self.ct_tool_discovery.setText(tr("discovery", lang))
        self.ct_tool_resample.setText(tr("resample", lang))

    def _on_ct_tool_discovery(self, value: bool) -> None:
        self._config.ct_tool_discovery = value
        self.changed.emit()

    def _on_ct_tool_resample(self, value: bool) -> None:
        self._config.ct_tool_resample = value
        self.changed.emit()
