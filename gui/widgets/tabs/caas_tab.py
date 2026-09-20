#!/usr/bin/env python3
# caas_tab.py — CAAS / CT module tab (contrast selection: discovery, resample).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
The discovery/resample sub-steps are separate boolean checkboxes on the
config (ct_tool_discovery/resample, see gui/models/modules.py) rather than
FieldSpec entries, since ModuleTabWidget's field kinds don't cover "2 checkboxes
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
        FieldSpec(name="include_b0", label="Include main hypothesis (b0)", kind="bool", importance="default"),
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
        # Borderline default/optional: a draw-budget safety cap (auto-escalated 50%
        # up to twice if the pool falls short) rather than a knob that itself
        # redefines the null — but too low a value can silently truncate the
        # permulation pool, so it is not purely cosmetic either.
        FieldSpec(
            name="max_tries",
            label="Max permulation tries",
            placeholder="draw budget; raised 50% up to twice if the pool falls short",
            importance="optional",
        ),
        FieldSpec(
            name="caas_perms_fop",
            label="Mirror FOP harvest in the CAAS permulation null",
            kind="bool",
            importance="default",
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
        FieldSpec(name="ct_discovery_batch_size", label="Discovery genes per task", importance="optional"),
        FieldSpec(
            name="ct_perm_replay_batch_size",
            label="Permulation-null replay genes per batch",
            importance="optional",
        ),
        FieldSpec(
            name="ct_disambig_perms_batch_size",
            label="Permulation-null disambiguation genes per batch",
            importance="optional",
        ),
        Section("Publishing norms (debug)"),
        FieldSpec(name="publish_intermediates", label="Publish intermediate files", kind="bool", importance="optional"),
        FieldSpec(name="export_groups", label="Export groups (DEBUG)", kind="bool", importance="optional"),
        FieldSpec(name="export_perm_discovery", label="Export permuted discovery (DEBUG)", kind="bool", importance="optional"),
        Section("Contrast-selection tuning (conf/common.config)"),
        FieldSpec(name="pss_top_pct", label="PSS top percentile candidate gate", importance="default"),
        FieldSpec(name="max_contrasts", label="Max contrasts (0 = dynamic)", importance="default"),
        FieldSpec(name="min_contrasts", label="Minimum foreground contrasts", importance="default"),
    ),
)


class CaasTab(ModuleTabWidget):
    def __init__(self, config: CaasConfig, parent=None):
        super().__init__(SPEC, config, parent)

        # ct_tool discovery/resample: 2 checkboxes jointly building --ct_tool.
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
