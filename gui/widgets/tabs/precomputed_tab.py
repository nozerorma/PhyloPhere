#!/usr/bin/env python3
# precomputed_tab.py — Precomputed Run tab: base path + one reuse checkbox per stage.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Replaces the old one-global-path-per-field design (gui/models/precomputed.py's
module docstring has the full story — a single global path can't work once a batch
run has more than one phenotype, since TRAIT varies per row). Every checkbox here
both supplies a precomputed input AND turns off the module that would otherwise
recompute it live — done by directly toggling that module's own enable checkbox on
its tab (self._module_tabs / self._postproc_checkbox below), not by poking its
config field, so the other tab's own enabled-state visuals (its essential/advanced
groups graying out) stay in sync for free.

CT/CAAS gets a general checkbox + 3 specific ones (discovery/resample/bootstrap),
mirroring the CAAS tab's own 3-checkbox --ct_tool pattern — ct.nf treats each of
discovery_from/resample_from/bootstrap_from independently, so partial reuse (e.g.
discovery + bootstrap precomputed, resample recomputed) is meaningful. Every other
stage is one checkbox: the actual "which of several files" fan-out (e.g. Disambiguation
still needing 2 files) is an implementation detail handled entirely in
gui/generation/templates/run_single.sh.j2's path construction, not exposed here.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.precomputed import PrecomputedConfig
from gui.widgets.common.path_field import PathField


class PrecomputedTab(QWidget):
    changed = Signal()

    def __init__(
        self,
        config: PrecomputedConfig,
        *,
        caas_tab,
        disambiguation_tab,
        accumulation_tab,
        rer_tab,
        fade_tab,
        vep_tab,
        parent=None,
    ):
        """The 6 *_tab args are the already-built module tabs (see
        MainWindow._build_project_tabs, which constructs this tab last) — checking a
        box here calls straight into that tab's own enable_toggle, so its config
        field and visuals update through its existing logic rather than this tab
        reaching into another module's config directly."""
        super().__init__(parent)
        self._config = config
        self._module_tabs = {
            "ct": caas_tab,
            "disambiguation": disambiguation_tab,
            "accumulation": accumulation_tab,
            "rer": rer_tab,
            "fade": fade_tab,
            "vep": vep_tab,
        }

        self._group_boxes: list[tuple[QGroupBox, str]] = []
        self._checkbox_labels: list[tuple[QCheckBox, str]] = []

        outer = QVBoxLayout(self)

        self.note_label = QLabel(
            "One base path, reused for every checked box below: base_path/<TRAIT>/... "
            "(each phenotype's own subdirectory, matching a prior completed run's "
            "output layout). Check a box to feed that stage's already-computed output "
            "in instead of recomputing it — doing so also switches that stage off."
        )
        self.note_label.setWordWrap(True)
        outer.addWidget(self.note_label)

        self.base_path_field = PathField(mode="dir")
        self.base_path_field.set_text(self._config.base_path)
        self.base_path_field.textChanged.connect(self._on_base_path_changed)
        base_form = QFormLayout()
        self.base_path_label = QLabel("Base path")
        base_form.addRow(self.base_path_label, self.base_path_field)
        outer.addLayout(base_form)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll, stretch=1)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.addWidget(self._build_ct_section())
        content_layout.addWidget(self._build_simple_section(
            "Disambiguation", "use_disambiguation", self._on_disambiguation_toggled
        ))
        content_layout.addWidget(self._build_simple_section(
            # Post-processing has no separate enable toggle of its own to flip off
            # here — it's derived straight from Disambiguation's enable_toggle in
            # gui/generation/context.py, so checking this box just supplies the
            # precomputed files without any module-toggle cascade.
            "Post-processing", "use_postproc", lambda v: None
        ))
        content_layout.addWidget(self._build_simple_section(
            "Accumulation", "use_accumulation", lambda v: self._toggle_module("accumulation", v)
        ))
        content_layout.addWidget(self._build_simple_section(
            "RERconverge", "use_rer", lambda v: self._toggle_module("rer", v)
        ))
        content_layout.addWidget(self._build_simple_section(
            "FADE", "use_fade", lambda v: self._toggle_module("fade", v)
        ))
        content_layout.addWidget(self._build_simple_section(
            "VEP", "use_vep", lambda v: self._toggle_module("vep", v)
        ))
        content_layout.addStretch(1)
        scroll.setWidget(content)

    # ── Sections ─────────────────────────────────────────────────────────────

    def _build_ct_section(self) -> QGroupBox:
        box = QGroupBox("CT / CAAS")
        self._group_boxes.append((box, "CT / CAAS"))
        layout = QVBoxLayout(box)

        self.use_ct = QCheckBox("Use precomputed CT / CAAS output")
        self.use_ct.setChecked(self._config.use_ct)
        self.use_ct.toggled.connect(self._on_ct_toggled)
        layout.addWidget(self.use_ct)

        form = QFormLayout()
        self.use_discovery = QCheckBox("Discovery")
        self.use_discovery.setChecked(self._config.use_discovery)
        self.use_discovery.toggled.connect(lambda v: self._set_bool("use_discovery", v))
        self.use_resample = QCheckBox("Resample")
        self.use_resample.setChecked(self._config.use_resample)
        self.use_resample.toggled.connect(lambda v: self._set_bool("use_resample", v))
        self.use_bootstrap = QCheckBox("Bootstrap")
        self.use_bootstrap.setChecked(self._config.use_bootstrap)
        self.use_bootstrap.toggled.connect(lambda v: self._set_bool("use_bootstrap", v))
        form.addRow("", self.use_discovery)
        form.addRow("", self.use_resample)
        form.addRow("", self.use_bootstrap)
        layout.addLayout(form)

        return box

    def _build_simple_section(self, title: str, field_name: str, on_toggled) -> QGroupBox:
        box = QGroupBox(title)
        self._group_boxes.append((box, title))
        layout = QVBoxLayout(box)
        checkbox = QCheckBox(f"Use precomputed {title} output")
        checkbox.setChecked(getattr(self._config, field_name))
        checkbox.toggled.connect(lambda v, n=field_name: self._on_checkbox(n, v, on_toggled))
        self._checkbox_labels.append((checkbox, f"Use precomputed {title} output"))
        layout.addWidget(checkbox)
        setattr(self, f"_checkbox_{field_name}", checkbox)
        return box

    # ── Slots ────────────────────────────────────────────────────────────────

    def _on_base_path_changed(self, value: str) -> None:
        self._config.base_path = value
        self.changed.emit()

    def _set_bool(self, name: str, value: bool) -> None:
        setattr(self._config, name, value)
        self.changed.emit()

    def _on_checkbox(self, name: str, value: bool, on_toggled) -> None:
        self._set_bool(name, value)
        on_toggled(value)

    def _on_ct_toggled(self, value: bool) -> None:
        self._set_bool("use_ct", value)
        self._toggle_module("ct", value)
        if value:
            # Convenience default: checking the general box reuses all 3 files: the
            # common case is a full precomputed CT stage, not a partial one.
            for cb in (self.use_discovery, self.use_resample, self.use_bootstrap):
                cb.setChecked(True)

    def _on_disambiguation_toggled(self, value: bool) -> None:
        self._toggle_module("disambiguation", value)

    def _toggle_module(self, key: str, value: bool) -> None:
        tab = self._module_tabs[key]
        if tab.enable_toggle.isChecked() != (not value):
            tab.enable_toggle.setChecked(not value)

    # ── i18n ─────────────────────────────────────────────────────────────────

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "note_label"):
            self.note_label.setText(tr(
                "One base path, reused for every checked box below: base_path/<TRAIT>/... "
                "(each phenotype's own subdirectory, matching a prior completed run's "
                "output layout). Check a box to feed that stage's already-computed output "
                "in instead of recomputing it — doing so also switches that stage off.",
                lang,
            ))
        if hasattr(self, "base_path_label"):
            self.base_path_label.setText(tr("Base path", lang))
        for box, title in self._group_boxes:
            box.setTitle(tr(title, lang))
        for checkbox, orig_text in self._checkbox_labels:
            checkbox.setText(tr(orig_text, lang))
        if hasattr(self, "use_ct"):
            self.use_ct.setText(tr("Use precomputed CT / CAAS output", lang))
        if hasattr(self, "use_discovery"):
            self.use_discovery.setText(tr("Discovery", lang))
        if hasattr(self, "use_resample"):
            self.use_resample.setText(tr("Resample", lang))
        if hasattr(self, "use_bootstrap"):
            self.use_bootstrap.setText(tr("Bootstrap", lang))
