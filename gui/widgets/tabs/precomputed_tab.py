#!/usr/bin/env python3
# precomputed_tab.py — Precomputed Run tab: base path + one reuse checkbox per stage.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
PrecomputedTab: one base path plus one reuse checkbox per pipeline stage.

The per-phenotype path of each input is base_path/<TRAIT>/..., so one base path
serves a batch with several phenotypes (gui/models/precomputed.py describes the
layout). Each checkbox both supplies a precomputed input and turns off the module
that would otherwise recompute it. The module is switched off by toggling the enable
checkbox on its own tab (self._module_tabs), not by writing its config field, so
that tab's enabled-state display (essential and advanced groups graying out) stays
consistent.

CT/CAAS has a general checkbox plus two specific ones (discovery, resample), as
the CAAS tab has two checkboxes for --ct_tool: the pipeline takes --discovery_from
and --resample_from independently, so reusing only one of them is meaningful.
Every other stage has a single checkbox; the choice of which files that implies is
made in gui/generation/templates/run_single.sh.j2 and is not exposed here.

Imported by: gui/widgets/main_window.py
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
        """The six *_tab arguments are the already-built module tabs (see
        MainWindow._build_project_tabs, which constructs this tab after them).
        Checking a box here calls that tab's own enable_toggle, so its config field
        and display update through its existing logic instead of this tab writing to
        another module's config."""
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
            # Post-processing has no enable toggle of its own to switch off:
            # gui/generation/context.py derives it from Disambiguation's toggle and
            # from this checkbox, so no module tab needs to be toggled here.
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

    # ── Sections ──────────────────────────────────────────────────────────────

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
        form.addRow("", self.use_discovery)
        form.addRow("", self.use_resample)
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

    # ── Slots ─────────────────────────────────────────────────────────────────

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
            # Checking the general box ticks both sub-boxes, since the usual case is a
            # fully precomputed CT stage; either can then be unticked.
            for cb in (self.use_discovery, self.use_resample):
                cb.setChecked(True)

    def _on_disambiguation_toggled(self, value: bool) -> None:
        self._toggle_module("disambiguation", value)

    def _toggle_module(self, key: str, value: bool) -> None:
        tab = self._module_tabs[key]
        if tab.enable_toggle.isChecked() != (not value):
            tab.enable_toggle.setChecked(not value)

    # ── i18n ──────────────────────────────────────────────────────────────────

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
