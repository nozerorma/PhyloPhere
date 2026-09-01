#!/usr/bin/env python3
# module_tab.py — ModuleTabWidget: reusable base class for all 9 module tabs.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Enable toggle + blurb + disclaimer + essential-fields group + collapsed-by-
default advanced-fields disclosure + free-text raw-flags box — built once from a
ModuleTabSpec + FieldSpec list (see gui/widgets/common/specs.py) rather than 9
separately hand-rolled layouts.

Enable-toggle behavior: the essential-fields group (and the advanced-fields
disclosure, if the spec has any) is disabled when the module is off. The tab
itself is never hidden or tab-disabled — the disclaimer is exactly what a user
needs to read when a module is off, pointing at the Precomputed Run tab
(gui/widgets/tabs/precomputed_tab.py) for the standalone/resume input that
substitutes for this module's output.

No tab currently populates its own `fallback_fields` (every standalone/resume
input now lives on the Precomputed Run tab instead) — the mechanism stays
available for a future module-specific case, but the group box housing it is
only built when a spec's `fallback_fields` is actually non-empty, so today no
tab shows an empty, confusingly-enabled box for it.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import ModuleConfigBase
from gui.widgets.common.collapsible import CollapsibleSection
from gui.widgets.common.path_field import PathField
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec


class ModuleTabWidget(QWidget):
    changed = Signal()

    def __init__(self, spec: ModuleTabSpec, config: ModuleConfigBase, parent=None):
        super().__init__(parent)
        self._spec = spec
        self._config = config

        layout = QVBoxLayout(self)

        self.enable_toggle = QCheckBox(f"Enable {spec.title}")
        self.enable_toggle.setChecked(getattr(config, spec.enabled_field))
        self.enable_toggle.toggled.connect(self._on_enable_toggled)
        layout.addWidget(self.enable_toggle)

        self.blurb_label = QLabel(spec.blurb)
        self.blurb_label.setWordWrap(True)
        layout.addWidget(self.blurb_label)

        self.disclaimer_label = QLabel(f"If disabled: {spec.disclaimer}")
        self.disclaimer_label.setWordWrap(True)
        self.disclaimer_label.setStyleSheet(
            "QLabel { background-color: rgba(234, 179, 8, 0.12); color: #d97706; "
            "border: 1px solid rgba(234, 179, 8, 0.3); border-radius: 6px; padding: 8px 12px; font-weight: 500; }"
        )
        layout.addWidget(self.disclaimer_label)

        self.essential_group = QGroupBox("Essential fields")
        self._essential_form = QFormLayout(self.essential_group)
        layout.addWidget(self.essential_group)

        # Fine-tuning knobs most runs leave at their conf/*.config default, tucked
        # behind a collapsed-by-default disclosure so they don't read as equally
        # important to the essential fields above. Only built when the spec
        # actually lists some (a few small tabs — e.g. VEP — have none).
        self.advanced_section: CollapsibleSection | None = None
        self._advanced_form: QFormLayout | None = None
        if spec.advanced_fields:
            advanced_box = QGroupBox()
            advanced_box.setFlat(True)
            self._advanced_form = QFormLayout(advanced_box)
            n_adv = len([f for f in spec.advanced_fields if f.kind != "section"])
            self.advanced_section = CollapsibleSection(
                f"Advanced parameters ({n_adv})", advanced_box
            )
            layout.addWidget(self.advanced_section)

        # Only built when the spec actually lists module-specific fallback fields
        # (none do today — see this file's module docstring); avoids an empty,
        # always-enabled-when-disabled box with nothing in it.
        self.fallback_group: QGroupBox | None = None
        self._fallback_form: QFormLayout | None = None
        if spec.fallback_fields:
            self.fallback_group = QGroupBox("Precomputed input (used when disabled)")
            self._fallback_form = QFormLayout(self.fallback_group)
            layout.addWidget(self.fallback_group)

        self.extra_flags_edit = QPlainTextEdit(config.extra_flags)
        self.extra_flags_edit.setPlaceholderText(
            "Extra Nextflow flags appended verbatim to this module's NF_FLAGS block"
        )
        self.extra_flags_edit.setMaximumHeight(60)
        self.extra_flags_edit.textChanged.connect(self._on_extra_flags_changed)
        layout.addWidget(QLabel("Power-user override: raw extra Nextflow flags"))
        layout.addWidget(self.extra_flags_edit)

        layout.addStretch(1)

        self._field_widgets: dict[str, QWidget] = {}
        self._label_widgets: list[tuple[QLabel, str]] = []
        for f in spec.essential_fields:
            self._add_field(self._essential_form, f)
        for f in spec.advanced_fields:
            self._add_field(self._advanced_form, f)
        for f in spec.fallback_fields:
            self._add_field(self._fallback_form, f)

        self._sync_group_enabled_state()

    # ── Field widget construction ────────────────────────────────────────────

    def _add_field(self, form: QFormLayout, f: FieldSpec) -> None:
        if f.kind == "section":
            header_lbl = QLabel(f"// {f.label}")
            header_lbl.setStyleSheet(
                "QLabel { font-weight: bold; color: #475569; padding-top: 10px; padding-bottom: 2px; border-bottom: 1px solid rgba(0,0,0,0.12); margin-top: 4px; }"
            )
            form.addRow(header_lbl)
            self._label_widgets.append((header_lbl, f"// {f.label}"))
            return

        current_value = getattr(self._config, f.name)

        if f.kind == "bool":
            widget = QCheckBox(f.label)
            widget.setChecked(bool(current_value))
            widget.toggled.connect(lambda v, name=f.name: self._set_field(name, v))
            form.addRow(widget)
            self._label_widgets.append((widget, f.label))
        elif f.kind == "choice":
            widget = QComboBox()
            widget.addItems(list(f.choices))
            widget.setCurrentText(str(current_value))
            widget.currentTextChanged.connect(lambda v, name=f.name: self._set_field(name, v))
            label_lbl = QLabel(f.label)
            form.addRow(label_lbl, widget)
            self._label_widgets.append((label_lbl, f.label))
        elif f.kind in ("path_file", "path_dir"):
            widget = PathField(mode="file" if f.kind == "path_file" else "dir")
            widget.set_text(str(current_value))
            widget.textChanged.connect(lambda v, name=f.name: self._set_field(name, v))
            label_lbl = QLabel(f.label)
            form.addRow(label_lbl, widget)
            self._label_widgets.append((label_lbl, f.label))
        else:  # "str"
            widget = QLineEdit(str(current_value))
            widget.setPlaceholderText(f.placeholder)
            widget.textChanged.connect(lambda v, name=f.name: self._set_field(name, v))
            label_lbl = QLabel(f.label)
            form.addRow(label_lbl, widget)
            self._label_widgets.append((label_lbl, f.label))

        self._field_widgets[f.name] = widget

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        self.enable_toggle.setText(f"{tr('Enable', lang)} {tr(self._spec.title, lang)}")
        if hasattr(self, "blurb_label") and self.blurb_label:
            self.blurb_label.setText(tr(self._spec.blurb, lang))
        if hasattr(self, "disclaimer_label") and self.disclaimer_label:
            self.disclaimer_label.setText(f"{tr('If disabled', lang)}: {tr(self._spec.disclaimer, lang)}")
        if hasattr(self, "essential_group") and self.essential_group:
            self.essential_group.setTitle(tr("Essential fields", lang))
        if hasattr(self, "advanced_section") and self.advanced_section:
            self.advanced_section.set_title(f"{tr('Advanced parameters', lang)} ({len(self._spec.advanced_fields)})")
        if hasattr(self, "fallback_group") and self.fallback_group:
            self.fallback_group.setTitle(tr("Precomputed input (used when disabled)", lang))
        for label_widget, orig_text in getattr(self, "_label_widgets", []):
            if isinstance(label_widget, QCheckBox):
                label_widget.setText(tr(orig_text, lang))
            elif isinstance(label_widget, QLabel):
                label_widget.setText(tr(orig_text, lang))

    def _set_field(self, name: str, value) -> None:
        setattr(self._config, name, value)
        self.changed.emit()

    # ── Enable toggle ────────────────────────────────────────────────────────

    def _on_enable_toggled(self, value: bool) -> None:
        setattr(self._config, self._spec.enabled_field, value)
        self._sync_group_enabled_state()
        self.changed.emit()

    def _sync_group_enabled_state(self) -> None:
        enabled = getattr(self._config, self._spec.enabled_field)
        self.essential_group.setEnabled(enabled)
        if self.advanced_section is not None:
            self.advanced_section.setEnabled(enabled)
        if self.fallback_group is not None:
            self.fallback_group.setEnabled(not enabled)
        self.disclaimer_label.setVisible(not enabled)

    def _on_extra_flags_changed(self) -> None:
        self._config.extra_flags = self.extra_flags_edit.toPlainText()
        self.changed.emit()
