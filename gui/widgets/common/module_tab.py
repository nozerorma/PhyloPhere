#!/usr/bin/env python3
# module_tab.py — ModuleTabWidget: reusable base class for all 9 module tabs.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Enable toggle + blurb + disclaimer + essential-fields group + fallback-inputs group
+ free-text advanced-flags box — built once from a ModuleTabSpec + FieldSpec list
(see gui/widgets/common/specs.py) rather than 9 separately hand-rolled layouts.

Enable-toggle behavior: the essential-fields group is disabled when the module is
off, while the fallback-inputs group becomes enabled at that same moment (inverted
enabled states driven by one QCheckBox). The tab itself is never hidden or
tab-disabled — the disclaimer and fallback fields are exactly what a user needs to
read/fill when a module is off (see implementation plan §5).
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

        blurb_label = QLabel(spec.blurb)
        blurb_label.setWordWrap(True)
        layout.addWidget(blurb_label)

        self.disclaimer_label = QLabel(f"If disabled: {spec.disclaimer}")
        self.disclaimer_label.setWordWrap(True)
        self.disclaimer_label.setStyleSheet("color: #8a6d3b;")  # amber note, not an error
        layout.addWidget(self.disclaimer_label)

        self.essential_group = QGroupBox("Essential fields")
        self._essential_form = QFormLayout(self.essential_group)
        layout.addWidget(self.essential_group)

        self.fallback_group = QGroupBox("Precomputed input (used when disabled)")
        self._fallback_form = QFormLayout(self.fallback_group)
        layout.addWidget(self.fallback_group)

        self.extra_flags_edit = QPlainTextEdit(config.extra_flags)
        self.extra_flags_edit.setPlaceholderText(
            "Extra Nextflow flags appended verbatim to this module's NF_FLAGS block"
        )
        self.extra_flags_edit.setMaximumHeight(60)
        self.extra_flags_edit.textChanged.connect(self._on_extra_flags_changed)
        layout.addWidget(QLabel("Advanced: extra Nextflow flags"))
        layout.addWidget(self.extra_flags_edit)

        layout.addStretch(1)

        self._field_widgets: dict[str, QWidget] = {}
        for f in spec.essential_fields:
            self._add_field(self._essential_form, f)
        for f in spec.fallback_fields:
            self._add_field(self._fallback_form, f)

        self._sync_group_enabled_state()

    # ── Field widget construction ────────────────────────────────────────────

    def _add_field(self, form: QFormLayout, f: FieldSpec) -> None:
        current_value = getattr(self._config, f.name)

        if f.kind == "bool":
            widget = QCheckBox()
            widget.setChecked(bool(current_value))
            widget.toggled.connect(lambda v, name=f.name: self._set_field(name, v))
        elif f.kind == "choice":
            widget = QComboBox()
            widget.addItems(list(f.choices))
            widget.setCurrentText(str(current_value))
            widget.currentTextChanged.connect(lambda v, name=f.name: self._set_field(name, v))
        elif f.kind in ("path_file", "path_dir"):
            widget = PathField(mode="file" if f.kind == "path_file" else "dir")
            widget.set_text(str(current_value))
            widget.textChanged.connect(lambda v, name=f.name: self._set_field(name, v))
        else:  # "str"
            widget = QLineEdit(str(current_value))
            widget.setPlaceholderText(f.placeholder)
            widget.textChanged.connect(lambda v, name=f.name: self._set_field(name, v))

        form.addRow(f.label, widget)
        self._field_widgets[f.name] = widget

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
        self.fallback_group.setEnabled(not enabled)
        self.disclaimer_label.setVisible(not enabled)

    def _on_extra_flags_changed(self) -> None:
        self._config.extra_flags = self.extra_flags_edit.toPlainText()
        self.changed.emit()
