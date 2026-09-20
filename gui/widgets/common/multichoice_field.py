#!/usr/bin/env python3
# multichoice_field.py — row-of-checkboxes widget for a comma-separated string field.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QHBoxLayout, QWidget


class MultiChoiceField(QWidget):
    """One checkbox per choice, serialized to/from a comma-separated string
    (e.g. CAAStools' `patterns`, "1,2,3"). Selection order in the stored string
    always follows `choices`' own order, not click order.
    """

    valueChanged = Signal(str)

    def __init__(self, choices: tuple[str, ...], parent=None):
        super().__init__(parent)
        self._choices = list(choices)
        self._checkboxes: dict[str, QCheckBox] = {}

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for choice in self._choices:
            cb = QCheckBox(choice)
            cb.toggled.connect(self._on_toggled)
            layout.addWidget(cb)
            self._checkboxes[choice] = cb
        layout.addStretch(1)

    def _on_toggled(self, _checked: bool) -> None:
        self.valueChanged.emit(self.text())

    def text(self) -> str:
        return ",".join(c for c in self._choices if self._checkboxes[c].isChecked())

    def set_text(self, value: str) -> None:
        selected = {v.strip() for v in value.split(",") if v.strip()}
        for choice, cb in self._checkboxes.items():
            cb.blockSignals(True)
            cb.setChecked(choice in selected)
            cb.blockSignals(False)
