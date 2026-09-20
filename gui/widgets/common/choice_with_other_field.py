#!/usr/bin/env python3
# choice_with_other_field.py — curated dropdown + free-text "other" fallback.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLineEdit, QWidget


class ChoiceWithOtherField(QWidget):
    """QComboBox over curated display labels, plus a free-text field revealed
    only when the last label (the "other" sentinel, e.g. "Other") is picked.

    `choices` are display labels; `values` gives the stored value for every
    label except the last, index-aligned with `choices[:-1]` (e.g. string_species:
    choices=("Human", "Mouse", "Other"), values=("9606", "10090")). The stored
    value while "Other" is selected is whatever the user types into the free-text
    field — never the label "Other" itself.
    """

    valueChanged = Signal(str)

    def __init__(self, choices: tuple[str, ...], values: tuple[str, ...], parent=None):
        super().__init__(parent)
        self._choices = list(choices)
        self._values = list(values)
        self._other_label = self._choices[-1] if self._choices else "Other"

        self.combo = QComboBox()
        self.combo.addItems(self._choices)
        self.combo.currentTextChanged.connect(self._on_combo_changed)

        self.other_edit = QLineEdit()
        self.other_edit.setVisible(False)
        self.other_edit.textChanged.connect(lambda _: self.valueChanged.emit(self.text()))

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.combo)
        layout.addWidget(self.other_edit, stretch=1)

    def _on_combo_changed(self, label: str) -> None:
        self.other_edit.setVisible(label == self._other_label)
        self.valueChanged.emit(self.text())

    def text(self) -> str:
        label = self.combo.currentText()
        if label == self._other_label:
            return self.other_edit.text()
        try:
            return self._values[self._choices.index(label)]
        except (ValueError, IndexError):
            return label

    def set_text(self, value: str) -> None:
        if value in self._values:
            label = self._choices[self._values.index(value)]
            self.combo.blockSignals(True)
            self.combo.setCurrentText(label)
            self.combo.blockSignals(False)
            self.other_edit.setVisible(False)
        else:
            self.combo.blockSignals(True)
            self.combo.setCurrentText(self._other_label)
            self.combo.blockSignals(False)
            self.other_edit.setVisible(True)
            self.other_edit.blockSignals(True)
            self.other_edit.setText(value)
            self.other_edit.blockSignals(False)
