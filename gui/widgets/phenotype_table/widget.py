#!/usr/bin/env python3
# widget.py — QTableView + add/remove/overrides controls for the phenotype catalogue.
# PhyloPhere | gui/widgets/phenotype_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
The 3 per-row Scoring fallback overrides (scoring_rer_input, scoring_fade_summary_top/
bottom) aren't visible columns — they're edited via a small dialog opened from the
"Row overrides..." button, per the implementation plan (§4/§5): most rows won't need
them (only relevant when RER/FADE are disabled), so a dedicated dialog keeps the
main table from growing 3 more mostly-empty columns.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.runtime import PhenotypeRow
from gui.widgets.common.path_field import PathField
from gui.widgets.phenotype_table.model import PhenotypeTableModel


class RowOverridesDialog(QDialog):
    """Edits a single PhenotypeRow's per-row Scoring fallback overrides."""

    def __init__(self, row: PhenotypeRow, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Precomputed overrides — {row.trait or '(untitled)'}")
        self._row = row

        self.rer_input = PathField(mode="file")
        self.rer_input.set_text(row.scoring_rer_input)
        self.fade_top = PathField(mode="file")
        self.fade_top.set_text(row.scoring_fade_summary_top)
        self.fade_bottom = PathField(mode="file")
        self.fade_bottom.set_text(row.scoring_fade_summary_bottom)

        form = QFormLayout()
        form.addRow("RERconverge summary (scoring_rer_input)", self.rer_input)
        form.addRow("FADE summary — top (scoring_fade_summary_top)", self.fade_top)
        form.addRow("FADE summary — bottom (scoring_fade_summary_bottom)", self.fade_bottom)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def apply_to_row(self) -> None:
        self._row.scoring_rer_input = self.rer_input.text().strip()
        self._row.scoring_fade_summary_top = self.fade_top.text().strip()
        self._row.scoring_fade_summary_bottom = self.fade_bottom.text().strip()


class PhenotypeTableWidget(QWidget):
    """QTableView + Add row / Remove row / Row overrides... controls."""

    changed = Signal()

    def __init__(self, rows: list[PhenotypeRow], parent=None):
        super().__init__(parent)
        self.model = PhenotypeTableModel(rows, self)
        self.model.dataChanged.connect(self.changed)
        self.model.rowsInserted.connect(self.changed)
        self.model.rowsRemoved.connect(self.changed)

        self.table = QTableView(self)
        self.table.setModel(self.model)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)

        add_btn = QPushButton("Add row")
        remove_btn = QPushButton("Remove selected")
        overrides_btn = QPushButton("Row overrides...")
        add_btn.clicked.connect(self._add_row)
        remove_btn.clicked.connect(self._remove_selected)
        overrides_btn.clicked.connect(self._edit_overrides)

        button_row = QHBoxLayout()
        button_row.addWidget(add_btn)
        button_row.addWidget(remove_btn)
        button_row.addWidget(overrides_btn)
        button_row.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addLayout(button_row)
        layout.addWidget(self.table)

    def _add_row(self) -> None:
        self.model.insertRows(self.model.rowCount(), 1)

    def _selected_row_index(self) -> int | None:
        indexes = self.table.selectionModel().selectedRows()
        return indexes[0].row() if indexes else None

    def _remove_selected(self) -> None:
        row = self._selected_row_index()
        if row is not None:
            self.model.removeRows(row, 1)

    def _edit_overrides(self) -> None:
        row_index = self._selected_row_index()
        if row_index is None:
            return
        row = self.model._rows[row_index]
        dialog = RowOverridesDialog(row, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            dialog.apply_to_row()
            self.changed.emit()
