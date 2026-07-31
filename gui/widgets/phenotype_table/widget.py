#!/usr/bin/env python3
# widget.py — QTableView + add/remove controls for the phenotype catalogue.
# PhyloPhere | gui/widgets/phenotype_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QPushButton, QTableView, QVBoxLayout, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.runtime import PhenotypeRow
from gui.widgets.phenotype_table.model import PhenotypeTableModel


class PhenotypeTableWidget(QWidget):
    """QTableView + Add row / Remove row controls."""

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

        self.add_btn = QPushButton("Add row")
        self.remove_btn = QPushButton("Remove selected")
        self.add_btn.clicked.connect(self._add_row)
        self.remove_btn.clicked.connect(self._remove_selected)

        button_row = QHBoxLayout()
        button_row.addWidget(self.add_btn)
        button_row.addWidget(self.remove_btn)
        button_row.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.table)
        layout.addLayout(button_row)

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "add_btn"):
            self.add_btn.setText(tr("Add row", lang))
        if hasattr(self, "remove_btn"):
            self.remove_btn.setText(tr("Remove selected", lang))
        if hasattr(self, "model") and hasattr(self.model, "set_lang"):
            self.model.set_lang(lang)

    def _add_row(self) -> None:
        self.model.insertRows(self.model.rowCount(), 1)

    def _selected_row_index(self) -> int | None:
        indexes = self.table.selectionModel().selectedRows()
        return indexes[0].row() if indexes else None

    def _remove_selected(self) -> None:
        row = self._selected_row_index()
        if row is not None:
            self.model.removeRows(row, 1)
