#!/usr/bin/env python3
# widget.py — QTableView + add/remove/load-defaults controls for per-process resource overrides.
# PhyloPhere | gui/widgets/resource_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
ResourceOverrideTableWidget: table of per-process resource overrides with Add row,
Remove selected and load-defaults buttons.

The table starts empty: conf/resources.config is the source of defaults, and each
row is a deliberate deviation from it. The load button replaces every row with the
current defaults parsed from that file (see gui/resource_defaults.py), as a
starting point to edit. Selector Type is edited through a combo-box delegate
because it has only two legal values (see ResourceOverrideTableModel.setData).

Imported by: gui/widgets/tabs/resources_tab.py
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QPushButton,
    QStyledItemDelegate,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from PySide6.QtCore import Signal

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import resource_defaults
from gui.models.resources import ProcessResourceOverride
from gui.widgets.resource_table.model import ResourceOverrideTableModel

_DEFAULTS_BUTTON = "Load conf defaults into the table"


class SelectorTypeDelegate(QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        combo = QComboBox(parent)
        combo.addItems(["withName", "withLabel"])
        return combo

    def setEditorData(self, editor, index):
        editor.setCurrentText(index.data())

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText())


class ResourceOverrideTableWidget(QWidget):
    """QTableView + Add row / Remove selected / load-conf-defaults buttons."""

    changed = Signal()

    def __init__(self, rows: list[ProcessResourceOverride], parent=None):
        super().__init__(parent)
        self.model = ResourceOverrideTableModel(rows, self)
        self.model.dataChanged.connect(self.changed)
        self.model.rowsInserted.connect(self.changed)
        self.model.rowsRemoved.connect(self.changed)
        self.model.modelReset.connect(self.changed)

        self.table = QTableView(self)
        self.table.setModel(self.model)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setItemDelegateForColumn(0, SelectorTypeDelegate(self.table))

        self.add_btn = QPushButton("Add row")
        self.remove_btn = QPushButton("Remove selected")
        self.add_btn.clicked.connect(self._add_row)
        self.remove_btn.clicked.connect(self._remove_selected)

        button_row = QHBoxLayout()
        button_row.addWidget(self.add_btn)
        button_row.addWidget(self.remove_btn)
        button_row.addStretch(1)

        defaults_row = QHBoxLayout()
        self.defaults_btn = QPushButton(_DEFAULTS_BUTTON)
        self.defaults_btn.clicked.connect(self._load_defaults)
        defaults_row.addWidget(self.defaults_btn)
        defaults_row.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addLayout(defaults_row)
        layout.addLayout(button_row)
        layout.addWidget(self.table)

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "add_btn"):
            self.add_btn.setText(tr("Add row", lang))
        if hasattr(self, "remove_btn"):
            self.remove_btn.setText(tr("Remove selected", lang))
        if hasattr(self, "defaults_btn"):
            self.defaults_btn.setText(tr(_DEFAULTS_BUTTON, lang))

    def _add_row(self) -> None:
        self.model.insertRows(self.model.rowCount(), 1)

    def _remove_selected(self) -> None:
        indexes = self.table.selectionModel().selectedRows()
        for index in sorted((i.row() for i in indexes), reverse=True):
            self.model.removeRows(index, 1)

    def _load_defaults(self) -> None:
        self.model.replace_all(resource_defaults.load_defaults())
