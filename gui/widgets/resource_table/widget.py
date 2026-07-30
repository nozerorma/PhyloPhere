#!/usr/bin/env python3
# widget.py — QTableView + add/remove/preset controls for per-process resource overrides.
# PhyloPhere | gui/widgets/resource_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
The 3 preset buttons wholesale-replace every row from one of the checked-in
conf/resources.config.{slurm,local,local_lowspec} files (see
gui/resource_presets.py) — a starting point to hand-tune from, not a locked
choice. Selector Type is edited via a combo box delegate rather than free text,
since it only has two legal values (see ResourceOverrideTableModel.setData).
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
from gui import resource_presets
from gui.models.resources import ProcessResourceOverride
from gui.widgets.resource_table.model import ResourceOverrideTableModel

_PRESET_BUTTONS = [
    ("slurm", "Load SLURM defaults"),
    ("local", "Load local defaults (32cpu/64GB)"),
    ("local_lowspec", "Load local low-spec defaults (8cpu/16GB)"),
]


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
    """QTableView + Add row / Remove row / 3 preset-loading buttons."""

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

        preset_row = QHBoxLayout()
        self.preset_btns: list[tuple[QPushButton, str, str]] = []
        for name, label in _PRESET_BUTTONS:
            btn = QPushButton(label)
            btn.clicked.connect(lambda _checked=False, n=name: self._load_preset(n))
            preset_row.addWidget(btn)
            self.preset_btns.append((btn, name, label))
        preset_row.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addLayout(preset_row)
        layout.addLayout(button_row)
        layout.addWidget(self.table)

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "add_btn"):
            self.add_btn.setText(tr("Add row", lang))
        if hasattr(self, "remove_btn"):
            self.remove_btn.setText(tr("Remove selected", lang))
        for btn, name, orig_label in getattr(self, "preset_btns", []):
            btn.setText(tr(orig_label, lang))

    def _add_row(self) -> None:
        self.model.insertRows(self.model.rowCount(), 1)

    def _remove_selected(self) -> None:
        indexes = self.table.selectionModel().selectedRows()
        for index in sorted((i.row() for i in indexes), reverse=True):
            self.model.removeRows(index, 1)

    def _load_preset(self, name: str) -> None:
        self.model.replace_all(resource_presets.load_preset(name))
