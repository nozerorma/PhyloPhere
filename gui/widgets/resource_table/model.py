#!/usr/bin/env python3
# model.py — QAbstractTableModel for the Resources tab's per-process override table.
# PhyloPhere | gui/widgets/resource_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Backed directly by ResourcesConfig.process_overrides: list[ProcessResourceOverride]
— table edits flow straight into the serializable model, same pattern as
gui/widgets/phenotype_table/model.py. selector_type is edited via a combo
delegate set up by the widget, not free text (see widget.py).
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import fields

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.resources import ProcessResourceOverride

COLUMNS = [
    ("selector_type", "Type"),
    ("selector", "Process selector"),
    ("cpus", "cpus"),
    ("memory", "memory"),
]


class ResourceOverrideTableModel(QAbstractTableModel):
    def __init__(self, rows: list[ProcessResourceOverride], parent=None):
        super().__init__(parent)
        self._rows = rows  # shared reference into ResourcesConfig.process_overrides

    # ── Required overrides ─────────────────────────────────────────────────

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation != Qt.Orientation.Horizontal:
            if role == Qt.ItemDataRole.DisplayRole:
                return str(section + 1)
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return COLUMNS[section][1]
        return None

    def _field_name(self, column: int) -> str:
        return COLUMNS[column][0]

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
            row = self._rows[index.row()]
            return getattr(row, self._field_name(index.column()))
        return None

    def setData(self, index: QModelIndex, value, role=Qt.ItemDataRole.EditRole) -> bool:
        if role != Qt.ItemDataRole.EditRole or not index.isValid():
            return False
        row = self._rows[index.row()]
        field_name = self._field_name(index.column())
        if field_name == "selector_type" and value not in ("withName", "withLabel"):
            return False
        setattr(row, field_name, value)
        self.dataChanged.emit(index, index)
        return True

    def flags(self, index: QModelIndex):
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        return Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsEditable

    # ── Row add/remove/replace ──────────────────────────────────────────────

    def insertRows(self, row: int, count: int, parent=QModelIndex()) -> bool:
        self.beginInsertRows(parent, row, row + count - 1)
        for i in range(count):
            self._rows.insert(row + i, ProcessResourceOverride())
        self.endInsertRows()
        return True

    def removeRows(self, row: int, count: int, parent=QModelIndex()) -> bool:
        if row < 0 or row + count > len(self._rows):
            return False
        self.beginRemoveRows(parent, row, row + count - 1)
        del self._rows[row : row + count]
        self.endRemoveRows()
        return True

    def replace_all(self, new_rows: list[ProcessResourceOverride]) -> None:
        """Wholesale swap used by the preset-loading buttons."""
        self.beginResetModel()
        self._rows[:] = new_rows
        self.endResetModel()


assert {name for name, _ in COLUMNS} == {f.name for f in fields(ProcessResourceOverride)}, (
    "COLUMNS must track ProcessResourceOverride's fields"
)
