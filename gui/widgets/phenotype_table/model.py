#!/usr/bin/env python3
# model.py — QAbstractTableModel for the Runtime tab's phenotype catalogue.
# PhyloPhere | gui/widgets/phenotype_table/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Backed directly by RuntimeConfig.phenotype_rows: list[PhenotypeRow] — table edits
flow straight into the serializable model, no separate view-only state (see
implementation plan §4). CLASS 1 rows need SECONDARY/NTRAIT/CTRAIT/PRUNE/PRUNE_SEC; CLASS 2 rows need only
TRAIT (TRAIT_TYPE is optional) — irrelevant columns are disabled (not just hidden)
per-row via flags(), and missing-but-required cells are background-tinted so the
mistake is visible without opening a separate validation dialog.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import fields

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QColor

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.runtime import PhenotypeRow

COLUMNS = [
    ("trait_class", "CLASS"),
    ("trait", "TRAIT"),
    ("secondary", "SECONDARY"),
    ("n_trait", "NTRAIT"),
    ("c_trait", "CTRAIT"),
    ("prune", "PRUNE"),
    ("prune_secondary", "PRUNE_SEC"),
    ("trait_type", "TRAIT_TYPE"),
]

_CLASS1_ONLY = {"secondary", "n_trait", "c_trait", "prune", "prune_secondary"}
_CLASS2_ONLY: set[str] = set()  # no CLASS-2-only required field (TRAIT_TYPE is optional)
# CLASS 2 but never required — blank means auto-infer.
_CLASS2_OPTIONAL = {"trait_type"}
_REQUIRED_ALWAYS = {"trait"}

_HEADER_TOOLTIPS = {
    "trait_class": (
        "1 = trait file has n_trait/c_trait columns (sample size + observed cases, e.g. "
        "prevalence) -> Jeffreys CI composition.\n"
        "2 = a single index value per species, no n/c columns -> PSS pair selection "
        "(continuous) or coded fg/bg (ordinal); see TRAIT_TYPE."
    ),
    "prune": (
        "Optional. Filled in => this row is pruned automatically (joined with the Runtime "
        "tab's Prune-list directory). Left blank => this row runs unpruned. No separate "
        "on/off toggle exists — this column is the toggle."
    ),
    "prune_secondary": "Optional secondary prune list, same trigger rule as PRUNE.",
    "trait_type": (
        "CLASS 2, optional. Blank/auto = infer (2-5 integer levels => coded fg/bg).\n"
        "ordinal = force coded: highest level foreground, lowest background, any middle "
        "level intermediate (excluded from contrasts) — use for binary presence/absence "
        "or ordinal category phenotypes.\n"
        "continuous = force the Phylogenetic Shift Score (PSS) pair-selection path."
    ),
}

_MISSING_TINT = QColor(255, 210, 210)
_IRRELEVANT_TINT = QColor(235, 235, 235)


class PhenotypeTableModel(QAbstractTableModel):
    def __init__(self, rows: list[PhenotypeRow], parent=None):
        super().__init__(parent)
        self._rows = rows  # shared reference into RuntimeConfig.phenotype_rows
        self.lang = "en"

    def set_lang(self, lang: str) -> None:
        self.lang = lang
        self.headerDataChanged.emit(Qt.Orientation.Horizontal, 0, len(COLUMNS) - 1)

    # ── Required overrides ─────────────────────────────────────────────────

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            col_key = COLUMNS[section][1]
            from gui.i18n import tr
            return tr(col_key, getattr(self, "lang", "en"))
        if orientation != Qt.Orientation.Horizontal:
            if role == Qt.ItemDataRole.DisplayRole:
                return str(section + 1)
            return None
        if role == Qt.ItemDataRole.ToolTipRole:
            return _HEADER_TOOLTIPS.get(COLUMNS[section][0])
        return None

    def _field_name(self, column: int) -> str:
        return COLUMNS[column][0]

    def _is_irrelevant(self, row: PhenotypeRow, field_name: str) -> bool:
        if row.trait_class == 1:
            return field_name in _CLASS2_ONLY or field_name in _CLASS2_OPTIONAL
        if row.trait_class == 2:
            return field_name in _CLASS1_ONLY
        return False

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row = self._rows[index.row()]
        field_name = self._field_name(index.column())
        value = getattr(row, field_name)

        if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
            return value

        if role == Qt.ItemDataRole.BackgroundRole:
            if self._is_irrelevant(row, field_name):
                return _IRRELEVANT_TINT
            is_required = field_name in _REQUIRED_ALWAYS or (
                field_name in _CLASS2_ONLY and row.trait_class == 2
            )
            if is_required and not str(value).strip():
                return _MISSING_TINT

        return None

    def setData(self, index: QModelIndex, value, role=Qt.ItemDataRole.EditRole) -> bool:
        if role != Qt.ItemDataRole.EditRole or not index.isValid():
            return False
        row = self._rows[index.row()]
        field_name = self._field_name(index.column())
        if field_name == "trait_class":
            try:
                value = int(value)
            except (TypeError, ValueError):
                return False
            if value not in (1, 2):
                return False
        setattr(row, field_name, value)
        self.dataChanged.emit(
            self.index(index.row(), 0),
            self.index(index.row(), self.columnCount() - 1),
        )
        return True

    def flags(self, index: QModelIndex):
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        base = Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled
        row = self._rows[index.row()]
        field_name = self._field_name(index.column())
        if self._is_irrelevant(row, field_name):
            return Qt.ItemFlag.ItemIsSelectable
        return base | Qt.ItemFlag.ItemIsEditable

    # ── Row add/remove ─────────────────────────────────────────────────────

    def insertRows(self, row: int, count: int, parent=QModelIndex()) -> bool:
        self.beginInsertRows(parent, row, row + count - 1)
        for i in range(count):
            self._rows.insert(row + i, PhenotypeRow())
        self.endInsertRows()
        return True

    def removeRows(self, row: int, count: int, parent=QModelIndex()) -> bool:
        if row < 0 or row + count > len(self._rows):
            return False
        self.beginRemoveRows(parent, row, row + count - 1)
        del self._rows[row : row + count]
        self.endRemoveRows()
        return True


assert {name for name, _ in COLUMNS} == {
    f.name for f in fields(PhenotypeRow)
}, "COLUMNS must track PhenotypeRow's fields"
