#!/usr/bin/env python3
# remote_browse_dialog.py — Lightweight SSH-backed directory browser for PathField.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Stands in for QFileDialog when a remote host is configured (see remote_context.py)
— QFileDialog can only browse the local filesystem. Deliberately minimal: a path
bar, a listing, and double-click navigation. No SFTP/paramiko dependency — reuses
gui.remote's plain `ssh` + `find` based listing.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.remote import RemoteCheckError, list_remote_directory


class RemoteBrowseDialog(QDialog):
    """mode="file": double-click a file (or select + OK) to pick it.
    mode="dir": navigate freely; "Select this directory" picks the current path."""

    def __init__(self, host: str, mode: str, start_path: str = "/", parent=None):
        super().__init__(parent)
        self._host = host
        self._mode = mode
        self.selected_path: str = ""

        self.setWindowTitle(f"Browse {host}")
        self.resize(600, 450)

        self.path_edit = QLineEdit(start_path or "/")
        self.path_edit.returnPressed.connect(self._navigate_to_typed_path)
        go_btn = QPushButton("Go")
        go_btn.clicked.connect(self._navigate_to_typed_path)
        up_btn = QPushButton("Up")
        up_btn.clicked.connect(self._go_up)

        path_row = QHBoxLayout()
        path_row.addWidget(self.path_edit, stretch=1)
        path_row.addWidget(go_btn)
        path_row.addWidget(up_btn)

        self.listing = QListWidget()
        self.listing.itemDoubleClicked.connect(self._on_item_double_clicked)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
        )
        if mode == "dir":
            select_btn = QPushButton("Select this directory")
            select_btn.clicked.connect(self._select_current_directory)
            buttons.addButton(select_btn, QDialogButtonBox.ButtonRole.AcceptRole)
        else:
            select_btn = QPushButton("Select highlighted file")
            select_btn.clicked.connect(self._select_highlighted_file)
            buttons.addButton(select_btn, QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(path_row)
        layout.addWidget(self.listing, stretch=1)
        layout.addWidget(QLabel("Double-click a directory to navigate into it."))
        layout.addWidget(buttons)

        self._refresh()

    def _current_path(self) -> str:
        return self.path_edit.text().strip() or "/"

    def _refresh(self) -> None:
        path = self._current_path()
        try:
            entries = list_remote_directory(self._host, path)
        except RemoteCheckError as exc:
            QMessageBox.warning(self, "Remote browse failed", str(exc))
            return
        self.listing.clear()
        for name, is_dir in entries:
            item = QListWidgetItem(f"{name}/" if is_dir else name)
            item.setData(Qt.ItemDataRole.UserRole, (name, is_dir))
            self.listing.addItem(item)

    def _navigate_to_typed_path(self) -> None:
        self._refresh()

    def _go_up(self) -> None:
        path = self._current_path().rstrip("/")
        parent = path.rsplit("/", 1)[0] or "/"
        self.path_edit.setText(parent)
        self._refresh()

    def _on_item_double_clicked(self, item: QListWidgetItem) -> None:
        name, is_dir = item.data(Qt.ItemDataRole.UserRole)
        base = self._current_path().rstrip("/")
        new_path = f"{base}/{name}" if base else f"/{name}"
        if is_dir:
            self.path_edit.setText(new_path)
            self._refresh()
        elif self._mode == "file":
            self.selected_path = new_path
            self.accept()

    def _select_current_directory(self) -> None:
        self.selected_path = self._current_path()
        self.accept()

    def _select_highlighted_file(self) -> None:
        item = self.listing.currentItem()
        if item is None:
            return
        name, is_dir = item.data(Qt.ItemDataRole.UserRole)
        if is_dir:
            return
        base = self._current_path().rstrip("/")
        self.selected_path = f"{base}/{name}" if base else f"/{name}"
        self.accept()
