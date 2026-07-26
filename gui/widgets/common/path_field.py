#!/usr/bin/env python3
# path_field.py — QLineEdit + "Browse..." button compound widget.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QLineEdit, QPushButton, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.widgets.common import remote_context


class PathField(QWidget):
    """A text field paired with a file/directory browse button.

    `mode` is "file" or "dir" — controls which browser is opened. When
    remote_context.get_remote_host() is non-empty at click time, Browse opens an
    SSH-backed RemoteBrowseDialog instead of the local QFileDialog (see
    remote_context.py for why this is read at click time rather than construction
    time).
    """

    textChanged = Signal(str)

    def __init__(self, mode: str = "file", placeholder: str = "", parent=None):
        super().__init__(parent)
        if mode not in ("file", "dir"):
            raise ValueError(f"mode must be 'file' or 'dir', got {mode!r}")
        self._mode = mode

        self.line_edit = QLineEdit(self)
        self.line_edit.setPlaceholderText(placeholder)
        self.line_edit.textChanged.connect(self.textChanged)

        browse_btn = QPushButton("Browse...", self)
        browse_btn.clicked.connect(self._browse)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.line_edit, stretch=1)
        layout.addWidget(browse_btn)

    def _browse(self) -> None:
        host = remote_context.get_remote_host()
        if host:
            self._browse_remote(host)
        else:
            self._browse_local()

    def _browse_local(self) -> None:
        if self._mode == "dir":
            path = QFileDialog.getExistingDirectory(self, "Select directory", self.text())
        else:
            path, _ = QFileDialog.getOpenFileName(self, "Select file", self.text())
        if path:
            self.set_text(path)

    def _browse_remote(self, host: str) -> None:
        from gui.widgets.common.remote_browse_dialog import RemoteBrowseDialog

        current = self.text().strip()
        if current.startswith("/"):
            start_path = current.rsplit("/", 1)[0]
        else:
            # Empty field (or a relative value): start at the configured remote
            # root directory instead of always "/" — General tab's "Remote root
            # directory" lets a user set this once (e.g. "/scratch/mramon")
            # rather than navigating down from "/" on every single field.
            start_path = remote_context.get_remote_root_dir()
        dialog = RemoteBrowseDialog(host, self._mode, start_path or "/", parent=self)
        if dialog.exec() == dialog.DialogCode.Accepted and dialog.selected_path:
            self.set_text(dialog.selected_path)

    def text(self) -> str:
        return self.line_edit.text()

    def set_text(self, value: str) -> None:
        self.line_edit.setText(value)
