#!/usr/bin/env python3
# general_tab.py — General tab: core nextflow.config knobs + remote validation target.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.general import GeneralConfig
from gui.widgets.common import remote_context
from gui.widgets.common.path_field import PathField


class GeneralTab(QWidget):
    changed = Signal()

    def __init__(self, config: GeneralConfig, parent=None):
        super().__init__(parent)
        self._config = config
        remote_context.set_remote_host(config.remote_host)  # sync on tab (re)creation
        remote_context.set_remote_root_dir(config.remote_root_dir)

        layout = QVBoxLayout(self)
        layout.addWidget(self._build_setup_note())
        layout.addWidget(self._build_nextflow_group())
        layout.addWidget(self._build_remote_group())
        layout.addStretch(1)

    def _build_setup_note(self) -> QGroupBox:
        box = QGroupBox("First-time setup")
        layout = QVBoxLayout(box)
        note = QLabel(
            "This GUI runs inside the <code>phylophere</code> conda environment, so it can't "
            "install that environment itself (nothing to launch it with beforehand). If you "
            "haven't set it up yet, run this once from a terminal:<br>"
            "<code>./environment/install_env.sh</code>"
        )
        note.setWordWrap(True)
        note.setTextFormat(note.textFormat().RichText)
        layout.addWidget(note)
        return box

    def _build_nextflow_group(self) -> QGroupBox:
        box = QGroupBox("Core nextflow.config knobs")
        form = QFormLayout(box)

        self.repo_dir = PathField(mode="dir", placeholder="path to PhyloPhere checkout")
        self.repo_dir.set_text(self._config.repo_dir)
        self.repo_dir.textChanged.connect(self._on_repo_dir_changed)
        form.addRow("Repo directory", self.repo_dir)

        self.nextflow_plugins_dir = PathField(mode="dir")
        self.nextflow_plugins_dir.set_text(self._config.nextflow_plugins_dir)
        self.nextflow_plugins_dir.textChanged.connect(self._on_plugins_dir_changed)
        form.addRow("Nextflow plugins directory", self.nextflow_plugins_dir)
        form.addRow(
            "",
            QLabel(
                "Symlinked into every run's NXF_HOME to avoid session/cache-database "
                "lock collisions between concurrent phenotype runs."
            ),
        )

        self.seed = QLineEdit(self._config.seed)
        self.seed.textChanged.connect(self._on_seed_changed)
        form.addRow("Random seed", self.seed)

        self.reporting = QCheckBox("Generate HTML trait-analysis reports (--reporting)")
        self.reporting.setChecked(self._config.reporting)
        self.reporting.toggled.connect(self._on_reporting_changed)
        form.addRow(self.reporting)
        form.addRow(
            "",
            QLabel(
                "Data pruning (--prune_data) is per-phenotype now, not a toggle here — "
                "see the Runtime tab's phenotype catalogue."
            ),
        )

        return box

    def _build_remote_group(self) -> QGroupBox:
        box = QGroupBox("Remote validation && browsing")
        layout = QVBoxLayout(box)

        note = QLabel(
            "Dataset paths often live on a remote HPC cluster rather than this machine. Set a "
            "host here (SSH key-based auth must already work, no password prompt) and every "
            "path field's Browse button (a real file explorer over SSH — lists and navigates "
            "remote directories, not just a validity check) and the Validate Paths action will "
            "check that host instead of the local filesystem."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        form = QFormLayout()
        self.remote_host = QLineEdit(self._config.remote_host)
        self.remote_host.setPlaceholderText("user@cluster.example.edu (blank = local)")
        self.remote_host.textChanged.connect(self._on_remote_host_changed)
        form.addRow("Remote host", self.remote_host)

        self.remote_root_dir = QLineEdit(self._config.remote_root_dir)
        self.remote_root_dir.setPlaceholderText("/scratch/username (blank = start browsing at /)")
        self.remote_root_dir.textChanged.connect(self._on_remote_root_dir_changed)
        form.addRow("Remote root directory", self.remote_root_dir)
        layout.addLayout(form)

        return box

    # ── Slots ────────────────────────────────────────────────────────────────

    def _on_repo_dir_changed(self, value: str) -> None:
        self._config.repo_dir = value
        self.changed.emit()

    def _on_plugins_dir_changed(self, value: str) -> None:
        self._config.nextflow_plugins_dir = value
        self.changed.emit()

    def _on_seed_changed(self, value: str) -> None:
        self._config.seed = value
        self.changed.emit()

    def _on_reporting_changed(self, value: bool) -> None:
        self._config.reporting = value
        self.changed.emit()

    def _on_remote_host_changed(self, value: str) -> None:
        self._config.remote_host = value
        remote_context.set_remote_host(value)
        self.changed.emit()

    def _on_remote_root_dir_changed(self, value: str) -> None:
        self._config.remote_root_dir = value
        remote_context.set_remote_root_dir(value)
        self.changed.emit()
