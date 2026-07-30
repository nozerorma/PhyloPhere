#!/usr/bin/env python3
# general_tab.py — General tab: core nextflow.config knobs + remote validation target.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.i18n import LANGUAGES, tr
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
        self.setup_box = QGroupBox("First-time setup")
        layout = QVBoxLayout(self.setup_box)
        self.setup_note = QLabel(
            "This GUI runs inside the <code>phylophere</code> conda environment, so it can't "
            "install that environment itself (nothing to launch it with beforehand). If you "
            "haven't set it up yet, run this once from a terminal:<br>"
            "<code>./environment/install_env.sh</code>"
        )
        self.setup_note.setWordWrap(True)
        self.setup_note.setTextFormat(self.setup_note.textFormat().RichText)
        layout.addWidget(self.setup_note)
        return self.setup_box

    def _build_nextflow_group(self) -> QGroupBox:
        self.nextflow_box = QGroupBox("Core nextflow.config knobs")
        form = QFormLayout(self.nextflow_box)

        self.repo_dir = PathField(mode="dir", placeholder="path to PhyloPhere checkout")
        self.repo_dir.set_text(self._config.repo_dir)
        self.repo_dir.textChanged.connect(self._on_repo_dir_changed)
        self.repo_dir_label = QLabel("Repo directory")
        form.addRow(self.repo_dir_label, self.repo_dir)

        self.nextflow_plugins_dir = PathField(mode="dir")
        self.nextflow_plugins_dir.set_text(self._config.nextflow_plugins_dir)
        self.nextflow_plugins_dir.textChanged.connect(self._on_plugins_dir_changed)
        self.plugins_dir_label = QLabel("Nextflow plugins directory")
        form.addRow(self.plugins_dir_label, self.nextflow_plugins_dir)
        self.plugins_dir_note = QLabel(
            "Symlinked into every run's NXF_HOME to avoid session/cache-database "
            "lock collisions between concurrent phenotype runs."
        )
        form.addRow("", self.plugins_dir_note)

        self.seed = QLineEdit(self._config.seed)
        self.seed.textChanged.connect(self._on_seed_changed)
        self.seed_label = QLabel("Random seed")
        form.addRow(self.seed_label, self.seed)

        self.reporting = QCheckBox("Generate HTML trait-analysis reports (--reporting)")
        self.reporting.setChecked(self._config.reporting)
        self.reporting.toggled.connect(self._on_reporting_changed)
        form.addRow(self.reporting)

        return self.nextflow_box

    def _build_remote_group(self) -> QGroupBox:
        self.remote_box = QGroupBox("Remote validation && browsing")
        layout = QVBoxLayout(self.remote_box)

        self.remote_note = QLabel(
            "Dataset paths often live on a remote HPC cluster rather than this machine. Set a "
            "host here (SSH key-based auth must already work, no password prompt) and every "
            "path field's Browse button (a real file explorer over SSH — lists and navigates "
            "remote directories, not just a validity check) and the Validate Paths action will "
            "check that host instead of the local filesystem."
        )
        self.remote_note.setWordWrap(True)
        layout.addWidget(self.remote_note)

        form = QFormLayout()
        self.remote_host = QLineEdit(self._config.remote_host)
        self.remote_host.setPlaceholderText("user@cluster.example.edu (blank = local)")
        self.remote_host.textChanged.connect(self._on_remote_host_changed)
        self.remote_host_label = QLabel("Remote host")
        form.addRow(self.remote_host_label, self.remote_host)

        self.remote_root_dir = QLineEdit(self._config.remote_root_dir)
        self.remote_root_dir.setPlaceholderText("/scratch/username (blank = start browsing at /)")
        self.remote_root_dir.textChanged.connect(self._on_remote_root_dir_changed)
        self.remote_root_dir_label = QLabel("Remote root directory")
        form.addRow(self.remote_root_dir_label, self.remote_root_dir)
        layout.addLayout(form)

        return self.remote_box

    def retranslate(self, lang: str = "en") -> None:
        self.setup_box.setTitle(tr("First-time setup", lang))
        self.nextflow_box.setTitle(tr("Core nextflow.config knobs", lang))
        self.remote_box.setTitle(tr("Remote validation && browsing", lang))
        self.repo_dir_label.setText(tr("Repo directory", lang))
        self.plugins_dir_label.setText(tr("Nextflow plugins directory", lang))
        self.seed_label.setText(tr("Random seed", lang))
        self.reporting.setText(tr("Generate HTML trait-analysis reports (--reporting)", lang))
        self.remote_host_label.setText(tr("Remote host", lang))
        self.remote_root_dir_label.setText(tr("Remote root directory", lang))
        if hasattr(self, "setup_note"):
            self.setup_note.setText(tr(
                "This GUI runs inside the <code>phylophere</code> conda environment, so it can't "
                "install that environment itself (nothing to launch it with beforehand). If you "
                "haven't set it up yet, run this once from a terminal:<br>"
                "<code>./environment/install_env.sh</code>",
                lang,
            ))
        if hasattr(self, "plugins_dir_note"):
            self.plugins_dir_note.setText(tr(
                "Symlinked into every run's NXF_HOME to avoid session/cache-database "
                "lock collisions between concurrent phenotype runs.",
                lang,
            ))
        if hasattr(self, "remote_note"):
            self.remote_note.setText(tr(
                "Dataset paths often live on a remote HPC cluster rather than this machine. Set a "
                "host here (SSH key-based auth must already work, no password prompt) and every "
                "path field's Browse button (a real file explorer over SSH — lists and navigates "
                "remote directories, not just a validity check) and the Validate Paths action will "
                "check that host instead of the local filesystem.",
                lang,
            ))

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
