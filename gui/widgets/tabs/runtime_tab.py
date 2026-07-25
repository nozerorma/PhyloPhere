#!/usr/bin/env python3
# runtime_tab.py — Runtime tab: local/slurm, phenotype catalogue, dataset paths, Tower.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
The Seqera/Tower access token is deliberately NOT bound to RuntimeConfig (see
gui/models/runtime.py's module docstring) — it's transient widget state, written
straight to the repo's gitignored token.tk via gui/secrets_io.py on demand rather
than persisted into the JSON project file.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import secrets_io
from gui.models.runtime import RuntimeConfig
from gui.widgets.common.path_field import PathField
from gui.widgets.phenotype_table.widget import PhenotypeTableWidget


class RuntimeTab(QWidget):
    changed = Signal()

    def __init__(self, config: RuntimeConfig, repo_dir_getter, parent=None):
        """`repo_dir_getter` is a zero-arg callable returning the current General-tab
        repo_dir, so the Tower token can be written to the right token.tk without
        RuntimeTab needing a live reference to GeneralConfig."""
        super().__init__(parent)
        self._config = config
        self._repo_dir_getter = repo_dir_getter
        self._tower_token = ""  # transient, never persisted

        layout = QVBoxLayout(self)
        layout.addWidget(self._build_top_toggles_group())
        layout.addWidget(self._build_execution_group())
        layout.addWidget(self._build_tower_group())
        layout.addWidget(self._build_dataset_group())
        layout.addWidget(self._build_catalogue_group(), stretch=1)

    # ── Groups ───────────────────────────────────────────────────────────────

    def _build_top_toggles_group(self) -> QGroupBox:
        box = QGroupBox("Run behavior")
        form = QFormLayout(box)

        self.resume = QCheckBox("Resume (-resume)")
        self.resume.setChecked(self._config.resume)
        self.resume.toggled.connect(self._on_resume_changed)
        form.addRow(self.resume)

        self.toy_mode = QCheckBox("Toy mode (--toy_mode)")
        self.toy_mode.setChecked(self._config.toy_mode)
        self.toy_mode.toggled.connect(self._on_toy_mode_changed)
        form.addRow(self.toy_mode)

        self.toy_n = QLineEdit(self._config.toy_n)
        self.toy_n.textChanged.connect(self._on_toy_n_changed)
        form.addRow("Toy sample size (--toy_n)", self.toy_n)

        return box

    def _build_execution_group(self) -> QGroupBox:
        box = QGroupBox("Execution target")
        form = QFormLayout(box)

        self.runtime_type = QComboBox()
        self.runtime_type.addItems(["slurm", "local"])
        self.runtime_type.setCurrentText(self._config.runtime_type)
        self.runtime_type.currentTextChanged.connect(self._on_runtime_type_changed)
        form.addRow("Runtime", self.runtime_type)

        self.batched = QCheckBox("Batched SLURM array job (generates the SBATCH wrapper)")
        self.batched.setChecked(self._config.batched)
        self.batched.toggled.connect(self._on_batched_changed)
        form.addRow(self.batched)

        self.work_dir = PathField(mode="dir")
        self.work_dir.set_text(self._config.work_dir)
        self.work_dir.textChanged.connect(self._on_work_dir_changed)
        form.addRow("Work directory", self.work_dir)

        self.results_dir = PathField(mode="dir")
        self.results_dir.set_text(self._config.results_dir)
        self.results_dir.textChanged.connect(self._on_results_dir_changed)
        form.addRow("Results directory", self.results_dir)

        self.sbatch_job_name = QLineEdit(self._config.sbatch_job_name)
        self.sbatch_job_name.textChanged.connect(self._on_sbatch_job_name_changed)
        form.addRow("SBATCH job name", self.sbatch_job_name)

        self.sbatch_partition = QLineEdit(self._config.sbatch_partition)
        self.sbatch_partition.textChanged.connect(self._on_sbatch_partition_changed)
        form.addRow("SBATCH partition", self.sbatch_partition)

        self.sbatch_time = QLineEdit(self._config.sbatch_time)
        self.sbatch_time.textChanged.connect(self._on_sbatch_time_changed)
        form.addRow("SBATCH time (-t)", self.sbatch_time)

        self.sbatch_mail_user = QLineEdit(self._config.sbatch_mail_user)
        self.sbatch_mail_user.textChanged.connect(self._on_sbatch_mail_user_changed)
        form.addRow("SBATCH mail-user (blank = no email)", self.sbatch_mail_user)

        self.sbatch_array_concurrency = QLineEdit(self._config.sbatch_array_concurrency)
        self.sbatch_array_concurrency.textChanged.connect(self._on_concurrency_changed)
        form.addRow("Array concurrency cap (blank = uncapped)", self.sbatch_array_concurrency)

        return box

    def _build_tower_group(self) -> QGroupBox:
        box = QGroupBox("Seqera Cloud (Tower)")
        form = QFormLayout(box)

        self.use_tower = QCheckBox("Enable Seqera Cloud monitoring (-with-tower)")
        self.use_tower.setChecked(self._config.use_tower)
        self.use_tower.toggled.connect(self._on_use_tower_changed)
        form.addRow(self.use_tower)

        token_row = QHBoxLayout()
        self.tower_token = QLineEdit()
        self.tower_token.setEchoMode(QLineEdit.EchoMode.Password)
        self.tower_token.setPlaceholderText("Access token (never saved to the project file)")
        self.tower_token.textChanged.connect(self._on_token_text_changed)
        save_token_btn = QPushButton("Save to token.tk")
        save_token_btn.clicked.connect(self._save_token)
        token_row.addWidget(self.tower_token, stretch=1)
        token_row.addWidget(save_token_btn)
        form.addRow("Access token", token_row)

        return box

    def _build_dataset_group(self) -> QGroupBox:
        box = QGroupBox("Dataset paths (shared across every phenotype in this batch)")
        form = QFormLayout(box)

        self.alignment_dir = PathField(mode="dir")
        self.alignment_dir.set_text(self._config.alignment_dir)
        self.alignment_dir.textChanged.connect(self._on_alignment_dir_changed)
        form.addRow("Alignments (--alignment)", self.alignment_dir)

        self.tree_file = PathField(mode="file")
        self.tree_file.set_text(self._config.tree_file)
        self.tree_file.textChanged.connect(self._on_tree_file_changed)
        form.addRow("Species tree (--tree)", self.tree_file)

        self.trait_file = PathField(mode="file")
        self.trait_file.set_text(self._config.trait_file)
        self.trait_file.textChanged.connect(self._on_trait_file_changed)
        form.addRow("CLASS 1 trait file", self.trait_file)

        self.simple_trait_file = PathField(mode="file")
        self.simple_trait_file.set_text(self._config.simple_trait_file)
        self.simple_trait_file.textChanged.connect(self._on_simple_trait_file_changed)
        form.addRow("CLASS 2 trait file", self.simple_trait_file)

        self.prune_dir = PathField(mode="dir")
        self.prune_dir.set_text(self._config.prune_dir)
        self.prune_dir.textChanged.connect(self._on_prune_dir_changed)
        form.addRow("Prune-list directory", self.prune_dir)

        self.branch_trait = QLineEdit(self._config.branch_trait)
        self.branch_trait.textChanged.connect(self._on_branch_trait_changed)
        form.addRow("Branch trait", self.branch_trait)

        self.ali_sp_names = PathField(mode="file")
        self.ali_sp_names.set_text(self._config.ali_sp_names)
        self.ali_sp_names.textChanged.connect(self._on_ali_sp_names_changed)
        form.addRow("Alignment species names (optional)", self.ali_sp_names)

        self.tax_id_file = PathField(mode="file")
        self.tax_id_file.set_text(self._config.tax_id_file)
        self.tax_id_file.textChanged.connect(self._on_tax_id_file_changed)
        form.addRow("Taxonomy ID mapping (optional)", self.tax_id_file)

        return box

    def _build_catalogue_group(self) -> QGroupBox:
        box = QGroupBox("Phenotype catalogue")
        layout = QVBoxLayout(box)
        self.catalogue = PhenotypeTableWidget(self._config.phenotype_rows)
        self.catalogue.changed.connect(self.changed)
        layout.addWidget(self.catalogue)
        return box

    # ── Slots ────────────────────────────────────────────────────────────────

    def _on_resume_changed(self, value: bool) -> None:
        self._config.resume = value
        self.changed.emit()

    def _on_toy_mode_changed(self, value: bool) -> None:
        self._config.toy_mode = value
        self.changed.emit()

    def _on_toy_n_changed(self, value: str) -> None:
        self._config.toy_n = value
        self.changed.emit()

    def _on_runtime_type_changed(self, value: str) -> None:
        self._config.runtime_type = value
        self.changed.emit()

    def _on_batched_changed(self, value: bool) -> None:
        self._config.batched = value
        self.changed.emit()

    def _on_work_dir_changed(self, value: str) -> None:
        self._config.work_dir = value
        self.changed.emit()

    def _on_results_dir_changed(self, value: str) -> None:
        self._config.results_dir = value
        self.changed.emit()

    def _on_sbatch_job_name_changed(self, value: str) -> None:
        self._config.sbatch_job_name = value
        self.changed.emit()

    def _on_sbatch_partition_changed(self, value: str) -> None:
        self._config.sbatch_partition = value
        self.changed.emit()

    def _on_sbatch_time_changed(self, value: str) -> None:
        self._config.sbatch_time = value
        self.changed.emit()

    def _on_sbatch_mail_user_changed(self, value: str) -> None:
        self._config.sbatch_mail_user = value
        self.changed.emit()

    def _on_concurrency_changed(self, value: str) -> None:
        self._config.sbatch_array_concurrency = value
        self.changed.emit()

    def _on_use_tower_changed(self, value: bool) -> None:
        self._config.use_tower = value
        self.changed.emit()

    def _on_token_text_changed(self, value: str) -> None:
        self._tower_token = value  # transient only, not written to _config

    def _save_token(self) -> None:
        repo_dir = self._repo_dir_getter()
        if not repo_dir:
            return
        if self._tower_token:
            secrets_io.write_tower_token(repo_dir, self._tower_token)
        else:
            secrets_io.clear_tower_token(repo_dir)

    def _on_alignment_dir_changed(self, value: str) -> None:
        self._config.alignment_dir = value
        self.changed.emit()

    def _on_tree_file_changed(self, value: str) -> None:
        self._config.tree_file = value
        self.changed.emit()

    def _on_trait_file_changed(self, value: str) -> None:
        self._config.trait_file = value
        self.changed.emit()

    def _on_simple_trait_file_changed(self, value: str) -> None:
        self._config.simple_trait_file = value
        self.changed.emit()

    def _on_prune_dir_changed(self, value: str) -> None:
        self._config.prune_dir = value
        self.changed.emit()

    def _on_branch_trait_changed(self, value: str) -> None:
        self._config.branch_trait = value
        self.changed.emit()

    def _on_ali_sp_names_changed(self, value: str) -> None:
        self._config.ali_sp_names = value
        self.changed.emit()

    def _on_tax_id_file_changed(self, value: str) -> None:
        self._config.tax_id_file = value
        self.changed.emit()
