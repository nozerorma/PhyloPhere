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
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import secrets_io
from gui.models.runtime import RuntimeConfig
from gui.remote import RemoteCheckError
from gui.widgets.common import remote_context
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

    # ── Groups ───────────────────────────────────────────────────────────────

    def _build_top_toggles_group(self) -> QGroupBox:
        self.top_box = QGroupBox("Run behavior")
        form = QFormLayout(self.top_box)

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
        self.toy_n_label = QLabel("Toy sample size (--toy_n)")
        form.addRow(self.toy_n_label, self.toy_n)

        return self.top_box

    def _build_execution_group(self) -> QGroupBox:
        self.exec_box = QGroupBox("Execution target")
        form = QFormLayout(self.exec_box)

        self.runtime_type = QComboBox()
        self.runtime_type.addItems(["slurm", "local"])
        self.runtime_type.setCurrentText(self._config.runtime_type)
        self.runtime_type.currentTextChanged.connect(self._on_runtime_type_changed)
        self.runtime_type_label = QLabel("Runtime")
        form.addRow(self.runtime_type_label, self.runtime_type)

        self.batched = QCheckBox("Batched SLURM array job (generates the SBATCH wrapper)")
        self.batched.setChecked(self._config.batched)
        self.batched.toggled.connect(self._on_batched_changed)
        form.addRow(self.batched)

        self.script_base_name = QLineEdit(self._config.script_base_name)
        self.script_base_name.setPlaceholderText("phenotypes (default)")
        self.script_base_name.textChanged.connect(self._on_script_base_name_changed)
        self.script_base_name_label = QLabel("Generated script name")
        form.addRow(self.script_base_name_label, self.script_base_name)

        self.work_dir = PathField(mode="dir")
        self.work_dir.set_text(self._config.work_dir)
        self.work_dir.textChanged.connect(self._on_work_dir_changed)
        self.work_dir_label = QLabel("Work directory")
        form.addRow(self.work_dir_label, self.work_dir)

        self.results_dir = PathField(mode="dir")
        self.results_dir.set_text(self._config.results_dir)
        self.results_dir.textChanged.connect(self._on_results_dir_changed)
        self.results_dir_label = QLabel("Results directory")
        form.addRow(self.results_dir_label, self.results_dir)

        self.sbatch_job_name = QLineEdit(self._config.sbatch_job_name)
        self.sbatch_job_name.textChanged.connect(self._on_sbatch_job_name_changed)
        self.job_name_label = QLabel("SBATCH job name")
        form.addRow(self.job_name_label, self.sbatch_job_name)

        self.sbatch_partition = QLineEdit(self._config.sbatch_partition)
        self.sbatch_partition.textChanged.connect(self._on_sbatch_partition_changed)
        self.partition_label = QLabel("SBATCH partition")
        form.addRow(self.partition_label, self.sbatch_partition)

        self.sbatch_time = QLineEdit(self._config.sbatch_time)
        self.sbatch_time.textChanged.connect(self._on_sbatch_time_changed)
        self.time_label = QLabel("SBATCH time (-t)")
        form.addRow(self.time_label, self.sbatch_time)

        self.sbatch_mail_user = QLineEdit(self._config.sbatch_mail_user)
        self.sbatch_mail_user.textChanged.connect(self._on_sbatch_mail_user_changed)
        self.mail_user_label = QLabel("SBATCH mail-user (blank = no email)")
        form.addRow(self.mail_user_label, self.sbatch_mail_user)

        self.sbatch_array_concurrency = QLineEdit(self._config.sbatch_array_concurrency)
        self.sbatch_array_concurrency.textChanged.connect(self._on_concurrency_changed)
        self.concurrency_label = QLabel("Array concurrency cap (blank = uncapped)")
        form.addRow(self.concurrency_label, self.sbatch_array_concurrency)

        return self.exec_box

    def _build_tower_group(self) -> QGroupBox:
        self.tower_box = QGroupBox("Seqera Cloud (Tower)")
        form = QFormLayout(self.tower_box)

        self.use_tower = QCheckBox("Enable Seqera Cloud monitoring (-with-tower)")
        self.use_tower.setChecked(self._config.use_tower)
        self.use_tower.toggled.connect(self._on_use_tower_changed)
        form.addRow(self.use_tower)

        token_row = QHBoxLayout()
        self.tower_token = QLineEdit()
        self.tower_token.setEchoMode(QLineEdit.EchoMode.Password)
        self.tower_token.setPlaceholderText("Access token (never saved to the project file)")
        self.tower_token.textChanged.connect(self._on_token_text_changed)
        self.save_token_btn = QPushButton("Save to token.tk")
        self.save_token_btn.clicked.connect(self._save_token)
        token_row.addWidget(self.tower_token, stretch=1)
        token_row.addWidget(self.save_token_btn)
        self.token_label = QLabel("Access token")
        form.addRow(self.token_label, token_row)

        return self.tower_box

    def _build_dataset_group(self) -> QGroupBox:
        self.dataset_box = QGroupBox("Dataset paths (shared across every phenotype in this batch)")
        form = QFormLayout(self.dataset_box)

        self.alignment_dir = PathField(mode="dir")
        self.alignment_dir.set_text(self._config.alignment_dir)
        self.alignment_dir.textChanged.connect(self._on_alignment_dir_changed)
        self.alignment_label = QLabel("Alignments (--alignment)")
        form.addRow(self.alignment_label, self.alignment_dir)

        self.ali_format = QLineEdit(self._config.ali_format)
        self.ali_format.textChanged.connect(self._on_ali_format_changed)
        self.ali_format_label = QLabel("Alignment format (--ali_format)")
        form.addRow(self.ali_format_label, self.ali_format)

        self.tree_file = PathField(mode="file")
        self.tree_file.set_text(self._config.tree_file)
        self.tree_file.textChanged.connect(self._on_tree_file_changed)
        self.tree_file_label = QLabel("Species tree (--tree)")
        form.addRow(self.tree_file_label, self.tree_file)

        self.trait_file = PathField(mode="file")
        self.trait_file.set_text(self._config.trait_file)
        self.trait_file.textChanged.connect(self._on_trait_file_changed)
        self.trait_file_label = QLabel("CLASS 1 trait file")
        form.addRow(self.trait_file_label, self.trait_file)

        self.simple_trait_file = PathField(mode="file")
        self.simple_trait_file.set_text(self._config.simple_trait_file)
        self.simple_trait_file.textChanged.connect(self._on_simple_trait_file_changed)
        self.simple_trait_file_label = QLabel("CLASS 2 trait file")
        form.addRow(self.simple_trait_file_label, self.simple_trait_file)

        self.class_note = QLabel(
            "CLASS 1 = trait file has n_trait/c_trait columns (sample size + observed-case "
            "counts, e.g. disease prevalence) — contrast selection uses a Jeffreys "
            "confidence interval on the proportion, and supports a secondary trait + prune "
            "lists. CLASS 2 = a single index value per species with no n/c columns — "
            "contrast selection picks divergent species pairs via the Phylogenetic Shift "
            "Score (PSS, OU/BM), or coded foreground/background levels when the trait is "
            "ordinal (see TRAIT_TYPE). Each phenotype row below picks its own CLASS."
        )
        self.class_note.setWordWrap(True)
        form.addRow("", self.class_note)

        self.prune_dir = PathField(mode="dir")
        self.prune_dir.set_text(self._config.prune_dir)
        self.prune_dir.textChanged.connect(self._on_prune_dir_changed)
        self.prune_dir_label = QLabel("Prune-list directory")
        form.addRow(self.prune_dir_label, self.prune_dir)
        self.prune_note = QLabel(
            "Pruning is per phenotype, not a separate toggle: a CLASS 1 row with its PRUNE "
            "column filled in (joined with this directory) is pruned automatically; a row "
            "with PRUNE left blank runs unpruned. CLASS 2 rows don't prune."
        )
        self.prune_note.setWordWrap(True)
        form.addRow("", self.prune_note)

        self.branch_trait = QLineEdit(self._config.branch_trait)
        self.branch_trait.textChanged.connect(self._on_branch_trait_changed)
        self.branch_trait_label = QLabel("Branch trait")
        form.addRow(self.branch_trait_label, self.branch_trait)

        self.ali_sp_names = PathField(mode="file")
        self.ali_sp_names.set_text(self._config.ali_sp_names)
        self.ali_sp_names.textChanged.connect(self._on_ali_sp_names_changed)
        self.ali_sp_names_label = QLabel("Alignment species names (optional)")
        form.addRow(self.ali_sp_names_label, self.ali_sp_names)

        self.tax_id_file = PathField(mode="file")
        self.tax_id_file.set_text(self._config.tax_id_file)
        self.tax_id_file.textChanged.connect(self._on_tax_id_file_changed)
        self.tax_id_file_label = QLabel("Taxonomy ID mapping (optional)")
        form.addRow(self.tax_id_file_label, self.tax_id_file)

        self.sp_colname = QLineEdit(self._config.sp_colname)
        self.sp_colname.textChanged.connect(self._on_sp_colname_changed)
        self.sp_colname_label = QLabel("Species column name (--sp_colname)")
        form.addRow(self.sp_colname_label, self.sp_colname)

        self.clade_name = QLineEdit(self._config.clade_name)
        self.clade_name.textChanged.connect(self._on_clade_name_changed)
        self.clade_name_label = QLabel("Clade name (--clade_name)")
        form.addRow(self.clade_name_label, self.clade_name)

        self.taxon_of_interest = QLineEdit(self._config.taxon_of_interest)
        self.taxon_of_interest.textChanged.connect(self._on_taxon_of_interest_changed)
        self.taxon_of_interest_label = QLabel("Taxonomic level of interest (--taxon_of_interest)")
        form.addRow(self.taxon_of_interest_label, self.taxon_of_interest)

        return self.dataset_box

    def _build_catalogue_group(self) -> QGroupBox:
        self.catalogue_box = QGroupBox("Phenotype catalogue")
        layout = QVBoxLayout(self.catalogue_box)
        self.catalogue = PhenotypeTableWidget(self._config.phenotype_rows)
        self.catalogue.changed.connect(self.changed)
        layout.addWidget(self.catalogue)
        return self.catalogue_box

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

    def _on_script_base_name_changed(self, value: str) -> None:
        self._config.script_base_name = value
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
        remote_host = remote_context.get_remote_host()
        try:
            if self._tower_token:
                secrets_io.write_tower_token(repo_dir, self._tower_token, remote_host)
            else:
                secrets_io.clear_tower_token(repo_dir, remote_host)
        except (OSError, RemoteCheckError) as exc:
            target = remote_host or "the local filesystem"
            QMessageBox.critical(
                self,
                "Failed to save token",
                f"Could not write token.tk under repo dir on {target}:\n{repo_dir}\n\n{exc}",
            )

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

    def _on_ali_format_changed(self, value: str) -> None:
        self._config.ali_format = value
        self.changed.emit()

    def _on_sp_colname_changed(self, value: str) -> None:
        self._config.sp_colname = value
        self.changed.emit()

    def _on_clade_name_changed(self, value: str) -> None:
        self._config.clade_name = value
        self.changed.emit()

    def _on_taxon_of_interest_changed(self, value: str) -> None:
        self._config.taxon_of_interest = value
        self.changed.emit()

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "top_box"):
            self.top_box.setTitle(tr("Run behavior", lang))
        if hasattr(self, "exec_box"):
            self.exec_box.setTitle(tr("Execution target", lang))
        if hasattr(self, "tower_box"):
            self.tower_box.setTitle(tr("Seqera Cloud (Tower)", lang))
        if hasattr(self, "dataset_box"):
            self.dataset_box.setTitle(tr("Dataset paths (shared across every phenotype in this batch)", lang))
        if hasattr(self, "toy_n_label"):
            self.toy_n_label.setText(tr("Toy sample size (--toy_n)", lang))
        if hasattr(self, "runtime_type_label"):
            self.runtime_type_label.setText(tr("Runtime", lang))
        if hasattr(self, "work_dir_label"):
            self.work_dir_label.setText(tr("Work directory", lang))
        if hasattr(self, "results_dir_label"):
            self.results_dir_label.setText(tr("Results directory", lang))
        if hasattr(self, "job_name_label"):
            self.job_name_label.setText(tr("SBATCH job name", lang))
        if hasattr(self, "partition_label"):
            self.partition_label.setText(tr("SBATCH partition", lang))
        if hasattr(self, "time_label"):
            self.time_label.setText(tr("SBATCH time (-t)", lang))
        if hasattr(self, "mail_user_label"):
            self.mail_user_label.setText(tr("SBATCH mail-user (blank = no email)", lang))
        if hasattr(self, "concurrency_label"):
            self.concurrency_label.setText(tr("Array concurrency cap (blank = uncapped)", lang))
        if hasattr(self, "use_tower"):
            self.use_tower.setText(tr("Enable Seqera Cloud monitoring (-with-tower)", lang))
        if hasattr(self, "token_label"):
            self.token_label.setText(tr("Access token", lang))
        if hasattr(self, "save_token_btn"):
            self.save_token_btn.setText(tr("Save to token.tk", lang))
        if hasattr(self, "catalogue_box"):
            self.catalogue_box.setTitle(tr("Phenotype catalogue", lang))
        if hasattr(self, "catalogue"):
            self.catalogue.retranslate(lang)
        if hasattr(self, "resume"):
            self.resume.setText(tr("Resume (-resume)", lang))
        if hasattr(self, "toy_mode"):
            self.toy_mode.setText(tr("Toy mode (--toy_mode)", lang))
        if hasattr(self, "batched"):
            self.batched.setText(tr("Batched SLURM array job (generates the SBATCH wrapper)", lang))
        if hasattr(self, "script_base_name_label"):
            self.script_base_name_label.setText(tr("Generated script name", lang))
        if hasattr(self, "alignment_label"):
            self.alignment_label.setText(tr("Alignments (--alignment)", lang))
        if hasattr(self, "ali_format_label"):
            self.ali_format_label.setText(tr("Alignment format (--ali_format)", lang))
        if hasattr(self, "tree_file_label"):
            self.tree_file_label.setText(tr("Species tree (--tree)", lang))
        if hasattr(self, "trait_file_label"):
            self.trait_file_label.setText(tr("CLASS 1 trait file", lang))
        if hasattr(self, "simple_trait_file_label"):
            self.simple_trait_file_label.setText(tr("CLASS 2 trait file", lang))
        if hasattr(self, "class_note"):
            self.class_note.setText(tr(
                "CLASS 1 = trait file has n_trait/c_trait columns (sample size + observed-case "
                "counts, e.g. disease prevalence) — contrast selection uses a Jeffreys "
                "confidence interval on the proportion, and supports a secondary trait + prune "
                "lists. CLASS 2 = a single index value per species with no n/c columns — "
                "contrast selection picks divergent species pairs via the Phylogenetic Shift "
                "Score (PSS, OU/BM), or coded foreground/background levels when the trait is "
                "ordinal (see TRAIT_TYPE). Each phenotype row below picks its own CLASS.",
                lang,
            ))
        if hasattr(self, "prune_dir_label"):
            self.prune_dir_label.setText(tr("Prune-list directory", lang))
        if hasattr(self, "prune_note"):
            self.prune_note.setText(tr(
                "Pruning is per phenotype, not a separate toggle: a CLASS 1 row with its PRUNE "
                "column filled in (joined with this directory) is pruned automatically; a row "
                "with PRUNE left blank runs unpruned. CLASS 2 rows don't prune.",
                lang,
            ))
        if hasattr(self, "branch_trait_label"):
            self.branch_trait_label.setText(tr("Branch trait", lang))
        if hasattr(self, "ali_sp_names_label"):
            self.ali_sp_names_label.setText(tr("Alignment species names (optional)", lang))
        if hasattr(self, "tax_id_file_label"):
            self.tax_id_file_label.setText(tr("Taxonomy ID mapping (optional)", lang))
        if hasattr(self, "sp_colname_label"):
            self.sp_colname_label.setText(tr("Species column name (--sp_colname)", lang))
        if hasattr(self, "clade_name_label"):
            self.clade_name_label.setText(tr("Clade name (--clade_name)", lang))
        if hasattr(self, "taxon_of_interest_label"):
            self.taxon_of_interest_label.setText(tr("Taxonomic level of interest (--taxon_of_interest)", lang))
