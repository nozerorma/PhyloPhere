#!/usr/bin/env python3
# main_window.py — QMainWindow: tab assembly, project New/Open/Save, Generate Scripts.
# PhyloPhere | gui/widgets/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QTabWidget,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import project_io
from gui.generation.render import render_sbatch, render_single
from gui.generation.validate import path_entries, validate, validate_paths
from gui.models.about import AboutInfo
from gui.remote import RemoteCheckError, check_remote_paths
from gui.models.project import ProjectConfig
from gui.widgets.tabs.about_tab import AboutTab
from gui.widgets.tabs.accumulation_tab import AccumulationTab
from gui.widgets.tabs.caas_tab import CaasTab
from gui.widgets.tabs.disambiguation_tab import DisambiguationTab
from gui.widgets.tabs.enrichment_tab import EnrichmentTab
from gui.widgets.tabs.fade_tab import FadeTab
from gui.widgets.tabs.general_tab import GeneralTab
from gui.widgets.tabs.rerconverge_tab import RerConvergeTab
from gui.widgets.tabs.resources_tab import ResourcesTab
from gui.widgets.tabs.runtime_tab import RuntimeTab
from gui.widgets.tabs.scoring_tab import ScoringTab
from gui.widgets.tabs.vep_tab import VepTab


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.project = ProjectConfig()
        self.project_path: Path | None = None
        self._dirty = False

        self.tabs = QTabWidget(self)
        self.setCentralWidget(self.tabs)

        self.about_tab = AboutTab(AboutInfo())  # static, never rebuilt on New/Open
        self._build_project_tabs()
        self.tabs.addTab(self.about_tab, "About")

        self._build_menu()
        self._build_toolbar()
        self._update_title()

    def _build_project_tabs(self) -> None:
        """(Re)builds every tab bound to self.project, in display order, and inserts
        them before the always-last About tab. Called on init and after New/Open."""
        self.general_tab = GeneralTab(self.project.general)
        self.runtime_tab = RuntimeTab(
            self.project.runtime, repo_dir_getter=lambda: self.project.general.repo_dir
        )
        modules = self.project.modules
        self.caas_tab = CaasTab(modules.caas)
        self.disambiguation_tab = DisambiguationTab(modules.disambiguation)
        self.accumulation_tab = AccumulationTab(modules.accumulation)
        self.rerconverge_tab = RerConvergeTab(modules.rer)
        self.fade_tab = FadeTab(modules.fade)
        self.vep_tab = VepTab(modules.vep)
        self.scoring_tab = ScoringTab(modules.scoring)
        self.enrichment_tab = EnrichmentTab(modules.enrichment)
        self.resources_tab = ResourcesTab(self.project.resources)

        self._project_tabs = [
            (self.general_tab, "General"),
            (self.runtime_tab, "Runtime"),
            (self.caas_tab, "CAAS"),
            (self.disambiguation_tab, "Disambiguation"),
            (self.accumulation_tab, "Accumulation"),
            (self.rerconverge_tab, "RERconverge"),
            (self.fade_tab, "FADE"),
            (self.vep_tab, "VEP"),
            (self.scoring_tab, "Scoring"),
            (self.enrichment_tab, "Enrichment"),
            (self.resources_tab, "Resources"),
        ]
        for index, (tab, title) in enumerate(self._project_tabs):
            self.tabs.insertTab(index, tab, title)
            tab.changed.connect(self._mark_dirty)

    # ── Menu ─────────────────────────────────────────────────────────────────

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")

        new_action = QAction("&New Project", self)
        new_action.setShortcut(QKeySequence.StandardKey.New)
        new_action.triggered.connect(self.new_project)
        file_menu.addAction(new_action)

        open_action = QAction("&Open Project...", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self.open_project)
        file_menu.addAction(open_action)

        self.save_action = QAction("&Save Project", self)
        self.save_action.setShortcut(QKeySequence.StandardKey.Save)
        self.save_action.triggered.connect(self.save_project)
        file_menu.addAction(self.save_action)

        save_as_action = QAction("Save Project &As...", self)
        save_as_action.setShortcut(QKeySequence.StandardKey.SaveAs)
        save_as_action.triggered.connect(self.save_project_as)
        file_menu.addAction(save_as_action)

        file_menu.addSeparator()

        generate_action = QAction("&Generate Scripts...", self)
        generate_action.setShortcut("Ctrl+G")
        generate_action.triggered.connect(self.generate_scripts)
        file_menu.addAction(generate_action)

        self.validate_paths_action = QAction("&Validate Paths...", self)
        self.validate_paths_action.setShortcut("Ctrl+Shift+V")
        self.validate_paths_action.triggered.connect(self.validate_paths_dialog)
        file_menu.addAction(self.validate_paths_action)

    def _build_toolbar(self) -> None:
        toolbar = QToolBar("Main", self)
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        toolbar.addAction(self.save_action)
        toolbar.addAction(self.validate_paths_action)

    # ── Dirty tracking ───────────────────────────────────────────────────────

    def _mark_dirty(self) -> None:
        self._dirty = True
        self._update_title()

    def _update_title(self) -> None:
        name = self.project_path.name if self.project_path else "Untitled project"
        star = "*" if self._dirty else ""
        self.setWindowTitle(f"PhyloPhere Runner GUI — {name}{star}")

    def _rebind_tabs_to_project(self) -> None:
        """After New/Open replaces self.project, rebuild every tab bound to it."""
        index = self.tabs.currentIndex()
        for tab, _title in self._project_tabs:
            self.tabs.removeTab(self.tabs.indexOf(tab))
        self._build_project_tabs()
        self.tabs.setCurrentIndex(min(index, self.tabs.count() - 1))

    # ── Project actions ──────────────────────────────────────────────────────

    def new_project(self) -> None:
        self.project = ProjectConfig()
        self.project_path = None
        self._dirty = False
        self._rebind_tabs_to_project()
        self._update_title()

    def open_project(self) -> None:
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open project", "", "PhyloPhere GUI project (*.json)"
        )
        if path_str:
            self.open_project_from_path(Path(path_str))

    def open_project_from_path(self, path: Path) -> None:
        try:
            self.project = project_io.load_project(path)
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Failed to open project", str(exc))
            return
        self.project_path = path
        self._dirty = False
        self._rebind_tabs_to_project()
        self._update_title()

    def save_project(self) -> None:
        if self.project_path is None:
            self.save_project_as()
            return
        project_io.save_project(self.project_path, self.project)
        self._dirty = False
        self._update_title()

    def save_project_as(self) -> None:
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save project as", "", "PhyloPhere GUI project (*.json)"
        )
        if not path_str:
            return
        self.project_path = Path(path_str)
        project_io.save_project(self.project_path, self.project)
        self._dirty = False
        self._update_title()

    def generate_scripts(self) -> None:
        errors = validate(self.project)
        if errors:
            QMessageBox.warning(
                self, "Cannot generate scripts", "\n".join(f"• {e}" for e in errors)
            )
            return

        sbatch_text = render_sbatch(self.project)
        single_text = render_single(self.project)
        self._show_preview(sbatch_text, single_text)

    def validate_paths_dialog(self) -> None:
        """Checks every filled-in path field and reports which ones don't exist.
        Complements Generate Scripts' validate() call, which only checks
        required-ness, not actual existence. Checks the local filesystem, unless
        General > Remote host is set, in which case it checks that host over SSH
        (a single round-trip — see gui/remote.py) instead."""
        host = self.project.general.remote_host.strip()
        target = host or "this filesystem"

        if host:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            try:
                problems = check_remote_paths(host, path_entries(self.project))
            except RemoteCheckError as exc:
                QMessageBox.critical(self, "Remote validation failed", str(exc))
                return
            finally:
                QApplication.restoreOverrideCursor()
        else:
            problems = validate_paths(self.project)

        if not problems:
            QMessageBox.information(self, "Validate paths", f"All filled-in paths exist on {target}.")
            return
        QMessageBox.warning(
            self,
            "Validate paths — missing paths found",
            f"{len(problems)} path(s) not found on {target}:\n\n"
            + "\n".join(f"• {p}" for p in problems),
        )

    def _show_preview(self, sbatch_text: str, single_text: str) -> None:
        preview = QWidget()  # standalone top-level window, not a child of MainWindow
        preview.setWindowTitle("Generated scripts preview")
        layout = QVBoxLayout(preview)

        preview_tabs = QTabWidget(preview)
        for title, text in [
            ("SBATCH_run_phenotypes.sh", sbatch_text),
            ("run_phenotype_single.sh", single_text),
        ]:
            editor = QPlainTextEdit()
            editor.setReadOnly(True)
            editor.setPlainText(text)
            editor.setFont(_monospace_font())
            preview_tabs.addTab(editor, title)
        layout.addWidget(preview_tabs)

        preview.resize(900, 700)
        preview.show()
        self._preview_window = preview  # keep a reference alive


def _monospace_font():
    from PySide6.QtGui import QFontDatabase

    return QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
