#!/usr/bin/env python3
# main_window.py — QMainWindow: tab assembly, project New/Open/Save, Generate Scripts.
# PhyloPhere | gui/widgets/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import autosave_io, project_io, remote
from gui.generation.render import render_batch, render_single
from gui.generation.validate import path_entries, validate, validate_paths
from gui.i18n import LANGUAGES, tr
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
from gui.widgets.tabs.precomputed_tab import PrecomputedTab
from gui.widgets.tabs.rerconverge_tab import RerConvergeTab
from gui.widgets.tabs.resources_tab import ResourcesTab
from gui.widgets.tabs.runtime_tab import RuntimeTab
from gui.widgets.tabs.scoring_tab import ScoringTab
from gui.widgets.tabs.vep_tab import VepTab


TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "templates"

_SETTINGS_LANGUAGE_KEY = "gui/language"  # QSettings key — persists across restarts (see app.py's
                                          # setOrganizationName/setApplicationName for where this is stored)


def _scrollable(widget: QWidget) -> QScrollArea:
    """Wrap a tab's content in a QScrollArea — several module tabs now carry 20+
    fields (see the ct.config/enrichment.config coverage pass) and shouldn't force
    the main window to grow to fit the tallest one."""
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setWidget(widget)
    return area


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        restored = autosave_io.read_autosave()
        if restored is not None:
            self.project, self.project_path, self._dirty = restored
        else:
            self.project = ProjectConfig()
            self.project_path = None
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
        self.precomputed_tab = PrecomputedTab(
            self.project.precomputed,
            caas_tab=self.caas_tab,
            disambiguation_tab=self.disambiguation_tab,
            accumulation_tab=self.accumulation_tab,
            rer_tab=self.rerconverge_tab,
            fade_tab=self.fade_tab,
            vep_tab=self.vep_tab,
        )
        self.resources_tab = ResourcesTab(self.project.resources)

        # (inner widget, title) — the inner widget is what emits `changed` and what
        # holds the actual model reference; each gets wrapped in a QScrollArea only
        # for insertion into the QTabWidget (several module tabs now carry 20+
        # fields and shouldn't force the window to grow to fit the tallest one).
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
            (self.precomputed_tab, "Precomputed Run"),
            (self.resources_tab, "Resources"),
        ]
        self._project_tab_wrappers = []
        for index, (tab, title) in enumerate(self._project_tabs):
            wrapper = _scrollable(tab)
            self._project_tab_wrappers.append(wrapper)
            self.tabs.insertTab(index, wrapper, title)
            tab.changed.connect(self._mark_dirty)

    # ── Menu ─────────────────────────────────────────────────────────────────

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")

        self.new_project_action = QAction("&New Project", self)
        self.new_project_action.setShortcut(QKeySequence.StandardKey.New)
        self.new_project_action.triggered.connect(self.new_project)
        file_menu.addAction(self.new_project_action)

        open_action = QAction("&Open Project...", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self.open_project)
        file_menu.addAction(open_action)

        file_menu.addSeparator()

        self.save_template_action = QAction("&Save Project As...", self)
        self.save_template_action.setShortcut(QKeySequence.StandardKey.Save)
        self.save_template_action.triggered.connect(self.save_template)
        file_menu.addAction(self.save_template_action)

        self.load_template_action = QAction("&Load Template...", self)
        self.load_template_action.triggered.connect(self.load_template)
        file_menu.addAction(self.load_template_action)

        file_menu.addSeparator()

        self.validate_paths_action = QAction("&Validate Paths...", self)
        self.validate_paths_action.setShortcut("Ctrl+Shift+V")
        self.validate_paths_action.triggered.connect(self.validate_paths_dialog)
        file_menu.addAction(self.validate_paths_action)

        self.generate_action = QAction("&Generate Scripts...", self)
        self.generate_action.setShortcut("Ctrl+G")
        self.generate_action.triggered.connect(self.generate_scripts)
        file_menu.addAction(self.generate_action)

        file_menu.addSeparator()

        self.readme_action = QAction("&Readme", self)
        self.readme_action.triggered.connect(self.show_readme_dialog)
        file_menu.addAction(self.readme_action)

    def _build_toolbar(self) -> None:
        toolbar = QToolBar("Main", self)
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        # Language selection combo box to the LEFT of Save Project As + Separator
        self.lang_combo = QComboBox(self)
        for code, label in LANGUAGES.items():
            self.lang_combo.addItem(label, userData=code)

        # Restore the language picked in a previous run (falls back to "en" for a
        # first-ever launch or an unrecognized/stale saved code).
        saved_lang = QSettings().value(_SETTINGS_LANGUAGE_KEY, "en", type=str)
        restored_index = self.lang_combo.findData(saved_lang)
        current_lang = saved_lang if restored_index >= 0 else "en"
        self.lang_combo.setCurrentIndex(max(restored_index, 0))

        self.lang_combo.currentIndexChanged.connect(self._on_toolbar_language_changed)

        toolbar.addWidget(self.lang_combo)
        toolbar.addSeparator()

        toolbar.addAction(self.new_project_action)
        toolbar.addSeparator()

        toolbar.addAction(self.save_template_action)
        toolbar.addAction(self.load_template_action)
        toolbar.addSeparator()
        toolbar.addAction(self.validate_paths_action)
        toolbar.addAction(self.generate_action)
        toolbar.addSeparator()
        toolbar.addAction(self.readme_action)

        self._apply_language(current_lang)

    def _on_toolbar_language_changed(self, index: int) -> None:
        code = self.lang_combo.itemData(index)
        if code:
            QSettings().setValue(_SETTINGS_LANGUAGE_KEY, code)
            self._apply_language(code)

    def _apply_language(self, lang: str = "en") -> None:
        self._current_lang = lang
        for idx, (tab, orig_title) in enumerate(self._project_tabs):
            translated_title = tr(f"Tab: {orig_title}", lang)
            if translated_title == f"Tab: {orig_title}":
                translated_title = tr(orig_title, lang)
            self.tabs.setTabText(idx, translated_title)

        self.save_template_action.setText(tr("Save Project As...", lang))
        self.load_template_action.setText(tr("Load Template...", lang))
        self.validate_paths_action.setText(tr("Validate Paths...", lang))
        self.generate_action.setText(tr("Generate Scripts...", lang))
        self.readme_action.setText(tr("Readme", lang))

        if hasattr(self, "about_tab") and hasattr(self.about_tab, "retranslate"):
            self.about_tab.retranslate(lang)
            about_idx = self.tabs.count() - 1
            if about_idx >= 0:
                self.tabs.setTabText(about_idx, tr("Tab: About", lang))

        for tab, _ in self._project_tabs:
            if hasattr(tab, "retranslate"):
                tab.retranslate(lang)

    def show_readme_dialog(self) -> None:
        from gui.widgets.common.readme_dialog import ReadmeDialog
        lang = getattr(self, "_current_lang", "en")
        dialog = ReadmeDialog(self, repo_dir=self.project.general.repo_dir, lang=lang)
        dialog.exec()


    # ── Dirty tracking ───────────────────────────────────────────────────────

    def _mark_dirty(self) -> None:
        self._dirty = True
        self._update_title()
        self._write_autosave()

    def _write_autosave(self) -> None:
        """Best-effort snapshot of the live project state, restored on next
        launch by __init__ regardless of whether the user ever hit Save (see
        gui/autosave_io.py) — called after every mutation to self.project,
        self.project_path, or self._dirty."""
        autosave_io.write_autosave(self.project, self.project_path, self._dirty)

    def _update_title(self) -> None:
        name = (
            self.project.general.project_name.strip()
            or (self.project_path.name if self.project_path else "")
            or "Untitled project"
        )
        star = "*" if self._dirty else ""
        self.setWindowTitle(f"PhyloPhere Runner GUI — {name}{star}")

    def _rebind_tabs_to_project(self) -> None:
        """After New/Open replaces self.project, rebuild every tab bound to it."""
        index = self.tabs.currentIndex()
        for wrapper in self._project_tab_wrappers:
            self.tabs.removeTab(self.tabs.indexOf(wrapper))
        self._build_project_tabs()
        current_lang = "en"
        if hasattr(self, "lang_combo"):
            current_lang = self.lang_combo.currentData() or "en"
        self._apply_language(current_lang)
        self.tabs.setCurrentIndex(min(index, self.tabs.count() - 1))

    # ── Project actions ──────────────────────────────────────────────────────

    def new_project(self) -> None:
        self.project = ProjectConfig()
        self.project_path = None
        self._dirty = False
        self._rebind_tabs_to_project()
        self._update_title()
        self._write_autosave()

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
        self._write_autosave()

    def save_template(self) -> None:
        """Always prompts for a filename — a real "Save As", not a silent in-place
        overwrite (see the toolbar's New Project + this dialog's own prior
        confirmation-on-overwrite for the incident that motivated this: loading
        "primates_full_run", tweaking it for a one-off trait, and saving clobbered
        the original template with no warning). Pre-fills the suggested filename
        from General > Project name when set, so naming the project once feeds
        straight into where it gets saved. The native save dialog's own "file
        exists, overwrite?" prompt covers the case where you deliberately pick an
        existing path again — no separate custom confirmation needed."""
        start_dir = str(self.project_path.parent) if self.project_path else str(TEMPLATES_DIR)
        TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)
        suggested_name = self.project.general.project_name.strip()
        if not suggested_name:
            suggested_name = self.project_path.stem if self.project_path else "untitled"
        suggested_name = "".join(
            c if (c.isalnum() or c in "-_") else "_" for c in suggested_name
        )
        start_path = str(Path(start_dir) / f"{suggested_name}.json")

        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save Project As", start_path, "PhyloPhere GUI project (*.json)"
        )
        if not path_str:
            return
        self.project_path = Path(path_str)
        project_io.save_project(self.project_path, self.project)
        self._dirty = False
        self._update_title()
        self._write_autosave()

    def load_template(self) -> None:
        """Loads one of the bundled example projects under gui/templates/ (see
        save_template) — same underlying JSON format as a regular project file,
        just browsed from a fixed shared directory instead of wherever the user
        last saved to."""
        TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load template", str(TEMPLATES_DIR), "PhyloPhere GUI project (*.json)"
        )
        if path_str:
            self.open_project_from_path(Path(path_str))

    def generate_scripts(self) -> None:
        errors = validate(self.project)
        if errors:
            QMessageBox.warning(
                self, "Cannot generate scripts", "\n".join(f"• {e}" for e in errors)
            )
            return

        is_local = self.project.runtime.runtime_type == "local"
        disambig = self.project.modules.disambiguation
        base = self.project.runtime.script_base_name.strip()

        def script_names(tag: str) -> tuple[str, str]:
            """(batch_name, single_name) for one postproc-mode tag ("", "_exploratory",
            "_complete") — base_name blank reproduces the exact default filenames;
            set, it overrides the "phenotypes"/"phenotype" token everywhere so both
            names change together (they must, see render_batch's single_runner_filename)."""
            if not base:
                plural = "phenotypes"
                singular = "phenotype"
            else:
                plural = singular = base
            batch_name = f"run_{plural}_local{tag}.sh" if is_local else f"SBATCH_run_{plural}{tag}.sh"
            single_name = f"run_{singular}_single{tag}.sh"
            return batch_name, single_name

        scripts: list[tuple[str, str]] = []

        run_expl = getattr(disambig, "run_postproc_exploratory", False)
        run_filt = getattr(disambig, "run_postproc_filter", False)

        if disambig.enabled and (run_expl or run_filt):
            if run_expl:
                batch_name, single_name = script_names("_exploratory")
                batch_text = render_batch(
                    self.project, postproc_mode="exploratory", single_runner_filename=single_name
                )
                single_text = render_single(self.project, postproc_mode="exploratory")
                scripts.append((batch_name, batch_text))
                scripts.append((single_name, single_text))

            if run_filt:
                # reuse_exploratory=True (both selected) skips re-running CAAS/Disambiguation
                # live in the _complete pass — it reads their output straight from the
                # sibling _exploratory run's outdir instead (see run_single.sh.j2's
                # POSTPROC_MODE="filter" branch). Standalone (run_expl off), everything
                # in _complete stays live, same as before.
                reuse_exploratory = run_expl
                batch_name, single_name = script_names("_complete")
                batch_text = render_batch(
                    self.project,
                    postproc_mode="filter",
                    reuse_exploratory=reuse_exploratory,
                    single_runner_filename=single_name,
                )
                single_text = render_single(
                    self.project, postproc_mode="filter", reuse_exploratory=reuse_exploratory
                )
                scripts.append((batch_name, batch_text))
                scripts.append((single_name, single_text))

            if run_expl and run_filt:
                expl_name, _ = script_names("_exploratory")
                complete_name, _ = script_names("_complete")
                QMessageBox.information(
                    self,
                    "Two-step Post-processing workflow",
                    "Both Exploratory Sweep and Complete were selected, so two independent "
                    f"sets of scripts were generated.\n\n"
                    f"1. Run {expl_name} first and check its report — it sweeps "
                    "minlen_values × maxcaas_values so you can pick the cluster-filtering "
                    "values that work best for your data.\n\n"
                    f"2. Set Cluster min length / Cluster max CAAS value (filter mode) on the "
                    f"Disambiguation tab to your chosen values, regenerate, then run "
                    f"{complete_name} — it reuses {expl_name}'s CAAS/Disambiguation output "
                    "directly (no recomputation) and runs the full downstream workflow "
                    "(Accumulation/VEP/Scoring/Enrichment, plus RER/FADE if you've enabled "
                    "them) on top of it.\n\n"
                    "Enjoy!",
                )
        else:
            batch_name, single_name = script_names("")
            batch_text = render_batch(self.project, single_runner_filename=single_name)
            single_text = render_single(self.project)
            scripts.append((batch_name, batch_text))
            scripts.append((single_name, single_text))

        self._show_preview(scripts)

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

    def _show_preview(self, scripts: list[tuple[str, str]]) -> None:
        preview = QWidget()  # standalone top-level window, not a child of MainWindow
        preview.setWindowTitle("Generated scripts preview")
        layout = QVBoxLayout(preview)

        preview_tabs = QTabWidget(preview)
        for title, text in scripts:
            editor = QPlainTextEdit()
            editor.setReadOnly(True)
            editor.setPlainText(text)
            editor.setFont(_monospace_font())
            preview_tabs.addTab(editor, title)
        layout.addWidget(preview_tabs)

        button_row = QHBoxLayout()
        save_btn = QPushButton("Save to disk...")
        save_btn.clicked.connect(
            lambda: self._save_generated_scripts(scripts)
        )
        button_row.addStretch(1)
        button_row.addWidget(save_btn)
        layout.addLayout(button_row)

        preview.resize(900, 700)
        preview.show()
        self._preview_window = preview  # keep a reference alive

    def _save_generated_scripts(self, scripts: list[tuple[str, str]]) -> None:
        repo_dir = self.project.general.repo_dir.strip()
        if not repo_dir:
            QMessageBox.warning(
                self,
                "Cannot save scripts",
                "Set General > Repo dir first — generated scripts are written there.",
            )
            return

        remote_host = self.project.general.remote_host.strip()
        target = remote_host or "the local filesystem"
        filenames_str = ", ".join([f[0] for f in scripts])
        reply = QMessageBox.question(
            self,
            "Save generated scripts",
            f"Write {filenames_str} to:\n{repo_dir}\n"
            f"on {target}?\n\nExisting files with the same names will be overwritten.",
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        try:
            if remote_host:
                base = repo_dir.rstrip("/")
                for filename, text in scripts:
                    remote.write_remote_file(remote_host, f"{base}/{filename}", text, mode="755")
            else:
                repo_path = Path(repo_dir)
                for filename, text in scripts:
                    out_path = repo_path / filename
                    out_path.write_text(text)
                    out_path.chmod(0o755)
        except (OSError, RemoteCheckError) as exc:
            QMessageBox.critical(self, "Failed to save scripts", str(exc))
            return

        QMessageBox.information(
            self,
            "Scripts saved",
            f"Wrote {filenames_str} to {repo_dir} on {target}.",
        )


def _monospace_font():
    from PySide6.QtGui import QFontDatabase

    return QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
