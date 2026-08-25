#!/usr/bin/env python3
# regenerate_dialog.py — RegenerateReportsDialog: standalone "regenerate HTML
# reports from an existing output directory" dialog, independent of any
# project configuration. Scans a chosen directory for the pipeline's known
# report HTML files and builds ready-to-run shell scripts that re-render them
# from the files already produced there (see gui/generation/report_registry.py).
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from gui import remote
from gui.generation.report_registry import (
    UNSUPPORTED_NOTE,
    DetectedReport,
    build_script,
    detect_reports,
)
from gui.widgets.common import remote_context
from gui.widgets.common.path_field import PathField

# A remote outdir can hold a lot of files; a plain existence/glob check needs
# more room than the 15s default used for the small single-round-trip checks
# elsewhere (gui.remote.DEFAULT_TIMEOUT).
_REMOTE_SCAN_TIMEOUT = 60

# Same yellow translucent disclaimer style used for module-tab "if disabled"
# banners (gui/widgets/common/module_tab.py) — reused here for consistency
# rather than inventing a second warning-box style.
_DISCLAIMER_STYLESHEET = (
    "QLabel { background-color: rgba(234, 179, 8, 0.12); color: #d97706; "
    "border: 1px solid rgba(234, 179, 8, 0.3); border-radius: 6px; padding: 8px 12px; font-weight: 500; }"
)

_STATUS_COL, _REPORT_COL, _SLOTS_COL = 0, 1, 2


class RegenerateReportsDialog(QDialog):
    """Standalone: takes no ProjectConfig, only a directory. Emits nothing —
    the caller (MainWindow) reads self.generated_scripts after exec()."""

    def __init__(self, parent: QWidget | None = None, repo_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Regenerate HTML Reports")
        self.resize(760, 560)

        self._detected: list[DetectedReport] = []
        self.generated_scripts: list[tuple[str, str]] = []

        layout = QVBoxLayout(self)

        disclaimer = QLabel(
            "Standalone regeneration: re-renders HTML reports from files already "
            "produced in the chosen output directory. Does not rerun any pipeline "
            "computation, and ignores the current project's settings on the other tabs."
        )
        disclaimer.setWordWrap(True)
        disclaimer.setStyleSheet(_DISCLAIMER_STYLESHEET)
        layout.addWidget(disclaimer)

        remote_host = remote_context.get_remote_host()
        self.remote_status_label = QLabel(
            f"Remote host active: {remote_host} — scanning and browsing happen over SSH."
            if remote_host
            else "No remote host set (General tab) — scanning the local filesystem."
        )
        self.remote_status_label.setStyleSheet("QLabel { color: palette(mid); font-size: 11px; }")
        layout.addWidget(self.remote_status_label)

        form = QGroupBox("Source")
        form_layout = QVBoxLayout(form)

        outdir_row = QHBoxLayout()
        outdir_row.addWidget(QLabel("Output directory:"))
        self.outdir_field = PathField(mode="dir", placeholder="Pipeline --outdir to scan for html_reports/")
        outdir_row.addWidget(self.outdir_field, 1)
        form_layout.addLayout(outdir_row)

        repo_row = QHBoxLayout()
        repo_row.addWidget(QLabel("Repo directory:"))
        self.repo_field = PathField(mode="dir", placeholder="PhyloPhere checkout (for the Rmd sources)")
        self.repo_field.set_text(repo_dir)
        repo_row.addWidget(self.repo_field, 1)
        form_layout.addLayout(repo_row)

        opts_row = QHBoxLayout()
        opts_row.addWidget(QLabel("Trait name override:"))
        self.traitname_field = QLineEdit()
        self.traitname_field.setPlaceholderText("leave blank to auto-guess per report")
        opts_row.addWidget(self.traitname_field, 1)
        self.singularity_check = QCheckBox("Use Singularity/Apptainer entrypoint")
        opts_row.addWidget(self.singularity_check)
        form_layout.addLayout(opts_row)

        scan_row = QHBoxLayout()
        scan_row.addStretch(1)
        self.scan_button = QPushButton("Scan")
        self.scan_button.clicked.connect(self._scan)
        scan_row.addWidget(self.scan_button)
        form_layout.addLayout(scan_row)

        layout.addWidget(form)

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["", "Report", "Inputs"])
        self.table.horizontalHeader().setSectionResizeMode(_REPORT_COL, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(_SLOTS_COL, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table, 1)

        note = QLabel(UNSUPPORTED_NOTE)
        note.setWordWrap(True)
        note.setStyleSheet("QLabel { color: palette(mid); font-size: 11px; }")
        layout.addWidget(note)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        generate_btn = QPushButton("Generate Script(s)")
        generate_btn.clicked.connect(self._generate)
        button_row.addWidget(generate_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        button_row.addWidget(close_btn)
        layout.addLayout(button_row)

    # ── Actions ──────────────────────────────────────────────────────────

    def _scan(self) -> None:
        outdir_str = self.outdir_field.text().strip()
        if not outdir_str:
            QMessageBox.warning(self, "Cannot scan", "Set the output directory first.")
            return

        # Read at scan time, not construction time — matches PathField's own
        # convention (see remote_context.py) so a host set/cleared on the
        # General tab after this dialog opened is still respected.
        host = remote_context.get_remote_host()
        repo_dir = Path(self.repo_field.text().strip() or ".")

        if host:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            try:
                remote_entries = remote.list_all_remote(host, outdir_str, timeout=_REMOTE_SCAN_TIMEOUT)
            except remote.RemoteCheckError as exc:
                QMessageBox.critical(self, "Remote scan failed", str(exc))
                return
            finally:
                QApplication.restoreOverrideCursor()
            if not any(rel == "html_reports" or rel.startswith("html_reports/") for rel in remote_entries):
                QMessageBox.warning(
                    self,
                    "Cannot scan",
                    f"No html_reports/ directory found under:\n{outdir_str}\non {host}\n\n"
                    "This should be the pipeline's --outdir, not a module subdirectory.",
                )
                return
            self._detected = detect_reports(
                outdir_str, repo_dir, traitname_override=self.traitname_field.text(), remote_entries=remote_entries
            )
        else:
            if not (Path(outdir_str) / "html_reports").is_dir():
                QMessageBox.warning(
                    self,
                    "Cannot scan",
                    f"No html_reports/ directory found under:\n{outdir_str}\n\n"
                    "This should be the pipeline's --outdir, not a module subdirectory.",
                )
                return
            self._detected = detect_reports(outdir_str, repo_dir, traitname_override=self.traitname_field.text())

        self._populate_table()

        if not self._detected:
            QMessageBox.information(
                self, "Scan complete", "No known report HTML files found under html_reports/."
            )

    def _populate_table(self) -> None:
        self.table.setRowCount(len(self._detected))
        for row, report in enumerate(self._detected):
            check = QCheckBox()
            check.setChecked(report.runnable)
            check.setEnabled(report.runnable)
            cell = QWidget()
            cell_layout = QHBoxLayout(cell)
            cell_layout.setContentsMargins(4, 0, 4, 0)
            cell_layout.addWidget(check)
            cell_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.table.setCellWidget(row, _STATUS_COL, cell)
            self.table._checkboxes = getattr(self.table, "_checkboxes", {})
            self.table._checkboxes[row] = check

            name_item = QTableWidgetItem(report.display_name)
            name_item.setToolTip(str(report.html_path))
            self.table.setItem(row, _REPORT_COL, name_item)

            if report.runnable:
                slot_text = ", ".join(s.label for s in report.slots if s.path is not None)
                tooltip = "\n".join(f"{s.label}: {s.path}" for s in report.slots)
            else:
                slot_text = "MISSING: " + ", ".join(report.missing_required)
                tooltip = "\n".join(
                    f"{s.label}: {s.path if s.path else 'not found — browse to it manually'}" for s in report.slots
                )
            slots_item = QTableWidgetItem(slot_text)
            slots_item.setToolTip(tooltip)
            if not report.runnable:
                slots_item.setForeground(Qt.GlobalColor.darkRed)
            self.table.setItem(row, _SLOTS_COL, slots_item)

        self.table.resizeRowsToContents()

    def _generate(self) -> None:
        if not self._detected:
            QMessageBox.warning(self, "Nothing to generate", "Scan a directory first.")
            return

        checkboxes = getattr(self.table, "_checkboxes", {})
        selected = [
            report for row, report in enumerate(self._detected)
            if checkboxes.get(row) is not None and checkboxes[row].isChecked()
        ]
        if not selected:
            QMessageBox.warning(self, "Nothing selected", "Check at least one runnable report.")
            return

        repo_dir = Path(self.repo_field.text().strip() or ".")
        use_singularity = self.singularity_check.isChecked()

        scripts: list[tuple[str, str]] = []
        errors: list[str] = []
        for report in selected:
            if not report.runnable:
                errors.append(f"{report.display_name}: missing {', '.join(report.missing_required)}")
                continue
            try:
                text = build_script(report, repo_dir, use_singularity)
            except Exception as exc:  # noqa: BLE001 — surface any per-report build failure, keep going
                errors.append(f"{report.display_name}: {exc}")
                continue
            scripts.append((report.script_filename, text))

        if errors:
            QMessageBox.warning(
                self, "Some reports were skipped", "\n".join(f"• {e}" for e in errors)
            )

        if not scripts:
            return

        self.generated_scripts = scripts
        self.accept()
