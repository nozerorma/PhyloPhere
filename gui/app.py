#!/usr/bin/env python3
# app.py — QApplication bootstrap: exception hook, app metadata.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
import sys
import traceback
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication, QMessageBox

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.theme import apply_theme
from gui.widgets.main_window import MainWindow

_REPO_DIR = Path(__file__).resolve().parent.parent
_ICON_PATH = _REPO_DIR / "res" / "icon.png"  # app/taskbar icon — About tab uses res/logo.png instead


def _install_exception_hook() -> None:
    def handle(exc_type, exc_value, exc_tb) -> None:
        message = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
        print(message, file=sys.stderr)
        QMessageBox.critical(None, "Unexpected error", message)

    sys.excepthook = handle


def run(initial_project_path: Path | None = None) -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("PhyloPhere Runner GUI")
    app.setOrganizationName("IBE-UPF")
    # Without this, Wayland compositors (GNOME/KDE) can't match the running window
    # to install_gui_launcher.sh's phylophere-gui.desktop entry — the taskbar/task
    # switcher falls back to a generic placeholder icon instead of res/logo.png.
    # Must match that .desktop file's basename exactly (no .desktop suffix here).
    app.setDesktopFileName("phylophere-gui")
    if _ICON_PATH.is_file():
        app.setWindowIcon(QIcon(str(_ICON_PATH)))
    apply_theme(app)
    _install_exception_hook()

    window = MainWindow()
    if initial_project_path is not None:
        window.open_project_from_path(initial_project_path)
    window.resize(1100, 800)
    window.show()

    return app.exec()
