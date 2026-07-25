#!/usr/bin/env python3
# about_tab.py — About tab: logo, repo, authorship, license, attributions.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Plain QLabel/QTextBrowser only — no QWebEngineView anywhere in this app (see
implementation plan §7): nothing here needs a web view, and QtWebEngine is by far
the biggest source of PyInstaller+AppImage packaging pain (sandboxed helper
process, .pak resource files, --no-sandbox requirements under a read-only mount).
"""

# ── Standard library ──────────────────────────────────────────────────────────
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QLabel, QTextBrowser, QVBoxLayout, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.about import AboutInfo

_GUI_PACKAGE_ROOT = Path(__file__).resolve().parents[2]  # .../gui/widgets/tabs -> repo root


class AboutTab(QWidget):
    def __init__(self, info: AboutInfo, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)

        logo_path = _GUI_PACKAGE_ROOT / info.logo_relpath
        if logo_path.is_file():
            pixmap = QPixmap(str(logo_path))
            if not pixmap.isNull():
                logo_label = QLabel()
                logo_label.setPixmap(
                    pixmap.scaledToHeight(120, Qt.TransformationMode.SmoothTransformation)
                )
                logo_label.setAlignment(Qt.AlignmentFlag.AlignHCenter)
                layout.addWidget(logo_label)

        attributions_html = "".join(f"<li>{a}</li>" for a in info.attributions)
        html = f"""
        <h2>{info.name} <span style="font-weight:normal">v{info.version}</span></h2>
        <p>{info.description}</p>
        <p><b>Repository:</b> <a href="{info.repo_url}">{info.repo_url}</a></p>
        <p><b>Author:</b> {info.author} (&lt;{info.author_email}&gt;)</p>
        <p><b>Institution:</b> {info.institution}</p>
        <p><b>License:</b> {info.license_name}</p>
        <p><b>Built on:</b></p>
        <ul>{attributions_html}</ul>
        """

        browser = QTextBrowser()
        browser.setOpenExternalLinks(True)
        browser.setHtml(html)
        layout.addWidget(browser, stretch=1)
