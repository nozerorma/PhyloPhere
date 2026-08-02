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
from PySide6.QtWidgets import QHBoxLayout, QLabel, QTextBrowser, QVBoxLayout, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.about import AboutInfo

_GUI_PACKAGE_ROOT = Path(__file__).resolve().parents[3]  # .../gui/widgets/tabs -> repo root


class AboutTab(QWidget):
    def __init__(self, info: AboutInfo, parent=None):
        super().__init__(parent)
        self._info = info

        layout = QVBoxLayout(self)

        logo_path = _GUI_PACKAGE_ROOT / info.logo_relpath
        if logo_path.is_file():
            pixmap = QPixmap(str(logo_path))
            if not pixmap.isNull():
                logo_container = QWidget()
                logo_container.setStyleSheet(
                    "QWidget { background-color: #ffffff; border: 1px solid #d4d4d8; border-radius: 8px; }"
                )
                container_layout = QHBoxLayout(logo_container)
                container_layout.setContentsMargins(12, 12, 12, 12)

                logo_label = QLabel()
                logo_label.setStyleSheet("background-color: transparent; border: none;")
                logo_label.setPixmap(
                    pixmap.scaled(
                        512,
                        512,
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation,
                    )
                )
                logo_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
                container_layout.addWidget(logo_label)
                container_layout.addStretch(1)
                layout.addWidget(logo_container)

        self.browser = QTextBrowser()
        self.browser.setOpenExternalLinks(True)
        layout.addWidget(self.browser, stretch=1)
        self.retranslate("en")

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "browser") and hasattr(self, "_info"):
            attributions_html = "".join(f"<li>{tr(a, lang)}</li>" for a in self._info.attributions)
            description_text = tr(self._info.description, lang)
            repo_label = tr("Repository", lang)
            author_label = tr("Author", lang)
            institution_label = tr("Institution", lang)
            license_label = tr("License", lang)
            built_on_label = tr("Built on", lang)
            html = f"""
            <h2>{self._info.name} <span style="font-weight:normal">v{self._info.version}</span></h2>
            <p>{description_text}</p>
            <p><b>{repo_label}:</b> <a href="{self._info.repo_url}">{self._info.repo_url}</a></p>
            <p><b>{author_label}:</b> {self._info.author} (&lt;{self._info.author_email}&gt;)</p>
            <p><b>{institution_label}:</b> {self._info.institution}</p>
            <p><b>{license_label}:</b> {self._info.license_name}</p>
            <p><b>{built_on_label}:</b></p>
            <ul>{attributions_html}</ul>
            """
            self.browser.setHtml(html)

