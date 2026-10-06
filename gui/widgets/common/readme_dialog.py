#!/usr/bin/env python3
# readme_dialog.py — Interactive HTML README viewer dialog for PhyloPhere.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
ReadmeDialog: shows the repository README.md in a dialog with collapsible sections.

Renders through QWebEngineView when PySide6's QtWebEngine is installed (every
<details> section expanded) and falls back to a QTextBrowser otherwise, where
sections are toggled through "toggle:<id>" anchors. The Markdown is converted to
HTML by readme_html_builder.build_qtextbrowser_html in both cases. The README is
looked up in the repo directory from the General tab, then in the checkout that
contains this file, then in the working directory.

Imported by: gui/widgets/main_window.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QPushButton,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

# QtWebEngine is optional: without it the dialog uses a QTextBrowser.
try:
    from PySide6.QtWebEngineCore import QWebEnginePage
    from PySide6.QtWebEngineWidgets import QWebEngineView

    HAS_WEBENGINE = True
except ImportError:
    HAS_WEBENGINE = False

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.widgets.common.readme_html_builder import build_qtextbrowser_html

_REPO_DIR = Path(__file__).resolve().parents[3]  # .../gui/widgets/common -> repo root


class WebEngineNavPage(QWebEnginePage if HAS_WEBENGINE else object):
    """Custom WebEnginePage that opens external links in the default browser."""

    def acceptNavigationRequest(self, url, nav_type, is_main_frame):
        if nav_type == QWebEnginePage.NavigationType.NavigationTypeLinkClicked:
            url_str = url.toString()
            if url_str.startswith("http://") or url_str.startswith("https://"):
                QDesktopServices.openUrl(url)
                return False
        return super().acceptNavigationRequest(url, nav_type, is_main_frame)


class ReadmeDialog(QDialog):
    """Dialog displaying README.md with expandable <details> sections.

    The QTextBrowser fallback keeps the set of expanded section ids and re-renders
    on each toggle, so the dialog works with or without QtWebEngine.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        repo_dir: Path | str | None = None,
        lang: str = "en",
    ):
        super().__init__(parent)
        self.setWindowTitle("PhyloPhere — Documentation")
        self.resize(1000, 750)

        candidate_dirs = []
        if repo_dir and str(repo_dir).strip():
            candidate_dirs.append(Path(repo_dir))
        candidate_dirs.extend([_REPO_DIR, Path.cwd()])

        self._readme_path = None
        for d in candidate_dirs:
            p = (d / "README.md").resolve()
            if p.is_file():
                self._readme_path = p
                break

        if self._readme_path is None:
            self._readme_path = _REPO_DIR / "README.md"

        self._lang = lang
        self._expanded_sections: set[int] = set()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)

        self._readme_text = ""
        if self._readme_path and self._readme_path.is_file():
            try:
                self._readme_text = self._readme_path.read_text(encoding="utf-8")
            except Exception as e:
                self._readme_text = f"# Readme\n\nCould not load README.md: {e}"
        else:
            self._readme_text = "# Readme\n\nREADME.md file not found."

        if HAS_WEBENGINE:
            self.web_view = QWebEngineView(self)
            self.web_page = WebEngineNavPage(self.web_view)
            self.web_view.setPage(self.web_page)

            html_content, _ = build_qtextbrowser_html(
                self._readme_text, set(range(100))
            )
            base_url = QUrl.fromLocalFile(str(self._readme_path.parent) + "/")
            self.web_view.setHtml(html_content, base_url)
            layout.addWidget(self.web_view, stretch=1)
        else:
            # Fallback: QTextBrowser, with sections toggled through anchor clicks.
            self.browser = QTextBrowser(self)
            self.browser.setOpenLinks(False)
            self.browser.anchorClicked.connect(self._on_textbrowser_anchor_clicked)
            self._render_textbrowser_content()
            layout.addWidget(self.browser, stretch=1)

        # Footer: open in the system browser, close
        footer_layout = QHBoxLayout()
        footer_layout.setContentsMargins(0, 8, 0, 0)

        self.open_external_btn = QPushButton("Open in External Browser", self)
        self.open_external_btn.clicked.connect(self._open_external)

        self.close_btn = QPushButton("Close", self)
        self.close_btn.clicked.connect(self.accept)

        footer_layout.addWidget(self.open_external_btn)
        footer_layout.addStretch(1)
        footer_layout.addWidget(self.close_btn)
        layout.addLayout(footer_layout)

        self.retranslate(lang)

    def _render_textbrowser_content(self) -> None:
        if not hasattr(self, "browser"):
            return
        scroll_bar = self.browser.verticalScrollBar()
        scroll_pos = scroll_bar.value()

        html_content, _ = build_qtextbrowser_html(
            self._readme_text, self._expanded_sections
        )

        base_dir = (
            str(self._readme_path.parent) if self._readme_path else str(_REPO_DIR)
        )
        self.browser.setSearchPaths([base_dir])
        self.browser.setHtml(html_content)

        scroll_bar.setValue(scroll_pos)

    def _on_textbrowser_anchor_clicked(self, url: QUrl) -> None:
        url_str = url.toString()
        if url_str.startswith("toggle:"):
            sec_id = int(url_str.split(":")[1])
            if sec_id in self._expanded_sections:
                self._expanded_sections.remove(sec_id)
            else:
                self._expanded_sections.add(sec_id)
            self._render_textbrowser_content()
        elif url_str.startswith("#"):
            self.browser.scrollToAnchor(url_str.lstrip("#"))
        elif url.fragment():
            self.browser.scrollToAnchor(url.fragment())
        elif url_str.startswith("http://") or url_str.startswith("https://"):
            QDesktopServices.openUrl(url)

    def _open_external(self) -> None:
        if self._readme_path and self._readme_path.is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self._readme_path)))

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr

        self._lang = lang
        title_str = tr("Readme", lang)
        self.setWindowTitle(f"PhyloPhere — {title_str}")
        self.open_external_btn.setText(
            tr("Open in External Browser", lang)
            if tr("Open in External Browser", lang)
            != "Open in External Browser"
            else "Open in External Browser"
        )
        self.close_btn.setText(
            tr("Close", lang) if tr("Close", lang) != "Close" else "Close"
        )
