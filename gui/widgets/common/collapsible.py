#!/usr/bin/env python3
# collapsible.py — CollapsibleSection: a disclosure triangle + hideable content area.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
CollapsibleSection: a disclosure-triangle toggle that shows or hides a content widget.

Module tabs carry an essential field group plus a longer advanced one (fine-tuning
parameters most runs leave at their conf/*.config defaults, see ModuleTabWidget).
Qt has no built-in disclosure triangle, so this pairs a checkable QToolButton
(with an arrow) with a content QWidget whose visibility follows the button. The
section is collapsed by default, so a tab with many advanced fields does not
present them all as equally important.

Imported by: gui/widgets/common/module_tab.py
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QToolButton, QVBoxLayout, QWidget


class CollapsibleSection(QWidget):
    """Wrap `content` (any QWidget, typically a QGroupBox with a QFormLayout) so
    it's hidden behind a "▶ <title>" toggle button until clicked."""

    def __init__(self, title: str, content: QWidget, parent=None, *, expanded: bool = False):
        super().__init__(parent)
        self._content = content

        self.toggle_button = QToolButton()
        self.toggle_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.toggle_button.setStyleSheet(
            "QToolButton { border: 1px solid transparent; font-weight: bold; border-radius: 6px; padding: 6px 10px; font-size: 13px; }"
            "QToolButton:hover { background-color: rgba(148, 163, 184, 0.15); }"
        )
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.ArrowType.RightArrow)
        self.toggle_button.setCheckable(True)
        self.toggle_button.setText(title)
        self.toggle_button.toggled.connect(self._on_toggled)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toggle_button)
        layout.addWidget(content)

        self.toggle_button.setChecked(expanded)
        self._on_toggled(expanded)

    def _on_toggled(self, checked: bool) -> None:
        self.toggle_button.setArrowType(Qt.ArrowType.DownArrow if checked else Qt.ArrowType.RightArrow)
        self._content.setVisible(checked)

    def set_title(self, title: str) -> None:
        self.toggle_button.setText(title)
