#!/usr/bin/env python3
# collapsible.py — CollapsibleSection: a disclosure triangle + hideable content area.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Most module tabs now carry a "basic" field group plus a much longer "advanced"
one (fine-tuning knobs most runs leave at their conf/*.config defaults — see
ModuleTabWidget). Qt has no built-in disclosure-triangle widget, so this is a
small QToolButton (arrow + checkable) paired with a content QWidget whose
visibility follows the button — collapsed by default so a tab with 20 advanced
fields doesn't read as 20 equally-important fields.
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
        self.toggle_button.setStyleSheet("QToolButton { border: none; font-weight: bold; }")
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
