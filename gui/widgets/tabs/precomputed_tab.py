#!/usr/bin/env python3
# precomputed_tab.py — Precomputed Run tab: every global standalone/resume input, in one place.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Consolidates the "used when the producing module is disabled" fallback inputs that
used to be scattered one-per-tab (a disabled RER/FADE tab's own fallback box was
enabled but empty — the real fields live on PhenotypeRow per phenotype instead —
which looked like a dead end; see gui/models/precomputed.py's docstring for the
full story). Grouped by which module's absence each field substitutes for, so it
reads like "if module X is off, fill in its section here" rather than one flat
24-field form.

Per-phenotype fallbacks (scoring_rer_input, scoring_fade_summary_top/bottom) stay
on the Runtime tab's phenotype table — they vary per phenotype in a batch run, so
they're intentionally NOT here.
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QGroupBox, QLabel, QScrollArea, QVBoxLayout, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.precomputed import PrecomputedConfig
from gui.widgets.common.path_field import PathField

# (config field name, label, path kind: "file" | "dir")
_CAAS_FIELDS = [
    ("discovery_from", "Discovery output", "dir"),
    ("resample_from", "Resample output", "dir"),
    ("bootstrap_from", "Bootstrap output", "dir"),
]
_DISAMBIGUATION_FIELDS = [
    ("signification_from", "Signification output", "dir"),
    ("disambiguation_input", "Disambiguation input file", "file"),
    ("disambiguation_dir", "Disambiguation directory", "dir"),
    ("background_input", "Background input", "file"),
]
_ACCUMULATION_FIELDS = [
    ("accumulation_caas_input", "CAAS input file", "file"),
    ("accumulation_background_input", "Background input", "file"),
]
_RER_FIELDS = [
    ("rer_continuous_file", "Precomputed continuous-analysis RDS", "file"),
    ("rer_perms_file", "Precomputed permutation RDS", "file"),
]
_FADE_FIELDS = [
    ("fade_json_dir_top", "Precomputed *.FADE.json directory — top", "dir"),
    ("fade_json_dir_bottom", "Precomputed *.FADE.json directory — bottom", "dir"),
]
_VEP_FIELDS = [
    ("vep_caas_input", "Precomputed CAAS input", "file"),
]
_SCORING_FIELDS = [
    ("scoring_postproc_input", "Post-processing input", "file"),
    ("scoring_accum_dir", "Accumulation directory", "dir"),
    ("scoring_vep_primateai", "VEP PrimateAI scores", "file"),
    ("scoring_vep_cosmic", "VEP COSMIC scores", "file"),
    ("scoring_background_input", "Background input", "file"),
    ("caas_perms_file", "CAAS permulation RDS", "file"),
    ("scoring_fade_site_top", "FADE site-level fallback — top", "file"),
    ("scoring_fade_site_bottom", "FADE site-level fallback — bottom", "file"),
]
_ENRICHMENT_FIELDS = [
    ("accumulation_enrichment_gene_lists_input", "Accumulation gene lists input", "file"),
    ("posenrich_background_file", "POSENRICH background file", "file"),
]

_SECTIONS = [
    ("CAAS / Contrast Selection", "used when CAAS is disabled", _CAAS_FIELDS),
    ("Disambiguation", "used when Disambiguation is disabled", _DISAMBIGUATION_FIELDS),
    ("Accumulation", "used when Accumulation is disabled", _ACCUMULATION_FIELDS),
    ("RERconverge", "precomputed report input — renders the report without a live RER run", _RER_FIELDS),
    ("FADE", "precomputed report input — renders the report without a live FADE run", _FADE_FIELDS),
    ("VEP", "used when VEP is disabled", _VEP_FIELDS),
    ("Scoring", "used when Scoring is disabled, or when a producing module above is off", _SCORING_FIELDS),
    ("Enrichment", "used when Enrichment or an upstream module is disabled", _ENRICHMENT_FIELDS),
]


class PrecomputedTab(QWidget):
    changed = Signal()

    def __init__(self, config: PrecomputedConfig, parent=None):
        super().__init__(parent)
        self._config = config
        self._widgets: dict[str, PathField] = {}

        self._group_boxes: list[tuple[QGroupBox, QLabel, str, str]] = []
        self._form_labels: list[tuple[QLabel, str]] = []

        outer = QVBoxLayout(self)

        self.note_label = QLabel(
            "Every global standalone/resume input, in one place. Fill in a section only "
            "if you're skipping that module this run and feeding a prior run's output "
            "downstream instead. Per-phenotype fallbacks (RER/FADE summaries used by "
            "Scoring) stay on the Runtime tab's phenotype table, since they can differ "
            "per phenotype in a batch run."
        )
        self.note_label.setWordWrap(True)
        outer.addWidget(self.note_label)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll, stretch=1)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        for title, subtitle, fields in _SECTIONS:
            content_layout.addWidget(self._build_section(title, subtitle, fields))
        content_layout.addStretch(1)
        scroll.setWidget(content)

    def _build_section(self, title: str, subtitle: str, fields: list[tuple[str, str, str]]) -> QGroupBox:
        box = QGroupBox(title)
        layout = QVBoxLayout(box)

        subtitle_label = QLabel(subtitle)
        subtitle_label.setWordWrap(True)
        subtitle_label.setStyleSheet("color: #6b6b6b; font-style: italic;")
        layout.addWidget(subtitle_label)

        self._group_boxes.append((box, subtitle_label, title, subtitle))

        form = QFormLayout()
        for name, label_text, kind in fields:
            widget = PathField(mode=kind)
            widget.set_text(getattr(self._config, name))
            widget.textChanged.connect(lambda v, n=name: self._set_field(n, v))
            lbl = QLabel(label_text)
            form.addRow(lbl, widget)
            self._form_labels.append((lbl, label_text))
            self._widgets[name] = widget
        layout.addLayout(form)

        return box

    def _set_field(self, name: str, value: str) -> None:
        setattr(self._config, name, value)
        self.changed.emit()

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "note_label"):
            self.note_label.setText(tr(
                "Every global standalone/resume input, in one place. Fill in a section only "
                "if you're skipping that module this run and feeding a prior run's output "
                "downstream instead. Per-phenotype fallbacks (RER/FADE summaries used by "
                "Scoring) stay on the Runtime tab's phenotype table, since they can differ "
                "per phenotype in a batch run.",
                lang
            ))
        for box, sub_lbl, title, subtitle in getattr(self, "_group_boxes", []):
            box.setTitle(tr(title, lang))
            sub_lbl.setText(tr(subtitle, lang))
        for lbl, orig_text in getattr(self, "_form_labels", []):
            lbl.setText(tr(orig_text, lang))
