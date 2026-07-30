#!/usr/bin/env python3
# resources_tab.py — Resources tab: local/slurm profile-level resource caps.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.resources import ResourcesConfig
from gui.widgets.resource_table.widget import ResourceOverrideTableWidget


class ResourcesTab(QWidget):
    changed = Signal()

    def __init__(self, config: ResourcesConfig, parent=None):
        super().__init__(parent)
        self._config = config

        layout = QVBoxLayout(self)
        self.top_note = QLabel(
            "These map to nextflow.config's --max_cpus/--max_memory/--max_time overrides on "
            "the generated `nextflow run` invocation."
        )
        layout.addWidget(self.top_note)
        layout.addWidget(self._build_group("Local defaults", "local"))
        layout.addWidget(self._build_group("SLURM defaults", "slurm"))

        self.overrides_note = QLabel(
            "Advanced: per-process cpus/memory overrides, rendered into a config file "
            "loaded via `-c` alongside the generated run scripts (see run_single.sh.j2) — "
            "Nextflow only reads these from a config file, never from --flags. Start from "
            "a preset below and hand-edit as needed; conf/resources.config's own defaults "
            "apply to any process left out of this table."
        )
        layout.addWidget(self.overrides_note)
        self.overrides_table = ResourceOverrideTableWidget(config.process_overrides)
        self.overrides_table.changed.connect(self.changed)
        layout.addWidget(self.overrides_table)

    def _build_group(self, title: str, prefix: str) -> QGroupBox:
        box = QGroupBox(title)
        form = QFormLayout(box)

        cpus_attr = f"{prefix}_max_cpus"
        mem_attr = f"{prefix}_max_memory"
        time_attr = f"{prefix}_max_time"

        cpus_field = QLineEdit(getattr(self._config, cpus_attr))
        cpus_field.textChanged.connect(lambda v: self._on_changed(cpus_attr, v))
        cpus_lbl = QLabel("Max CPUs (--max_cpus)")
        form.addRow(cpus_lbl, cpus_field)

        mem_field = QLineEdit(getattr(self._config, mem_attr))
        mem_field.textChanged.connect(lambda v: self._on_changed(mem_attr, v))
        mem_lbl = QLabel("Max memory (--max_memory)")
        form.addRow(mem_lbl, mem_field)

        time_field = QLineEdit(getattr(self._config, time_attr))
        time_field.textChanged.connect(lambda v: self._on_changed(time_attr, v))
        time_lbl = QLabel("Max time (--max_time)")
        form.addRow(time_lbl, time_field)

        setattr(self, f"{prefix}_box", box)
        setattr(self, f"{prefix}_cpus_lbl", cpus_lbl)
        setattr(self, f"{prefix}_mem_lbl", mem_lbl)
        setattr(self, f"{prefix}_time_lbl", time_lbl)
        setattr(self, f"{prefix}_cpus_field", cpus_field)
        setattr(self, f"{prefix}_mem_field", mem_field)
        setattr(self, f"{prefix}_time_field", time_field)

        return box

    def _on_changed(self, attr: str, value: str) -> None:
        setattr(self._config, attr, value)
        self.changed.emit()

    def retranslate(self, lang: str = "en") -> None:
        from gui.i18n import tr
        if hasattr(self, "local_box"):
            self.local_box.setTitle(tr("Local defaults", lang))
        if hasattr(self, "slurm_box"):
            self.slurm_box.setTitle(tr("SLURM defaults", lang))
        for prefix in ("local", "slurm"):
            cpus_lbl = getattr(self, f"{prefix}_cpus_lbl", None)
            if cpus_lbl:
                cpus_lbl.setText(tr("Max CPUs (--max_cpus)", lang))
            mem_lbl = getattr(self, f"{prefix}_mem_lbl", None)
            if mem_lbl:
                mem_lbl.setText(tr("Max memory (--max_memory)", lang))
            time_lbl = getattr(self, f"{prefix}_time_lbl", None)
            if time_lbl:
                time_lbl.setText(tr("Max time (--max_time)", lang))
        if hasattr(self, "overrides_table"):
            self.overrides_table.retranslate(lang)
        if hasattr(self, "top_note"):
            self.top_note.setText(tr(
                "These map to nextflow.config's --max_cpus/--max_memory/--max_time overrides on "
                "the generated `nextflow run` invocation.",
                lang,
            ))
        if hasattr(self, "overrides_note"):
            self.overrides_note.setText(tr(
                "Advanced: per-process cpus/memory overrides, rendered into a config file "
                "loaded via `-c` alongside the generated run scripts (see run_single.sh.j2) — "
                "Nextflow only reads these from a config file, never from --flags. Start from "
                "a preset below and hand-edit as needed; conf/resources.config's own defaults "
                "apply to any process left out of this table.",
                lang,
            ))
