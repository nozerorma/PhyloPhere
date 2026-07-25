#!/usr/bin/env python3
# resources_tab.py — Resources tab: local/slurm profile-level resource caps.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFormLayout, QGroupBox, QLabel, QLineEdit, QVBoxLayout, QWidget

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.resources import ResourcesConfig


class ResourcesTab(QWidget):
    changed = Signal()

    def __init__(self, config: ResourcesConfig, parent=None):
        super().__init__(parent)
        self._config = config

        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "These map to nextflow.config's --max_cpus/--max_memory/--max_time overrides on "
                "the generated `nextflow run` invocation. conf/resources.config's per-process "
                "labels (the pipeline's internal resource requests) are not exposed here — they "
                "were confirmed identical between the local and slurm profiles and are pipeline "
                "internals, not something to tune per run."
            )
        )
        layout.addWidget(self._build_group("Local defaults", "local"))
        layout.addWidget(self._build_group("SLURM defaults", "slurm"))
        layout.addStretch(1)

    def _build_group(self, title: str, prefix: str) -> QGroupBox:
        box = QGroupBox(title)
        form = QFormLayout(box)

        cpus_attr = f"{prefix}_max_cpus"
        mem_attr = f"{prefix}_max_memory"
        time_attr = f"{prefix}_max_time"

        cpus_field = QLineEdit(getattr(self._config, cpus_attr))
        cpus_field.textChanged.connect(lambda v: self._on_changed(cpus_attr, v))
        form.addRow("Max CPUs (--max_cpus)", cpus_field)

        mem_field = QLineEdit(getattr(self._config, mem_attr))
        mem_field.textChanged.connect(lambda v: self._on_changed(mem_attr, v))
        form.addRow("Max memory (--max_memory)", mem_field)

        time_field = QLineEdit(getattr(self._config, time_attr))
        time_field.textChanged.connect(lambda v: self._on_changed(time_attr, v))
        form.addRow("Max time (--max_time)", time_field)

        setattr(self, f"{prefix}_cpus_field", cpus_field)
        setattr(self, f"{prefix}_mem_field", mem_field)
        setattr(self, f"{prefix}_time_field", time_field)

        return box

    def _on_changed(self, attr: str, value: str) -> None:
        setattr(self._config, attr, value)
        self.changed.emit()
