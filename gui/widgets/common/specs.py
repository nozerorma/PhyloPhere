#!/usr/bin/env python3
# specs.py — FieldSpec/ModuleTabSpec: declarative field lists for ModuleTabWidget.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
FieldSpec and ModuleTabSpec: declarative description of a module tab.

Each module tab file declares one ModuleTabSpec holding its FieldSpec lists, and
ModuleTabWidget builds the layout from it, so the structure shared by all module
tabs lives in one place instead of in per-tab layout code.

Imported by: gui/widgets/common/module_tab.py and the module tab files in
gui/widgets/tabs/
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field
from typing import Literal

FieldKind = str  # "bool" | "str" | "path_file" | "path_dir" | "choice" | "multichoice" | "choice_with_other" | "section"

# required: no usable default exists downstream and the pipeline fails without it
#   (gui/generation/validate.py is the reference for this tier).
# default: has a working default, but changing it affects statistical/scientific
#   validity in a way that needs real understanding (seed, FDR thresholds,
#   permulation counts, model-selection parameters, ...).
# optional: auxiliary paths, standalone/precomputed overrides, cosmetic/report
#   parameters — safe to leave blank or default with no scientific consequence.
Importance = Literal["required", "default", "optional"]


@dataclass(frozen=True, kw_only=True)
class FieldSpec:
    name: str = ""  # attribute name on the module's config dataclass (empty for sections)
    label: str
    kind: FieldKind = "str"
    importance: Importance = "optional"  # see Importance above
    choices: tuple[str, ...] = ()  # "choice"/"multichoice": raw stored values. "choice_with_other":
    # display labels, with the LAST entry being the free-text "other" sentinel (see choice_other_values).
    choice_other_values: tuple[str, ...] = ()  # "choice_with_other" only: stored value for each of
    # `choices` except the last (the "other" sentinel) — index-aligned with choices[:-1]. The stored
    # value when "other" is picked is whatever the user types in the revealed free-text field.
    editable: bool = False  # "choice" only: render as an editable QComboBox pre-populated with
    # `choices` as presets, but still accepting free numeric/text entry (e.g. min_divergent_fraction)
    # instead of a strict closed dropdown.
    placeholder: str = ""
    help: str = ""  # tooltip text; shown on both the label and the input widget


def Section(label: str) -> FieldSpec:
    """Build a section sub-header entry, mirroring the comment blocks of conf/*.config."""
    return FieldSpec(name="", label=label, kind="section")


@dataclass(frozen=True, kw_only=True)
class ModuleTabSpec:
    title: str
    blurb: str
    disclaimer: str  # shown when the module is disabled — what downstream needs
    essential_fields: tuple[FieldSpec, ...] = field(default_factory=tuple)
    # Fine-tuning parameters most runs leave at their conf/*.config default, shown
    # in a collapsed "Advanced parameters" disclosure rather than next to
    # essential_fields (see gui/widgets/common/collapsible.py).
    advanced_fields: tuple[FieldSpec, ...] = field(default_factory=tuple)
    fallback_fields: tuple[FieldSpec, ...] = field(default_factory=tuple)
    enabled_field: str = "enabled"  # attribute name for the top enable toggle
