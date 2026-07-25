#!/usr/bin/env python3
# specs.py — FieldSpec/ModuleTabSpec: declarative field lists for ModuleTabWidget.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Each of the 9 module tab files instantiates one ModuleTabSpec with its curated
field list (see implementation plan §5) rather than hand-rolling a QVBoxLayout —
this is the mechanism that keeps ~90% of module-tab structure shared.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field

FieldKind = str  # "bool" | "str" | "path_file" | "path_dir" | "choice"


@dataclass(frozen=True, kw_only=True)
class FieldSpec:
    name: str  # attribute name on the module's config dataclass
    label: str
    kind: FieldKind = "str"
    choices: tuple[str, ...] = ()  # only used when kind == "choice"
    placeholder: str = ""


@dataclass(frozen=True, kw_only=True)
class ModuleTabSpec:
    title: str
    blurb: str
    disclaimer: str  # shown when the module is disabled — what downstream needs
    essential_fields: tuple[FieldSpec, ...] = field(default_factory=tuple)
    fallback_fields: tuple[FieldSpec, ...] = field(default_factory=tuple)
    enabled_field: str = "enabled"  # attribute name for the top enable toggle
