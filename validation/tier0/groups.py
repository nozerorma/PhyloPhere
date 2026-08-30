"""Amino-acid grouping schemes, vendored from the pipeline.

Copied verbatim from
``subworkflows/CT_DISAMBIGUATION/local/src/biochem/grouping.py`` (which is itself
kept aligned with ``subworkflows/CT/local/modules/caap_id.py``). Vendored rather
than imported so Tier 0 has no import-path dependency on the pipeline tree; if
the pipeline schemes change, update this file and note it in PROVENANCE.md.

Used only to build the ``grouped_caap`` planted mechanism: a convergent shift to
a shared physicochemical *class* with the residue left free inside it, so the
strict ``US`` scheme sees divergent residues while ``GS*`` sees convergence.
"""

from __future__ import annotations

US = {a: a for a in "ACDEFGHIKLMNPQRSTVWY"}

GS1 = {
    "C": "t", "V": "t", "A": "n", "G": "n", "P": "n", "S": "n",
    "N": "p", "D": "p", "Q": "p", "E": "p", "R": "b", "H": "b", "K": "b",
    "I": "h", "L": "h", "M": "h", "F": "h", "W": "h", "Y": "h", "T": "o",
}
GS2 = {
    "C": "c", "A": "s", "G": "s", "V": "s", "D": "a", "E": "a",
    "N": "n", "Q": "n", "H": "n", "W": "n", "R": "b", "K": "b",
    "I": "h", "L": "h", "F": "h", "P": "h", "Y": "x", "M": "x", "T": "x", "S": "x",
}
GS3 = {
    "C": "c", "A": "n", "G": "n", "P": "n", "S": "n", "T": "n",
    "N": "s", "D": "s", "Q": "s", "E": "s", "R": "b", "H": "b", "K": "b",
    "I": "l", "L": "l", "M": "l", "V": "l", "F": "g", "W": "g", "Y": "g",
}
GS4 = {
    "C": "c", "A": "h", "I": "h", "L": "h", "V": "h", "S": "o", "T": "o",
    "N": "p", "Q": "p", "D": "a", "E": "a", "R": "b", "H": "b",
    "G": "g", "P": "r", "K": "k", "M": "m", "F": "f", "W": "y", "Y": "y",
}

SCHEMES: dict[str, dict[str, str]] = {"US": US, "GS1": GS1, "GS2": GS2, "GS3": GS3, "GS4": GS4}

GS_SCHEMES = ("GS1", "GS2", "GS3", "GS4")


def group_members(scheme: str, label: str) -> list[str]:
    """Residues that share group ``label`` under ``scheme`` (sorted)."""
    return sorted(a for a, g in SCHEMES[scheme].items() if g == label)


def groups_of(scheme: str, *, min_size: int = 2) -> dict[str, list[str]]:
    """All groups of ``scheme`` with at least ``min_size`` members."""
    out: dict[str, list[str]] = {}
    for a, g in SCHEMES[scheme].items():
        out.setdefault(g, []).append(a)
    return {g: sorted(v) for g, v in out.items() if len(v) >= min_size}
