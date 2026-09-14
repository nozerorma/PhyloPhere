#!/usr/bin/env python3
"""Shared support-string formatting for cross-hypothesis/cross-scheme pooling.

Vendored here (not imported) because CT_DISAMBIGUATION and CT_POSTPROC run in
separate Nextflow process work dirs and cannot share a `local/src` tree.
"""

from __future__ import annotations

from typing import Dict


def fmt_support(counts: Dict[str, int]) -> str:
    """'L:3,S:2'-style string: count-descending, then alphabetical tiebreak.

    Mirrors CT_POSTPROC's residue_descriptors._fmt_support convention.
    """
    counts = {k: v for k, v in counts.items() if v}
    if not counts:
        return ""
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{k}:{v}" for k, v in ordered)
