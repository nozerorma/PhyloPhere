#!/usr/bin/env python3
# support_fmt.py — Format residue support tallies as 'L:3,S:2' strings.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/convergence/

"""
Support-string formatting for cross-hypothesis pooling.

Imported by: src/convergence/disambiguate_single.py, src/convergence/fop_pool.py
Inputs: a {residue or tag: count} dict
Outputs: a comma-joined 'key:count' string (empty when there is no count)
"""

from __future__ import annotations

from typing import Dict


def fmt_support(counts: Dict[str, int]) -> str:
    """'L:3,S:2'-style string: count-descending, then alphabetical tiebreak. Zero counts are dropped."""
    counts = {k: v for k, v in counts.items() if v}
    if not counts:
        return ""
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{k}:{v}" for k, v in ordered)
