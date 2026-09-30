"""Untrimmed alignment columns of the trimmed-alignment positions, read from the trimmer's MAP tables.

Discovery positions are 0-based columns of the trimmed alignment. The MAP table of a gene has one row per
column of the untrimmed alignment (``ori_codon_col``, 1-based) with ``status`` (``selected`` | ``removed``) and,
for the selected columns, ``prot_ali_col``: the 1-based column of the trimmed alignment. A position ``p`` therefore
sits at untrimmed column ``ori_of_prot[p + 1]``. Removed columns cannot hold a CAAS; measured in untrimmed columns
they count towards a window's span and not towards its CAAS count (see :func:`core.postproc.ctrain`).

Files are named ``<gene>[.<version>].<species><tail>`` and a gene is found by the part of the name before the first
``.``, whatever its reference species. Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional, Tuple

__all__ = ["read_map", "column_map", "index_files", "find_file", "gene_columns"]


def read_map(path: str) -> Tuple[List[bool], Dict[int, int]]:
    """(removed flag per untrimmed column, trimmed column -> untrimmed column), both 1-based.

    Raises ValueError when the rows are not numbered 1..n, a status is neither ``selected`` nor ``removed``, or the
    selected columns are not numbered 1..n in order.
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    removed, ori_of_prot = [], {}
    for i, r in enumerate(rows, 1):
        if int(r["ori_codon_col"]) != i:
            raise ValueError("ori_codon_col is not 1..n")
        status = r["status"]
        if status not in ("selected", "removed"):
            raise ValueError("status outside {selected, removed}")
        removed.append(status == "removed")
        if status == "selected":
            if int(r["prot_ali_col"]) != len(ori_of_prot) + 1:
                raise ValueError("prot_ali_col of the selected columns is not 1..n in order")
            ori_of_prot[len(ori_of_prot) + 1] = i
    return removed, ori_of_prot


def column_map(path: str) -> Dict[int, int]:
    """0-based trimmed position -> untrimmed 1-based column."""
    return {p - 1: c for p, c in read_map(path)[1].items()}


def index_files(directory: str, tail: str) -> Dict[str, Optional[str]]:
    """Gene key (file name before the first '.') -> path, for the files ending in ``tail``.

    A key that several files share maps to None.
    """
    out: Dict[str, Optional[str]] = {}
    for name in os.listdir(directory):
        if name.endswith(tail):
            key = name[: -len(tail)].split(".")[0]
            out[key] = None if key in out else os.path.join(directory, name)
    return out


def find_file(index: Dict[str, Optional[str]], gene: str, tail: str) -> str:
    """Path of the gene's file; FileNotFoundError when it has none, ValueError when it has several."""
    path = index.get(gene, False)
    if path is False:
        raise FileNotFoundError(f"no file ending in {tail} for {gene}")
    if path is None:
        raise ValueError(f"several files ending in {tail} for {gene}")
    return path


def gene_columns(index: Dict[str, Optional[str]], gene: str, tail: str) -> Optional[Dict[int, int]]:
    """:func:`column_map` of the gene's MAP file, or None when the gene has no file.

    A gene with several files or an unreadable or inconsistent file raises, so a wrong map is never used silently.
    """
    try:
        path = find_file(index, gene, tail)
    except FileNotFoundError:
        return None
    return column_map(path)
