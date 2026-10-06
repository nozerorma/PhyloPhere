#!/usr/bin/env python3
# residue_descriptors.py — Residues of the extant contrast species at each CAAS position.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Extant-species residue tally for the disambiguation table: for each (Gene, Position), how
many foreground (top) and background (bottom) contrast species carry each residue.

Imported by: prepare_postproc_input.py (add_species_tally), which runs in the input-prep
             step (CAAS_PREPARE_POSTPROC_INPUT), before the cluster and gene filters
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd

# ── Extant-species residue tally (alignment-based) ────────────────────────────
# Counts the species of the foreground / background contrast that carry each residue
# at the alignment column (an "X of N species" check on the call). It does not depend on
# the hypothesis.

SPECIES_TALLY_COLUMNS = (
    "top_species_residues",     # e.g. "N:15,A:8" -- FG contrast species, count-desc
    "bottom_species_residues",  # e.g. "A:52,S:2" -- BG contrast species
    "n_top_species",            # FG contrast species with a non-gap residue here
    "n_bottom_species",         # BG contrast species with a non-gap residue here
)

_TALLY_EMPTY = {"top_species_residues": "", "bottom_species_residues": "",
                "n_top_species": "", "n_bottom_species": ""}

_GAP_CHARS = set("-.")


def _read_species_list(path) -> List[str]:
    """Species of a one-per-line file (a "species" header line is skipped); [] when absent or unreadable."""
    if not path:
        return []
    try:
        with open(path) as fh:
            return [ln.strip() for ln in fh if ln.strip() and ln.strip() != "species"]
    except OSError:
        return []


def _index_alignment_dir(alignment_dir) -> Dict[str, str]:
    """gene symbol -> alignment file path. Flat directory; the gene is the file name up to the
    first '.', the convention of CT_ACCUMULATION's concatenate.py."""
    import glob
    import os
    out: Dict[str, str] = {}
    for f in glob.glob(os.path.join(str(alignment_dir), "*")):
        if os.path.isfile(f):
            out.setdefault(os.path.basename(f).split(".")[0], f)
    return out


def _fmt_counts(counts: Counter) -> str:
    """"aa:n" pairs joined by commas, most frequent first (ties by residue)."""
    if not counts:
        return ""
    return ",".join(f"{aa}:{n}" for aa, n in
                    sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def add_species_tally(
    df: pd.DataFrame,
    alignment_dir,
    fg_species,
    bg_species,
    *,
    ali_format: str = "fasta",
    gene_col: str = "Gene",
    position_col: str = "Position",
) -> pd.DataFrame:
    """Add ``SPECIES_TALLY_COLUMNS`` to ``df``.

    ``Position`` is a 0-based alignment-column index. ``fg_species`` / ``bg_species`` are the
    full foreground / background contrast species, as lists or as files of
    ``selection/species_sets/{top,bottom}_species.txt``.

    Safe to call without inputs: when the alignment directory is missing, a species list is
    empty, the gene file is absent, the column is out of range or Bio.AlignIO is unavailable,
    the row keeps "" in these columns, so the output schema is stable.
    """
    out = df.copy()
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = ""

    fg = fg_species if isinstance(fg_species, (list, set, tuple)) else _read_species_list(fg_species)
    bg = bg_species if isinstance(bg_species, (list, set, tuple)) else _read_species_list(bg_species)
    fg_set, bg_set = set(fg), set(bg)
    if not fg_set or not bg_set or alignment_dir is None:
        return out
    if gene_col not in out.columns or position_col not in out.columns:
        return out

    try:
        from Bio import AlignIO
    except ImportError:
        return out

    files = _index_alignment_dir(alignment_dir)
    if not files:
        return out

    results: Dict[Tuple[str, object], Dict[str, str]] = {}
    for gene, g in out.groupby(gene_col, sort=False):
        path = files.get(str(gene))
        if path is None:
            continue
        try:
            aln = AlignIO.read(path, ali_format)
        except Exception:
            continue
        seqs = {rec.id: str(rec.seq) for rec in aln}
        aln_len = aln.get_alignment_length()
        fg_here = [seqs[s] for s in fg_set if s in seqs]
        bg_here = [seqs[s] for s in bg_set if s in seqs]
        if not fg_here or not bg_here:
            continue
        for pos in g[position_col].unique():
            try:
                idx = int(pos)
            except (TypeError, ValueError):
                continue
            if idx < 0 or idx >= aln_len:
                continue
            top_c: Counter = Counter()
            bot_c: Counter = Counter()
            for s in fg_here:
                aa = s[idx].upper()
                if aa not in _GAP_CHARS:
                    top_c[aa] += 1
            for s in bg_here:
                aa = s[idx].upper()
                if aa not in _GAP_CHARS:
                    bot_c[aa] += 1
            results[(gene, pos)] = {
                "top_species_residues": _fmt_counts(top_c),
                "bottom_species_residues": _fmt_counts(bot_c),
                "n_top_species": str(sum(top_c.values())),
                "n_bottom_species": str(sum(bot_c.values())),
            }

    if not results:
        return out
    keys = list(zip(out[gene_col], out[position_col]))
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = [results.get(k, _TALLY_EMPTY)[c] for k in keys]
    return out
