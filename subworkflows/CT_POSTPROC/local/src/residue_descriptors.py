#!/usr/bin/env python3
# residue_descriptors.py — Residues of the extant contrast species at each CAAS position.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Extant-species residue tally for the disambiguation table: for each (Gene, Position), how
many foreground (top) and background (bottom) contrast species carry each residue, and how
many contrast pairs keep the same residue on both sides.

Imported by: prepare_postproc_input.py (add_species_tally), which runs in the input-prep
             step (CAAS_PREPARE_POSTPROC_INPUT), before the cluster and gene filters
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd

# ── Extant-species residue tally (alignment-based) ────────────────────────────
# Counts, per (Gene, Position), the contrast species of the hypotheses that call the position
# (the `participating_hypotheses` of its rows, unioned over schemes and sides) that carry each
# residue at the alignment column (an "X of N species" check on the call). The species come from
# the contrast design (contrast_hypotheses_pairs.tsv: species1 = top / foreground, species2 =
# bottom / background of each pair), so the tally needs neither FADE nor the species_sets of
# the selection step.

SPECIES_TALLY_COLUMNS = (
    "top_species_residues",     # e.g. "N:15,A:8" -- FG contrast species, count-desc
    "bottom_species_residues",  # e.g. "A:52,S:2" -- BG contrast species
    "n_top_species",            # FG contrast species with a non-gap residue here
    "n_bottom_species",         # BG contrast species with a non-gap residue here
    "n_conserved_pairs",        # e.g. "3/6": pairs with the same residue on both sides / pairs counted
)

_TALLY_EMPTY = {c: "" for c in SPECIES_TALLY_COLUMNS}

_GAP_CHARS = set("-.")


def _read_hyp_pairs(path) -> Dict[str, List[Tuple[str, str]]]:
    """{hypothesis: [(top_species, bottom_species), ...]} from contrast_hypotheses_pairs.tsv.

    {} when the file is absent, unreadable or lacks hypothesis_id / species1 / species2.
    """
    if not path:
        return {}
    try:
        pairs = pd.read_csv(path, sep="\t", dtype=str)
    except (OSError, ValueError, pd.errors.EmptyDataError):
        return {}
    if not {"hypothesis_id", "species1", "species2"} <= set(pairs.columns):
        return {}
    out: Dict[str, List[Tuple[str, str]]] = {}
    for h, a, b in zip(pairs["hypothesis_id"], pairs["species1"], pairs["species2"]):
        out.setdefault(str(h).strip(), []).append((str(a).strip(), str(b).strip()))
    return out


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


def _fmt_conserved(per_hyp: List[Tuple[int, int]]) -> str:
    """"k/n" over the hypotheses of a position; "k1-k2/n" when the count differs between them."""
    if not per_hyp:
        return ""
    ks = [k for k, _ in per_hyp]
    n = max(n for _, n in per_hyp)
    return f"{min(ks)}/{n}" if min(ks) == max(ks) else f"{min(ks)}-{max(ks)}/{n}"


def tally_position(seqs: Dict[str, str], idx: int, hyps: List[str],
                   hyp_pairs: Dict[str, List[Tuple[str, str]]]) -> Dict[str, str]:
    """Tally of one alignment column over the contrast pairs of `hyps`.

    Top / bottom counts are over the distinct species of those hypotheses. A pair is conserved
    when both of its species have a non-gap residue here and it is the same on both sides; the
    count is taken per hypothesis (a hypothesis is a set of independent pairs).
    """
    top_sp, bot_sp, per_hyp = set(), set(), []
    for h in hyps:
        k = n = 0
        for a, b in hyp_pairs.get(h, []):
            top_sp.add(a)
            bot_sp.add(b)
            ra = seqs[a][idx].upper() if a in seqs else "-"
            rb = seqs[b][idx].upper() if b in seqs else "-"
            if ra in _GAP_CHARS or rb in _GAP_CHARS:
                continue
            n += 1
            k += ra == rb
        if n:
            per_hyp.append((k, n))
    top_c = Counter(seqs[s][idx].upper() for s in top_sp if s in seqs and seqs[s][idx] not in _GAP_CHARS)
    bot_c = Counter(seqs[s][idx].upper() for s in bot_sp if s in seqs and seqs[s][idx] not in _GAP_CHARS)
    return {
        "top_species_residues": _fmt_counts(top_c),
        "bottom_species_residues": _fmt_counts(bot_c),
        "n_top_species": str(sum(top_c.values())),
        "n_bottom_species": str(sum(bot_c.values())),
        "n_conserved_pairs": _fmt_conserved(per_hyp),
    }


def add_species_tally(
    df: pd.DataFrame,
    alignment_dir,
    hyp_pairs_file,
    *,
    ali_format: str = "fasta",
    gene_col: str = "Gene",
    position_col: str = "Position",
    hyp_col: str = "participating_hypotheses",
) -> pd.DataFrame:
    """Add ``SPECIES_TALLY_COLUMNS`` to ``df``.

    ``Position`` is a 0-based alignment-column index. The hypotheses of a position are the union
    of ``participating_hypotheses`` over its rows (all hypotheses of the design when the column
    is absent).

    Safe to call without inputs: when the alignment directory or the pairs file is missing, the
    gene file is absent, the column is out of range or Bio.AlignIO is unavailable, the row keeps
    "" in these columns, so the output schema is stable.
    """
    out = df.copy()
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = ""

    hyp_pairs = _read_hyp_pairs(hyp_pairs_file)
    if not hyp_pairs or alignment_dir is None:
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

    # Hypotheses of each (Gene, Position): union over its rows, in design order.
    order = {h: i for i, h in enumerate(hyp_pairs)}
    by_gene: Dict[str, Dict[object, set]] = {}
    for gene, pos, hs in zip(out[gene_col], out[position_col],
                             out[hyp_col] if hyp_col in out.columns else [""] * len(out)):
        found = {x.strip() for x in str(hs).split(",") if x.strip() in hyp_pairs}
        by_gene.setdefault(gene, {}).setdefault(pos, set()).update(found or hyp_pairs)

    results: Dict[Tuple[str, object], Dict[str, str]] = {}
    for gene, positions in by_gene.items():
        path = files.get(str(gene))
        if path is None:
            continue
        try:
            aln = AlignIO.read(path, ali_format)
        except Exception:
            continue
        seqs = {rec.id: str(rec.seq) for rec in aln}
        aln_len = aln.get_alignment_length()
        for pos, hs in positions.items():
            try:
                idx = int(pos)
            except (TypeError, ValueError):
                continue
            if 0 <= idx < aln_len:
                results[(gene, pos)] = tally_position(seqs, idx, sorted(hs, key=order.get), hyp_pairs)

    keys = list(zip(out[gene_col], out[position_col]))
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = [results.get(k, _TALLY_EMPTY)[c] for k in keys]
    return out
