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

import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION/local"))
from src.core.ct_io import strict_io
from src.utils.io_utils import _scan_alignment_dir, read_alignment

# ── Extant-species residue tally (alignment-based) ────────────────────────────
# Counts, per (Gene, Position), the contrast species of the hypotheses that call the position
# (the `participating_hypotheses` of its rows, unioned over sides) that carry each
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

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _read_hyp_pairs(path) -> Dict[str, List[Tuple[str, str]]]:
    """Read complete independent pairs and validate the shared CT design contract."""
    if path is None:
        return {}
    pairs = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    required = {'hypothesis_id', 'pair', 'species1', 'species2'}
    if not required <= set(pairs):
        raise ValueError(f'Missing paired-design columns: {path}')
    out = {}
    for h, frame in pairs.groupby('hypothesis_id', sort=False):
        if not h.startswith('H') or not h[1:].isdigit() or int(h[1:]) < 1:
            raise ValueError(f'Invalid hypothesis: {h}')
        design, seen = {}, set()
        for row in frame.to_dict('records'):
            pair = strict_io.integer(row['pair'], 'pair', 1)
            a, b = row['species1'], row['species2']
            if (pair in design or a == b or any(not x or x != x.strip() or ',' in x or any(c.isspace() for c in x) for x in [a,b])
                    or a in seen or b in seen):
                raise ValueError(f'Invalid or duplicate contrast pair in {h}')
            design[pair] = (a, b)
            seen.update([a,b])
            if 'pss_score' in row:
                import math
                weight = float(row['pss_score'])
                if not math.isfinite(weight) or weight < 0:
                    raise ValueError('Invalid PSS weight')
        if set(design) != set(range(1, len(design)+1)):
            raise ValueError('Pair IDs must be contiguous')
        out[h] = [design[i] for i in sorted(design)]
    if not out or {int(h[1:]) for h in out} != set(range(1,len(out)+1)) or len({len(v) for v in out.values()}) != 1:
        raise ValueError('Empty, incomplete or unequal-size hypothesis design')
    return out


def _index_alignment_dir(alignment_dir):
    directory = Path(alignment_dir)
    if not directory.is_dir():
        raise ValueError(f'Alignment directory not found: {directory}')
    return _scan_alignment_dir.__wrapped__(directory)


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
    ns = [n for _, n in per_hyp]
    if len(set(ns)) > 1:
        return ";".join(f"{k}/{n}" for k, n in per_hyp)
    n = ns[0]
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
            if ra not in _VALID_AA or rb not in _VALID_AA:
                continue
            n += 1
            k += ra == rb
        if n:
            per_hyp.append((k, n))
    top_c = Counter(seqs[s][idx].upper() for s in top_sp if s in seqs and seqs[s][idx].upper() in _VALID_AA)
    bot_c = Counter(seqs[s][idx].upper() for s in bot_sp if s in seqs and seqs[s][idx].upper() in _VALID_AA)
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
    taxid_mapping=None,
) -> pd.DataFrame:
    """Add ``SPECIES_TALLY_COLUMNS`` to ``df``.

    ``Position`` is a 0-based alignment-column index. The hypotheses of a position are the union
    of ``participating_hypotheses`` over its rows (all hypotheses of the design when the column
    is absent).

    Optional absent inputs leave empty descriptors. Supplied invalid designs,
    missing/ambiguous alignments, missing species and out-of-range positions fail;
    unknown residues do not contribute to species or conserved-pair counts.
    """
    out = df.copy()
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = ""

    hyp_pairs = _read_hyp_pairs(hyp_pairs_file)
    if not hyp_pairs or alignment_dir is None:
        return out
    if gene_col not in out.columns or position_col not in out.columns:
        return out

    mapping = {}
    if taxid_mapping:
        from src.phylo.species_mapping import read_taxid_mapping
        mapping = read_taxid_mapping(Path(taxid_mapping))
        if len(set(mapping.values())) != len(mapping):
            raise ValueError('Ambiguous taxon mapping')
    files = _index_alignment_dir(alignment_dir)
    if not files:
        return out

    # Hypotheses of each (Gene, Position): union over its rows, in design order.
    order = {h: i for i, h in enumerate(hyp_pairs)}
    by_gene: Dict[str, Dict[object, set]] = {}
    for gene, pos, hs in zip(out[gene_col], out[position_col],
                             out[hyp_col] if hyp_col in out.columns else [""] * len(out)):
        found = {x.strip() for x in str(hs).split(',') if x.strip() and x.strip() != 'nan'}
        if found - hyp_pairs.keys():
            raise ValueError(f'Unknown participating hypotheses: {found - hyp_pairs.keys()}')
        # Empty provenance on diagnostic rows means all design hypotheses.
        by_gene.setdefault(gene, {}).setdefault(pos, set()).update(found or hyp_pairs)

    results: Dict[Tuple[str, object], Dict[str, str]] = {}
    for gene, positions in by_gene.items():
        path = files.get(str(gene))
        if path is None:
            raise FileNotFoundError(f'Alignment missing for {gene}')
        aln = read_alignment(path, ali_format)
        seqs = {rec.id: str(rec.seq) for rec in aln}
        if len(seqs) != len(aln):
            raise ValueError(f'Duplicate alignment species: {gene}')
        required_species = {sp for hs in positions.values() for h in hs for pair in hyp_pairs[h] for sp in pair}
        for species, taxid in mapping.items():
            if str(taxid) in seqs:
                if species in seqs and seqs[species] != seqs[str(taxid)]:
                    raise ValueError(f'Conflicting species mapping: {species}')
                seqs[species] = seqs[str(taxid)]
        if required_species - seqs.keys():
            raise ValueError(f'Contrast species missing from alignment {gene}: {required_species - seqs.keys()}')
        aln_len = aln.get_alignment_length()
        for pos, hs in positions.items():
            idx = int(pos)
            if not 0 <= idx < aln_len:
                raise ValueError(f'Position outside alignment: {gene}:{idx}')
            results[(gene, pos)] = tally_position(seqs, idx, sorted(hs, key=order.get), hyp_pairs)

    keys = list(zip(out[gene_col], out[position_col]))
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = [results.get(k, _TALLY_EMPTY)[c] for k in keys]
    return out
