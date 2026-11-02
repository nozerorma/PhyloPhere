#!/usr/bin/env python3
# residue_descriptors.py — Residues of the extant contrast species at each CAAS position.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Extant-species residue tally for each (Gene, Position, hypothesis): how
many foreground (top) and background (bottom) contrast species carry each residue, and how
many contrast pairs keep the same residue on both sides.

Imported by: prepare_postproc_input.py (add_species_tally), which runs in the input-prep
             step (CAAS_PREPARE_POSTPROC_INPUT), before the cluster and gene filters
"""

from __future__ import annotations

from collections import Counter
import logging
from typing import Dict, List, Tuple

import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION/local"))
from src.core.ct_io import strict_io
from src.utils.io_utils import _scan_alignment_dir, read_alignment

logger = logging.getLogger(__name__)

# ── Extant-species residue tally (alignment-based) ────────────────────────────
# Each descriptor cell contains hypothesis-labelled values, e.g. H1=A:2;H2=T:2.
# Only the row's participating hypotheses are described; species are never pooled
# across hypotheses or across rows/sides. Empty provenance leaves empty descriptors.

SPECIES_TALLY_COLUMNS = (
    "top_species_residues",     # e.g. "N:15,A:8" -- FG contrast species, count-desc
    "bottom_species_residues",  # e.g. "A:52,S:2" -- BG contrast species
    "n_top_species",            # FG contrast species with a non-gap residue here
    "n_bottom_species",         # BG contrast species with a non-gap residue here
    "n_conserved_pairs",        # e.g. "3/6": pairs with the same residue on both sides / pairs counted
)

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


def tally_position(seqs: Dict[str, str], idx: int, hypothesis: str,
                   hyp_pairs: Dict[str, List[Tuple[str, str]]]) -> Dict[str, str]:
    """Tally one hypothesis at one alignment column, using available valid residues.

    A conserved pair requires valid amino acids in both members. An available member
    of an incomplete pair still contributes to its side's residue count.
    """
    top_c, bot_c = Counter(), Counter()
    k = n = 0
    for a, b in hyp_pairs[hypothesis]:
        ra = seqs[a][idx].upper() if a in seqs else "-"
        rb = seqs[b][idx].upper() if b in seqs else "-"
        if ra in _VALID_AA:
            top_c[ra] += 1
        if rb in _VALID_AA:
            bot_c[rb] += 1
        if ra in _VALID_AA and rb in _VALID_AA:
            n += 1
            k += ra == rb
    return {
        "top_species_residues": _fmt_counts(top_c),
        "bottom_species_residues": _fmt_counts(bot_c),
        "n_top_species": str(sum(top_c.values())),
        "n_bottom_species": str(sum(bot_c.values())),
        "n_conserved_pairs": f"{k}/{n}",
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

    ``Position`` is a 0-based alignment-column index. Each cell contains semicolon-
    separated ``Hn=value`` entries for that row's participating hypotheses, in design
    order. Missing/empty provenance leaves descriptors empty; no design-wide fallback
    is inferred. Counts are computed independently per hypothesis.

    Optional absent inputs leave empty descriptors. Supplied invalid designs,
    missing/ambiguous alignments and out-of-range positions fail. Missing species
    are warned about and, like unknown residues, do not contribute to species or
    conserved-pair counts; available members of incomplete pairs still contribute
    to their side's species tally.
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

    # Cache computations per (Gene, Position, hypothesis); retain each row's provenance.
    order = {h: i for i, h in enumerate(hyp_pairs)}
    by_gene: Dict[str, Dict[object, set]] = {}
    row_hypotheses = []
    for gene, pos, hs in zip(out[gene_col], out[position_col],
                             out[hyp_col] if hyp_col in out.columns else [""] * len(out)):
        found = {x.strip() for x in str(hs).split(',') if x.strip() and x.strip() != 'nan'}
        if found - hyp_pairs.keys():
            raise ValueError(f'Unknown participating hypotheses: {found - hyp_pairs.keys()}')
        row_hypotheses.append(sorted(found, key=order.get))
        if found:
            by_gene.setdefault(gene, {}).setdefault(pos, set()).update(found)

    results: Dict[Tuple[str, object, str], Dict[str, str]] = {}
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
        missing_species = required_species - seqs.keys()
        if missing_species:
            logger.warning(
                'Contrast species missing from alignment %s: %s; '
                'tallying available residues and excluding incomplete pairs',
                gene, ', '.join(sorted(missing_species)),
            )
        aln_len = aln.get_alignment_length()
        for pos, hs in positions.items():
            idx = int(pos)
            if not 0 <= idx < aln_len:
                raise ValueError(f'Position outside alignment: {gene}:{idx}')
            for h in sorted(hs, key=order.get):
                results[(gene, pos, h)] = tally_position(seqs, idx, h, hyp_pairs)

    keys = list(zip(out[gene_col], out[position_col]))
    for c in SPECIES_TALLY_COLUMNS:
        out[c] = [";".join(f"{h}={results[(gene, pos, h)][c]}" for h in hs)
                  for (gene, pos), hs in zip(keys, row_hypotheses)]
    return out
