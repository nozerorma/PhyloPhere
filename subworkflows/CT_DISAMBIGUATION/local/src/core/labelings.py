# labelings.py — Readers of foreground/background labelings, trait files and the design pairs with their PSS.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Labelings and their hypothesis weights: one reader for the observed design and the null.

Imported by: contract_main.py, observed_b0_main.py, explain_positions.py, src/core/driver.py, src/core/meta.py,
src/utils/gene_wrapper.py
Inputs: resample_*.tab, fop_labelings.tab, traitfile*.tab, fop_pairs.tsv, contrast_hypotheses_pairs.tsv
Outputs: in-memory Labeling objects, trait pairs and PSS maps

A labeling is one foreground/background assignment: a resample cycle (`b_12`), or one of the
Dunn-independent hypotheses of a cycle (`b_12~H3`, the FOP fan-out). The real (observed) labeling
is the cycle `b_0`, with its hypotheses read from the run's own traitfile_H*.tab.

Three file formats describe them; they are read here and nowhere else:

  * resample_*.tab / fop_labelings.tab   `tag <TAB> fg_csv <TAB> bg_csv`, no header
    (fg[k] and bg[k] are pair k+1, as in a trait file);
  * traitfile*.tab                       `species <TAB> 1|0 <TAB> pair` (1 = foreground);
  * fop_pairs.tsv (null) and contrast_hypotheses_pairs.tsv (observed): the species and the PSS of the pair
    of every (hypothesis, domain), `cycle` being absent from the observed file (its cycle is b_0).

Hypothesis ids are normalized to their `H<n>` token everywhere, so the observed and null sides
key the design pairs identically.
"""
from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

OBSERVED_CYCLE = "b_0"

_HYP = re.compile(r"H\d+")


def hyp_id(value: str) -> str:
    """'b_12~H3' / 'traitfile_H3.tab' / 'H3' -> 'H3'; no hypothesis token -> 'H1' (single contrast)."""
    m = _HYP.search(str(value))
    return m.group(0) if m else "H1"


def base_cycle(tag: str) -> str:
    """'b_12~H3' -> 'b_12'; a plain 'b_12' passes through."""
    return str(tag).split("~", 1)[0]


def trait_pairs_from(fg, bg) -> Dict[int, List[Tuple[str, str]]]:
    """The single-contrast trait_pairs the disambiguation consumes: pair k is (fg[k], bg[k]).

    The contrast key is never inspected downstream when there is exactly one contrast.
    """
    if not fg or len(fg) != len(bg) or len(set(fg) | set(bg)) != len(fg) + len(bg):
        raise ValueError('Invalid paired labeling')
    return {1: list(zip(fg, bg))}


def design_pairs_from(trait_pairs: Dict[int, List[Tuple[str, str]]]) -> Dict[Tuple[str, int], Tuple[str, str, Optional[float]]]:
    """The design of `trait_pairs` as `{(H<contrast>, domain): (top species, bottom species, None)}`: the pairs of a design that has
    no fop_pairs / contrast_hypotheses_pairs table, so without PSS (every pair weighs the same). Domains are numbered from 1."""
    return {(f"H{c}", i): (top, bottom, None) for c, pairs in trait_pairs.items() for i, (top, bottom) in enumerate(pairs, 1)}


@dataclass(frozen=True)
class Labeling:
    tag: str
    fg: Tuple[str, ...]
    bg: Tuple[str, ...]

    @property
    def base(self) -> str:
        return base_cycle(self.tag)

    @property
    def hyp(self) -> str:
        return hyp_id(self.tag.split("~", 1)[1]) if "~" in self.tag else "H1"

    def trait_pairs(self) -> Dict[int, List[Tuple[str, str]]]:
        """The single-contrast shape the disambiguation consumes: pair k is (fg[k], bg[k])."""
        return trait_pairs_from(self.fg, self.bg)


def read_labelings(path) -> Dict[str, Labeling]:
    """Every labeling of a resample directory (or file), keyed by tag.

    In a directory `fop_labelings.tab` (tags `<base>~H<m>`) takes precedence over `resample_*.tab`.
    """
    p = Path(path)
    if p.is_dir():
        tabs = sorted(p.glob("fop_labelings.tab")) or sorted(p.glob("resample_*.tab"))
    else:
        tabs = [p]
    from src.core.ct_io import strict_io
    return {tag: Labeling(tag, fg, bg) for tag, (fg, bg) in strict_io.read_labelings(tabs).items()}


def read_trait_pairs(trait_file_path: Path) -> Dict[int, List[Tuple[str, str]]]:
    """Read the exact paired design using the same validation as CT preparation."""
    from src.core.ct_io import strict_io
    design = strict_io.read_design(trait_file_path)
    return {int(h[1:]): [(pair['1'], pair['0']) for pair in pairs.values()]
            for h, pairs in design.items()}


PssMap = Dict[str, Dict[Tuple[str, int], Tuple[str, str, float]]]


def read_pss(path, default_cycle: str = OBSERVED_CYCLE, labelings=None) -> PssMap:
    """Read finite, nonnegative PSS values without repairing supplied metadata.

    With labelings, enforce exact active coverage and species orientation using
    CT's canonical metadata validator. Metadata for inactive cycles is ignored.
    """
    import math
    from src.core.ct_io import strict_io
    if not path or Path(path).name.startswith('NO_'):
        return {}
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Pair metadata not found: {path}')
    out: PssMap = {}
    with path.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'hypothesis_id', 'pair', 'species1', 'species2', 'pss_score'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f'Missing pair metadata columns in {path}')
        observed = 'cycle' not in reader.fieldnames
        for row in reader:
            h = row['hypothesis_id']
            if not re.fullmatch(r'H[1-9][0-9]*', h):
                raise ValueError(f'Invalid hypothesis ID: {h!r}')
            domain = strict_io.integer(row['pair'], 'metadata pair', 1)
            weight = float(row['pss_score'])
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f'Invalid PSS weight in {path}: {row}')
            cycle = row.get('cycle') or default_cycle
            if not re.fullmatch(r'b_(0|[1-9][0-9]*)', cycle):
                raise ValueError(f'Invalid cycle: {cycle!r}')
            key = (h, domain)
            if key in out.setdefault(cycle, {}):
                raise ValueError(f'Duplicate pair metadata: {cycle}, {key}')
            if not row['species1'] or not row['species2']:
                raise ValueError(f'Missing pair species in {path}')
            out[cycle][key] = (row['species1'], row['species2'], weight)
    if labelings is not None:
        validated = strict_io.pair_metadata(path, labelings, observed=observed)
        out = {}
        for (cycle, hyp, domain), weight in validated.items():
            tag = f'{cycle}~{hyp}'
            fg, bg = labelings.get(tag, labelings.get(cycle, ((), ())))
            out.setdefault(cycle, {})[(hyp, domain)] = (fg[domain-1], bg[domain-1], weight)
    return out


def design_max_pairs(trait_file: Path) -> int:
    """The common, validated number of domains of the observed design."""
    return len(next(iter(read_trait_pairs(trait_file).values())))


def observed_pss(path, trait_pairs=None):
    labelings = None if trait_pairs is None else {
        f'b_0~H{h}': (tuple(f for f, _ in pairs), tuple(b for _, b in pairs))
        for h, pairs in trait_pairs.items()}
    return read_pss(path, labelings=labelings).get(OBSERVED_CYCLE) or None
