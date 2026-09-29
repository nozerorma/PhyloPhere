"""Labelings and their hypothesis weights: one reader for the observed design and the null.

A labeling is one foreground/background assignment: a resample cycle (`b_12`), or one of the
Dunn-independent hypotheses of a cycle (`b_12~H3`, the FOP fan-out). The real (observed) labeling
is the cycle `b_0`, with its hypotheses read from the run's own traitfile_H*.tab.

Three file formats describe them; they are read here and nowhere else:

  * resample_*.tab / fop_labelings.tab   `tag <TAB> fg_csv <TAB> bg_csv`, no header
    (fg[k] and bg[k] are pair k+1, as in a trait file);
  * traitfile*.tab                       `species <TAB> 1|0 <TAB> pair` (1 = foreground);
  * fop_pairs.tsv (null) and contrast_hypotheses_pairs.tsv (observed): the PSS weight of every
    (hypothesis, domain), `cycle` being absent from the observed file (its cycle is b_0).

Hypothesis ids are normalized to their `H<n>` token everywhere, so the observed and null sides
key the PSS weights identically.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    return {1: list(zip(fg, bg))}


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
    out: Dict[str, Labeling] = {}
    for tab in tabs:
        with open(tab, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) < 3:
                    continue
                tag = row[0].strip()
                fg = tuple(s for s in row[1].split(",") if s.strip())
                bg = tuple(s for s in row[2].split(",") if s.strip())
                if tag and fg and bg:
                    out[tag] = Labeling(tag, fg, bg)
    return out


def read_design(path, cycle: str = OBSERVED_CYCLE) -> Dict[str, Labeling]:
    """The observed labelings: one per traitfile_H<n>.tab of a directory, or the single trait file.

    A directory yields `b_0~H<n>` for each hypothesis file; a single file yields the plain `b_0`.
    Pairs are ordered by pair id and take the first species of each side, as the disambiguation's
    trait-file parser does.
    """
    p = Path(path)
    if p.is_dir():
        files = sorted(p.glob("traitfile_H*.tab"), key=lambda f: int(hyp_id(f.name)[1:]))
        tagged = [(f, f"{cycle}~{hyp_id(f.name)}") for f in files]
    else:
        tagged = [(p, cycle)]
    out: Dict[str, Labeling] = {}
    for f, tag in tagged:
        fg: Dict[int, str] = {}
        bg: Dict[int, str] = {}
        with open(f, newline="", encoding="utf-8-sig") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) != 3 or not row[2].strip().isdigit():
                    continue
                side = fg if row[1].strip() == "1" else bg if row[1].strip() == "0" else None
                if side is not None:
                    side.setdefault(int(row[2]), row[0].strip())
        pairs = sorted(set(fg) & set(bg))
        if pairs:
            out[tag] = Labeling(tag, tuple(fg[k] for k in pairs), tuple(bg[k] for k in pairs))
    return out


PssMap = Dict[str, Dict[Tuple[str, int], float]]


def read_pss(path, default_cycle: str = OBSERVED_CYCLE) -> PssMap:
    """PSS weights `{cycle: {(H<n>, domain): pss}}` from fop_pairs.tsv or contrast_hypotheses_pairs.tsv.

    Rows lacking a cycle column (the observed file) belong to `default_cycle`. Missing files, `NO_*`
    sentinels and files without the needed columns give an empty map (the pooler then weights every
    node equally).
    """
    if not path or str(path).startswith("NO_") or not Path(path).is_file():
        return {}
    out: PssMap = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"hypothesis_id", "pair", "pss_score"} <= set(reader.fieldnames):
            return {}
        for row in reader:
            raw = (row.get("hypothesis_id") or "").strip()
            if not _HYP.search(raw):
                continue
            try:
                domain = int(float(row["pair"]))
                pss = float(row["pss_score"])
            except (TypeError, ValueError):
                continue
            cycle = (row.get("cycle") or "").strip() or default_cycle
            out.setdefault(cycle, {})[(hyp_id(raw), domain)] = pss
    return out


def observed_pss(path) -> Optional[Dict[Tuple[str, int], float]]:
    """PSS weights of the observed cycle, or None when there are none."""
    return read_pss(path).get(OBSERVED_CYCLE) or None
