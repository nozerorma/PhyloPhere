# labelings.py — Readers of foreground/background labelings, trait files and PSS weights.
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
  * fop_pairs.tsv (null) and contrast_hypotheses_pairs.tsv (observed): the PSS weight of every
    (hypothesis, domain), `cycle` being absent from the observed file (its cycle is b_0).

Hypothesis ids are normalized to their `H<n>` token everywhere, so the observed and null sides
key the PSS weights identically.
"""
from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
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
    """The observed labelings: one per traitfile_H<n>.tab of a directory (`b_0~H<n>`), or the single trait
    file (`b_0`). Built on :func:`read_trait_pairs`, so it visits the hypotheses in the same order as the
    observed disambiguation and takes the same first species of each side per pair."""
    p = Path(path)
    out: Dict[str, Labeling] = {}
    for contrast, pairs in read_trait_pairs(p).items():
        if not pairs:
            continue
        tag = f"{cycle}~H{contrast}" if p.is_dir() else cycle
        out[tag] = Labeling(tag, tuple(f for f, _ in pairs), tuple(b for _, b in pairs))
    return out


def read_trait_pairs(
    trait_file_path: Path,
) -> Dict[int, List[Tuple[str, str]]]:
    """
    Parse trait file or directory of trait files and return species pairs grouped by contrast.

    Expected tab-separated format (only supported format):
    - No header
    - Exactly 3 columns per row: species, trait, pair
    - Returns {contrast: [(high_species, low_species), ...]} with pairs sorted by
      numeric pair_id where possible
    - Ignores rows with missing fields or invalid trait values
    - If a directory is provided, all *.tab files are parsed with contrast derived
      from filename (e.g. traitfile_H5.tab -> contrast 5) or sequential index.

    Contrasts come back in file-name order (H1, H10, H100, H11, ...), not numeric order: the order in which
    hypotheses are visited enters floating-point sums downstream, so it is part of the result.
    """
    files_to_read: List[Path] = []
    if trait_file_path.is_dir():
        # Specifically match H_n hypothesis files (traitfile_H1.tab, traitfile_H2.tab, ...)
        h_files = sorted(trait_file_path.glob("traitfile_H*.tab"))
        if h_files:
            files_to_read = h_files
        else:
            files_to_read = [
                f for f in sorted(trait_file_path.glob("*.tab"))
                if f.name != "traitfile_fop.tab"
            ]
        if not files_to_read:
            files_to_read = [f for f in sorted(trait_file_path.glob("*")) if f.is_file()]
    else:
        files_to_read = [trait_file_path]

    if not files_to_read:
        logger.warning("No trait files found in %s", trait_file_path)
        return {}

    import re

    # Structure: contrast -> pair_id -> {'high': [species], 'low': [species]}
    by_contrast_and_pair: Dict[int, Dict[str, Dict[str, list]]] = defaultdict(
        lambda: defaultdict(lambda: {"high": [], "low": []})
    )

    def _to_int(value: str) -> Optional[int]:
        try:
            return int(str(value).strip())
        except (ValueError, TypeError):
            return None

    try:
        for file_idx, fpath in enumerate(files_to_read, start=1):
            h_match = re.search(r"H(\d+)", fpath.name)
            default_contrast = int(h_match.group(1)) if h_match else file_idx

            with open(fpath, "r", encoding="utf-8-sig") as f:
                reader = csv.reader(f, delimiter="\t")
                rows_seen = 0

                def _consume_row(cells: List[str], row_num: int, contrast_num: int) -> None:
                    if not cells or all(not str(c).strip() for c in cells):
                        return

                    if len(cells) != 3:
                        logger.debug(
                            "Skipping malformed trait row %d in %s (expected 3 columns): %s",
                            row_num,
                            fpath.name,
                            cells,
                        )
                        return

                    species = cells[0].strip()
                    trait_val = _to_int(cells[1])
                    pair_id = cells[2].strip()

                    if not species or trait_val is None or not pair_id:
                        logger.debug(
                            "Skipping malformed trait row %d in %s: %s",
                            row_num,
                            fpath.name,
                            cells,
                        )
                        return

                    if trait_val == 1:
                        by_contrast_and_pair[contrast_num][pair_id]["high"].append(species)
                    elif trait_val == 0:
                        by_contrast_and_pair[contrast_num][pair_id]["low"].append(species)
                    else:
                        logger.debug(
                            "Skipping row %d in %s with non-binary trait value (%s): %s",
                            row_num,
                            fpath.name,
                            trait_val,
                            cells,
                        )

                for row_num, row in enumerate(reader, start=1):
                    rows_seen += 1
                    _consume_row(row, row_num, default_contrast)

        # Build final structure: contrast -> list of pairs
        contrast_to_pairs = {}
        for contrast_num, pairs_dict in by_contrast_and_pair.items():
            pairs_list = []
            # Sort by pair_id numerically to maintain order (pair 1, pair 2, pair 3)
            for pair_id in sorted(
                pairs_dict.keys(), key=lambda x: int(x) if x.isdigit() else x
            ):
                group = pairs_dict[pair_id]
                high_species = group["high"]
                low_species = group["low"]

                if high_species and low_species:
                    # Take first species from each side
                    pairs_list.append((high_species[0], low_species[0]))

                    if len(high_species) > 1 or len(low_species) > 1:
                        logger.warning(
                            "Contrast %d, pair %s has multiple species per side: "
                            "high=%s, low=%s. Using first from each.",
                            contrast_num,
                            pair_id,
                            high_species,
                            low_species,
                        )
                else:
                    logger.debug(
                        "Skipping incomplete pair %s in contrast %d: high=%s, low=%s",
                        pair_id,
                        contrast_num,
                        high_species,
                        low_species,
                    )

            contrast_to_pairs[contrast_num] = pairs_list
            logger.debug("Contrast %d: %d pairs loaded", contrast_num, len(pairs_list))

        total_pairs = sum(len(p) for p in contrast_to_pairs.values())
        logger.info(
            "Loaded %d contrasts with %d total pairs across %d traitfile(s) from %s",
            len(contrast_to_pairs),
            total_pairs,
            len(files_to_read),
            trait_file_path,
        )
        return contrast_to_pairs

    except Exception as e:
        logger.error("Failed to load trait pairs: %s", e, exc_info=True)
        return {}


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


def design_max_pairs(trait_file: Path) -> int:
    """Largest pair id of the observed design (a trait file, or a directory of traitfile_H*.tab): the number of
    pair columns the master CSV schema needs. 1 when there is none."""
    max_pair_id = 0
    files_to_read = []
    if trait_file.is_dir():
        h_files = sorted(trait_file.glob("traitfile_H*.tab"))
        if h_files:
            files_to_read = h_files
        else:
            files_to_read = [
                f for f in sorted(trait_file.glob("*.tab"))
                if f.name != "traitfile_fop.tab"
            ]
        if not files_to_read:
            files_to_read = [f for f in sorted(trait_file.glob("*")) if f.is_file()]
    elif trait_file.is_file():
        files_to_read = [trait_file]

    for fpath in files_to_read:
        try:
            with open(fpath, "r", encoding="utf-8-sig") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if not row or len(row) < 3:
                        continue
                    pair_val = row[2]
                    if pair_val is None or str(pair_val).strip() == "":
                        continue
                    try:
                        pair_id = int(str(pair_val).strip())
                        max_pair_id = max(max_pair_id, pair_id)
                    except Exception:
                        continue
        except Exception as e:
            logger.warning(f"Could not read traitfile {fpath}: {e}")
    return max_pair_id or 1


def observed_pss(path) -> Optional[Dict[Tuple[str, int], float]]:
    """PSS weights of the observed cycle, or None when there are none."""
    return read_pss(path).get(OBSERVED_CYCLE) or None
