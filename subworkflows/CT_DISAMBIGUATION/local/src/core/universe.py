# universe.py — The positional universe of the observed design, as the replay of a tile and the contract read it.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
`<gene>\\t<position,position,...|NULL>` lines written by `ct perm-replay --export_b0_universe` (one `*.b0.universe` file
per gene): the positions b_0 tested that pass the CAAS test under some orientation of the pairs of some observed
hypothesis (subworkflows/CT/local/modules/universe.py). Positions are zero based columns of the trimmed alignment.

Imported by: src/utils/gene_wrapper.py (the null keeps only these positions), src/core/contract.py (background.output)
"""
import hashlib
from pathlib import Path
from typing import Dict, Iterable, Tuple

UNIVERSE_SUFFIX = ".b0.universe"
STRICT_SUFFIX = ".b0.universe_strict"   # same lines, judged with no conserved pair allowed (the strict analysis)


def parse_line(line: str) -> Tuple[str, Tuple[int, ...]]:
    """(gene, ascending positions) of one universe line; NULL or an empty field is no position."""
    gene, tab, text = line.rstrip("\n").partition("\t")
    if not tab or not gene or gene.strip() != gene:
        raise ValueError(f"Invalid universe line: {line!r}")
    values = [] if text in ("", "NULL") else text.split(",")
    if any(not v.isdigit() or int(v) >= 2 ** 31 for v in values):
        raise ValueError(f"Invalid zero-based universe position for {gene}")
    positions = tuple(sorted(int(v) for v in values))
    if len(set(positions)) != len(positions):
        raise ValueError(f"Duplicate universe position for {gene}")
    return gene, positions


def read_files(files: Iterable[Path]) -> Dict[str, Tuple[int, ...]]:
    """{gene: positions} of universe files; a gene in two files is an error."""
    universe: Dict[str, Tuple[int, ...]] = {}
    for path in files:
        for line in Path(path).read_text().splitlines():
            if not line.strip():
                continue
            gene, positions = parse_line(line)
            if gene in universe:
                raise ValueError(f"Duplicate gene in the positional universe: {gene}")
            universe[gene] = positions
    return universe


def read_dir(directory: Path, suffix: str = UNIVERSE_SUFFIX) -> Dict[str, Tuple[int, ...]]:
    """The universe of every `*<suffix>` file under a directory."""
    files = sorted(Path(directory).glob(f"*{suffix}"))
    if not files:
        raise ValueError(f"No positional universe files (*{suffix}) in {directory}")
    return read_files(files)


def digest(gene: str, positions: Iterable[int]) -> str:
    """SHA-256 of the canonical line of a gene (ascending positions, NULL when empty): the identity of its universe
    in the tile provenance and in background.output (SCORING recomputes it from the file it is given)."""
    ordered = sorted(positions)
    return hashlib.sha256(f"{gene}\t{','.join(map(str, ordered)) if ordered else 'NULL'}".encode()).hexdigest()
