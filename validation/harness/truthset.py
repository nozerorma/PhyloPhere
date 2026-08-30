"""Truth-set schema + loaders.

Two kinds of truth set, both TSV with a header and ``#``-comment lines:

site-level  (truthsets/tier1/*.sites.tsv)
    gene            reference gene symbol / ortholog id used in the fixture
    ref_seq         sequence whose numbering ``position`` refers to
    position        1-based residue index in ``ref_seq`` (ungapped)
    ref_aa,alt_aa   ancestral / derived residue where known ("." if unknown)
    tier            evidence class: mutagenesis | manual-alignment | reanalysis
    source          short citation key (full ref in the sibling README)

gene-level  (truthsets/tier1/*.genes.tsv, tier2/tier3 *.genes.tsv)
    gene            reference gene symbol / ortholog id
    direction       accelerated | fg-shift | either
    tier            evidence class
    source          short citation key

``position`` in a site truth set is in ``ref_seq`` coordinates; map it into the
fixture alignment's column space with ``map_site_to_alignment`` before comparing
to pipeline output (which reports alignment columns / gene positions).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Site:
    gene: str
    ref_seq: str
    position: int
    ref_aa: str
    alt_aa: str
    tier: str
    source: str

    @property
    def key(self) -> tuple[str, int]:
        return (self.gene, self.position)


@dataclass(frozen=True)
class GeneCall:
    gene: str
    direction: str
    tier: str
    source: str


def _rows(path: Path):
    with path.open() as fh:
        reader = csv.DictReader(
            (ln for ln in fh if not ln.lstrip().startswith("#")), delimiter="\t"
        )
        for row in reader:
            if any((v or "").strip() for v in row.values()):
                yield {k: (v or "").strip() for k, v in row.items()}


def load_sites(path: str | Path) -> list[Site]:
    return [
        Site(
            gene=r["gene"], ref_seq=r["ref_seq"], position=int(r["position"]),
            ref_aa=r.get("ref_aa", "."), alt_aa=r.get("alt_aa", "."),
            tier=r.get("tier", ""), source=r.get("source", ""),
        )
        for r in _rows(Path(path))
    ]


def load_genes(path: str | Path) -> list[GeneCall]:
    return [
        GeneCall(
            gene=r["gene"], direction=r.get("direction", "either"),
            tier=r.get("tier", ""), source=r.get("source", ""),
        )
        for r in _rows(Path(path))
    ]


def map_site_to_alignment(position: int, ref_row: str) -> int:
    """1-based ungapped ``position`` in ``ref_row`` -> 1-based alignment column.

    ``ref_row`` is the gapped aligned sequence of ``Site.ref_seq`` taken from the
    fixture alignment.
    """
    seen = 0
    for col, ch in enumerate(ref_row, start=1):
        if ch not in "-.":
            seen += 1
            if seen == position:
                return col
    raise ValueError(f"position {position} beyond ungapped length {seen} of ref row")
