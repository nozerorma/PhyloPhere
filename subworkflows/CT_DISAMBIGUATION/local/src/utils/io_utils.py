# io_utils.py — Locate and read the per-gene alignment files.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
I/O utilities: find the alignment of a gene in a directory tree and read it with Biopython.

Imported by: src/core/driver.py (find_gene_alignment), src/asr/asr_single.py (read_alignment)
"""

import functools
from pathlib import Path
from typing import Dict, Optional, Set
import logging

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

logger = logging.getLogger(__name__)


def _is_supported_alignment_path(path: Path) -> bool:
    """True for a regular file whose extension is a PHYLIP or FASTA one, or that has none."""
    suffix = path.suffix.lower()
    return suffix in {".phy", ".phylip", ".aln", ".fa", ".fasta", ".fas", ".clustal", ""} and path.is_file()


@functools.lru_cache(maxsize=32)
def _scan_alignment_dir(alignment_dir: Path) -> Dict[str, Path]:
    """Index complete IDs, stripping only the alignment suffix; reject ambiguity."""
    result = {}
    for path in sorted(alignment_dir.glob('**/*')):
        if not _is_supported_alignment_path(path):
            continue
        gene = path.stem if path.suffix else path.name
        if gene in result:
            raise ValueError(f'Ambiguous alignment for {gene}: {result[gene]}, {path}')
        result[gene] = path
    return result


def read_alignment_manifest(path, alignment_dir):
    import csv
    result = {}
    with Path(path).open(newline='') as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if not row:
                continue
            if len(row) != 2 or not row[0] or row[0] in result:
                raise ValueError(f'Invalid/duplicate alignment manifest row: {row}')
            target = Path(alignment_dir) / row[1]
            if not target.is_file():
                raise FileNotFoundError(target)
            result[row[0]] = target
    return result


def infer_alignment_format(alignment_file: Path, format: str = "auto") -> str:
    """Infer a BioPython alignment format from the filename when requested."""
    if format != "auto":
        return format

    with Path(alignment_file).open() as fh:
        first = next((line.strip() for line in fh if line.strip()), '')
    if first.startswith('>'):
        return 'fasta'
    if first.upper().startswith(('CLUSTAL', 'MUSCLE')):
        return 'clustal'
    return 'phylip-relaxed'


def find_gene_alignment(
    alignment_dir: Optional[Path],
    gene: str,
    ensembl_genes: Optional[Set[str]] = None,
) -> Path:
    """Find the alignment file for a specific gene.

    :param alignment_dir: Directory containing alignment files. Supports recursive search for
        PHYLIP-like and FASTA extensions.
    :type alignment_dir: Optional[Path]
    :param gene: Complete gene identifier.
    :type gene: str
    :param ensembl_genes: Optional set of Ensembl gene IDs to limit allowed genes.
    :type ensembl_genes: Optional[Set[str]]
    :returns: Path to the matching alignment file.
    :rtype: Path
    :raises FileNotFoundError: If no suitable alignment file is found or gene not allowed by ensembl_genes.
    """
    if ensembl_genes is not None and gene not in ensembl_genes:
        raise FileNotFoundError(
            f"Gene '{gene}' not present in provided Ensembl list; refusing to match alignments"
        )

    if alignment_dir and alignment_dir.exists():
        # Complete identifiers must match exactly.
        path = _scan_alignment_dir(alignment_dir).get(gene)
        if path is not None:
            return path

    raise FileNotFoundError(f"No alignment file found for gene {gene}")


def read_alignment(
    alignment_file: Path,
    format: str = "auto",
) -> MultipleSeqAlignment:
    """Read a multiple sequence alignment file using BioPython.

    :param alignment_file: Path to alignment file (Path object or string)
    :type alignment_file: Path
    :param format: Alignment format (default: auto)
    :type format: str
    :returns: BioPython MultipleSeqAlignment object
    :rtype: MultipleSeqAlignment
    :raises FileNotFoundError: If alignment file doesn't exist
    :raises ValueError: If alignment cannot be parsed or format invalid
    """
    # Accept a str as well as a Path
    if isinstance(alignment_file, str):
        alignment_file = Path(alignment_file)

    if not alignment_file.exists():
        raise FileNotFoundError(f"Alignment file not found: {alignment_file}")

    try:
        alignment = AlignIO.read(alignment_file, infer_alignment_format(alignment_file, format))
        logger.debug(
            f"Read alignment from {alignment_file}: {len(alignment)} sequences, {alignment.get_alignment_length()} positions"
        )
        return alignment
    except Exception as e:
        raise ValueError(f"Failed to read alignment {alignment_file}: {e}")
