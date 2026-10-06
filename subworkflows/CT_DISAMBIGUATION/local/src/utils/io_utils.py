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
    return suffix in {".phy", ".phylip", ".aln", ".fa", ".fasta", ""} and path.is_file()


@functools.lru_cache(maxsize=32)
def _scan_alignment_dir(alignment_dir: Path) -> Dict[str, Path]:
    """Map every alignment file under `alignment_dir` to its gene prefix, scanning once per directory.

    The scan is memoized because find_gene_alignment runs once per gene (and per chunk of a
    gene), and a recursive glob plus a stat per entry is expensive on a network file system
    when the directory holds thousands of alignments. The key is the file name up to the
    first dot; the first path in sorted order wins when several files share a prefix.
    """
    candidates = sorted(
        path for path in alignment_dir.glob("**/*") if _is_supported_alignment_path(path)
    )
    by_prefix: Dict[str, Path] = {}
    for path in candidates:
        prefix = path.name.split(".", 1)[0]
        # First match in sorted order wins.
        by_prefix.setdefault(prefix, path)
    return by_prefix


def infer_alignment_format(alignment_file: Path, format: str = "auto") -> str:
    """Infer a BioPython alignment format from the filename when requested."""
    if format != "auto":
        return format

    suffix = alignment_file.suffix.lower()
    if suffix in {".fa", ".fasta"}:
        return "fasta"

    return "phylip-relaxed"


def find_gene_alignment(
    alignment_dir: Optional[Path],
    gene: str,
    ensembl_genes: Optional[Set[str]] = None,
) -> Path:
    """Find the alignment file for a specific gene.

    :param alignment_dir: Directory containing alignment files. Supports recursive search for
        PHYLIP-like and FASTA extensions.
    :type alignment_dir: Optional[Path]
    :param gene: Gene name (prefix used to match file name before the first dot).
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
        # Require exact prefix match before the first dot to avoid partials (e.g., MHS1 vs MHS12)
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
