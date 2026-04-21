# loaders.py — Load and validate pipeline-standard data formats for CT_DISAMBIGUATION.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/data/

"""
Loaders: Centralised I/O for all data types consumed by the disambiguation pipeline.

Defines typed loaders for gene trees, CAAS tables, and alignment blocks so that
parsing details stay in one place and callers receive validated, normalised objects.

Imported by: src/convergence/convergence.py, src/asr/reconstruct.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from pathlib import Path
from typing import Optional

# ── Third-party ───────────────────────────────────────────────────────────────
import pandas as pd
from Bio import Phylo

# ── Package-internal ──────────────────────────────────────────────────────────
from ..utils.logger import get_logger

logger = get_logger(__name__)


# ── Public API ────────────────────────────────────────────────────────────────


def load_caas_table(path: Path) -> pd.DataFrame:
    """Read a CAAS discovery TSV and normalise column names to snake_case.

    The discovery output format changed between caastools v1 and v2;
    this function handles both variants transparently.
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]

    required = {"trait", "gene", "position"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name}: missing required columns {missing}")

    return df


def load_gene_tree(path: Path, gene: Optional[str] = None) -> Phylo.BaseTree.Tree:
    """Parse a Newick gene tree.

    Args:
        path: Path to the Newick file.
        gene: Gene name used only in error messages for clearer diagnostics.
    """
    try:
        trees = list(Phylo.parse(str(path), "newick"))
    except Exception as exc:
        label = f"gene {gene}" if gene else str(path)
        raise IOError(f"Failed to parse tree for {label}: {exc}") from exc

    if not trees:
        raise IOError(f"No trees found in {path}")
    return trees[0]
