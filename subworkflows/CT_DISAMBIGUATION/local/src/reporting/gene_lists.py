"""
Gene List Export Utilities
==========================

Exports per-pattern and per-change gene lists for CAAS convergence analysis.
Generates plain-text lists and TSV tables for all positions and significance-filtered
subsets using explicit, deterministic rules.

Output Structure
----------------
The exporter creates two organizational schemes under ``gene_reports/``:

1. **by_pattern**: Gene lists organized by convergence pattern type::

    gene_reports/all/by_pattern/<pattern>/gene_list.txt
    gene_reports/significant/by_pattern/<pattern>/gene_list.txt

2. **by_change**: TSV tables organized by which side changed::

    gene_reports/all/by_change/top.tsv
    gene_reports/all/by_change/bottom.tsv
    gene_reports/all/by_change/both.tsv
    gene_reports/significant/by_change/top.tsv
    gene_reports/significant/by_change/bottom.tsv
    gene_reports/significant/by_change/both.tsv

Each TSV has columns: ``gene``, ``position``, ``pattern``

Significance Rule
-----------------
``is_significant`` AND ``is_stable``

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-08
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple, Union

__all__ = [
    "export_gene_lists",
    "extract_gene_sets",
    "extract_change_records",
]

logger = logging.getLogger(__name__)

# Canonical ordering used to keep outputs deterministic; any additional patterns
# discovered at runtime are appended alphabetically afterward.
DEFAULT_PATTERN_ORDER: Sequence[str] = (
    "convergent",
    "divergent",
    "parallel_convergence",
    "parallel_divergence",
    "parallel_mixed",
    "parallel_codivergent",
    "codivergent",
    "codivergent_top",
    "codivergent_bottom",
    "ambiguous",
    "no_change",
)


def _coerce_bool(value: Any) -> bool:
    """
    Normalize truthy/falsey values that may arrive as strings or numbers.

    Args:
        value: Input value to normalize.

    Returns:
        Boolean interpretation of the value.
    """

    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return False


def _safe_pattern_name(pattern: str) -> str:
    """
    Sanitize pattern strings for filesystem-safe directory names.

    Args:
        pattern: Raw pattern label.

    Returns:
        Normalized pattern label safe for directory creation.
    """

    pattern = pattern.replace("/", "_").replace("\\", "_")
    pattern = pattern.replace(" ", "_")
    return pattern


def _normalize_results(
    results: Union[Sequence[Mapping[str, Any]], Mapping[str, Sequence[Mapping[str, Any]]]]
) -> List[MutableMapping[str, Any]]:
    """
    Flatten results into a list while ensuring each entry has a ``gene`` field.

    Accepts either a list of result dicts or a mapping of ``gene -> [results]``.

    Args:
        results: CAAS analysis results.

    Returns:
        List of mutable mappings with guaranteed ``gene`` and ``pattern_type`` keys.

    Raises:
        ValueError: If required fields are missing.
    """

    normalized: List[MutableMapping[str, Any]] = []

    if isinstance(results, Mapping):
        for gene, entries in results.items():
            for entry in entries:
                if not isinstance(entry, MutableMapping):
                    raise ValueError("Each result entry must be a mapping.")
                item = dict(entry)
                item.setdefault("gene", gene)
                normalized.append(item)
    else:
        for entry in results:
            if not isinstance(entry, MutableMapping):
                raise ValueError("Each result entry must be a mapping.")
            item = dict(entry)
            normalized.append(item)

    for item in normalized:
        if "gene" not in item or not item["gene"]:
            raise ValueError("Result entry is missing required 'gene' field.")
        if "pattern_type" not in item or not item["pattern_type"]:
            raise ValueError("Result entry is missing required 'pattern_type' field.")

    return normalized


def _build_pattern_order(entries: Iterable[Mapping[str, Any]], patterns: Optional[Sequence[str]]) -> List[str]:
    """
    Build deterministic pattern order combining defaults, user-specified, and observed.

    Args:
        entries: Iterable of result mappings.
        patterns: Optional explicit pattern order.

    Returns:
        Ordered list of pattern labels to export (excluding the ``all`` catch-all).
    """

    observed = {str(item.get("pattern_type", "unknown")) for item in entries}
    base_order = list(patterns) if patterns else list(DEFAULT_PATTERN_ORDER)

    # Append any missing observed patterns alphabetically to keep deterministic.
    missing = sorted([p for p in observed if p not in base_order])
    return base_order + missing


def _significance_passes(entry: Mapping[str, Any]) -> bool:
    """
    Evaluate the significance rule for exporting genes.

    Rule: is_significant AND is_stable
    Missing fields default to False to enforce conservative filtering.

    Args:
        entry: Result entry.

    Returns:
        True if entry qualifies for the "significant" gene lists.
    """

    is_significant = _coerce_bool(entry.get("is_significant"))
    is_stable = _coerce_bool(entry.get("is_stable"))
    return bool(is_significant and is_stable)


def extract_gene_sets(
    results: Union[Sequence[Mapping[str, Any]], Mapping[str, Sequence[Mapping[str, Any]]]],
    patterns: Optional[Sequence[str]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Extract unique gene sets per pattern for all entries and significance-filtered ones.

    Args:
        results: CAAS results as list or mapping of gene -> entries.
        patterns: Optional explicit pattern ordering. Extra observed patterns are appended.

    Returns:
        Tuple of two dictionaries:
        - all_gene_sets: pattern -> set of genes (including ``all``)
        - significant_gene_sets: pattern -> set of genes (including ``all``)
    """

    entries = _normalize_results(results)
    pattern_order = _build_pattern_order(entries, patterns)

    all_gene_sets: Dict[str, Set[str]] = {p: set() for p in pattern_order}
    sig_gene_sets: Dict[str, Set[str]] = {p: set() for p in pattern_order}

    for item in entries:
        gene = str(item["gene"])
        pattern = str(item.get("pattern_type", "unknown"))

        # Populate all entries
        if pattern not in all_gene_sets:
            all_gene_sets[pattern] = set()
            sig_gene_sets.setdefault(pattern, set())
        all_gene_sets[pattern].add(gene)

        # Populate significant entries
        if _significance_passes(item):
            if pattern not in sig_gene_sets:
                sig_gene_sets[pattern] = set()
            sig_gene_sets[pattern].add(gene)

    # Add catch-all buckets
    all_gene_sets["all"] = {g for genes in all_gene_sets.values() for g in genes}
    sig_gene_sets["all"] = {g for genes in sig_gene_sets.values() for g in genes}

    return all_gene_sets, sig_gene_sets


def extract_change_records(
    results: Union[Sequence[Mapping[str, Any]], Mapping[str, Sequence[Mapping[str, Any]]]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Extract gene/position/pattern records grouped by change side (top, bottom, both).

    Args:
        results: CAAS results as list or mapping of gene -> entries.

    Returns:
        Tuple of two lists (all_records, significant_records), each containing dicts with
        keys: gene, position, pattern, change_side.
    """

    entries = _normalize_results(results)

    all_records: List[Dict[str, Any]] = []
    sig_records: List[Dict[str, Any]] = []

    for item in entries:
        gene = str(item["gene"])
        position = item.get("position_zero_based", item.get("position", ""))
        pattern = str(item.get("pattern_type", "unknown"))
        change_side = str(item.get("change_side", "none")).lower()

        # Skip entries with no meaningful change
        if change_side not in {"top", "bottom", "both"}:
            continue

        record = {
            "gene": gene,
            "position": position,
            "pattern": pattern,
            "change_side": change_side,
        }

        all_records.append(record)

        if _significance_passes(item):
            sig_records.append(record)

    return all_records, sig_records


def _write_change_tsv(target_dir: Path, records: List[Dict[str, Any]], change_side: str) -> Path:
    """
    Write a TSV file for a specific change side with gene/position/pattern columns.

    Args:
        target_dir: Directory where TSV will be created.
        records: List of record dicts filtered for the given change_side.
        change_side: The change side label (top, bottom, both).

    Returns:
        Path to the written TSV file.
    """

    target_dir.mkdir(parents=True, exist_ok=True)
    output_path = target_dir / f"{change_side}.tsv"

    # Filter records for this change side
    filtered = [r for r in records if r["change_side"] == change_side]

    # Sort deterministically by gene, then position
    sorted_records = sorted(filtered, key=lambda r: (r["gene"], r["position"]))

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "position", "pattern"], delimiter="\t")
        writer.writeheader()
        for record in sorted_records:
            writer.writerow({
                "gene": record["gene"],
                "position": record["position"],
                "pattern": record["pattern"],
            })

    return output_path


def _write_gene_list(target_dir: Path, genes: Set[str]) -> Path:
    """
    Write a gene list file with deterministic ordering.

    Args:
        target_dir: Directory where ``gene_list.txt`` will be created.
        genes: Set of gene names.

    Returns:
        Path to the written file.
    """

    target_dir.mkdir(parents=True, exist_ok=True)
    output_path = target_dir / "gene_list.txt"
    sorted_genes = sorted(genes)
    output_path.write_text("\n".join(sorted_genes) + ("\n" if sorted_genes else ""), encoding="utf-8")
    return output_path


def export_gene_lists(
    results: Union[Sequence[Mapping[str, Any]], Mapping[str, Sequence[Mapping[str, Any]]]],
    output_root: Path,
    patterns: Optional[Sequence[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Export per-pattern gene lists and per-change TSV tables for all and significant subsets.

    Creates two organizational schemes:
    1. by_pattern/<pattern>/gene_list.txt
    2. by_change/<side>.tsv (with gene/position/pattern columns)

    Args:
        results: CAAS results as list or mapping of gene -> entries.
        output_root: Base directory where ``gene_reports`` will be created.
        patterns: Optional explicit pattern ordering. Extra observed patterns are appended.

    Returns:
        Dictionary with keys ``all`` and ``significant``, each containing:
        - ``by_pattern``: pattern -> Path of gene_list.txt
        - ``by_change``: change_side -> Path of TSV file
    """

    entries = _normalize_results(results)
    pattern_order = _build_pattern_order(entries, patterns)

    all_sets, sig_sets = extract_gene_sets(entries, pattern_order)
    all_records, sig_records = extract_change_records(entries)

    output_paths: Dict[str, Dict[str, Any]] = {
        "all": {"by_pattern": {}, "by_change": {}},
        "significant": {"by_pattern": {}, "by_change": {}},
    }
    base_dir = Path(output_root) / "gene_reports"

    # Export by_pattern gene lists
    for subset_label, gene_sets in (("all", all_sets), ("significant", sig_sets)):
        subset_dir = base_dir / subset_label / "by_pattern"
        for pattern, genes in gene_sets.items():
            safe_pattern = _safe_pattern_name(pattern)
            target_dir = subset_dir / safe_pattern
            path = _write_gene_list(target_dir, genes)
            output_paths[subset_label]["by_pattern"][pattern] = path
        logger.info(
            "✓ Exported by_pattern gene lists for '%s' subset to %s", subset_label, subset_dir
        )

    # Export by_change TSV tables
    for subset_label, records in (("all", all_records), ("significant", sig_records)):
        subset_dir = base_dir / subset_label / "by_change"
        for change_side in ("top", "bottom", "both"):
            path = _write_change_tsv(subset_dir, records, change_side)
            output_paths[subset_label]["by_change"][change_side] = path
        logger.info(
            "✓ Exported by_change TSV tables for '%s' subset to %s", subset_label, subset_dir
        )

    return output_paths
