"""
Utils Module
============

Utilities for the CAAS aggregation pipeline including I/O, logging, and
database helpers.

Components:
    - logging_conf: Centralized logging configuration
    - io_utils: File readers/writers for alignments, trees, metadata
    - validation: Input validation and error handling

Usage Example
-------------
::

    from src.utils import find_gene_alignment, read_alignment

    path = find_gene_alignment(Path("alignments/"), "BRCA2")
    alignment = read_alignment(path)

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-09
"""


from .io_utils import (
    find_gene_alignment,
    read_alignment,
)
from .logger import (
    configure_logging,
    get_logger,
)
from .concurrency import (
    plan_concurrency,
    init_worker,
    codeml_slot,
)
from .disambiguation_db import (
    init_db,
    get_connection,
    insert_gene_alignment,
    insert_result,
    fetch_alignment_for_gene,
    iter_group_keys,
    iter_results_for_group,
)

from .gene_wrapper import (
    convert_biochem_result_to_dict,
    merge_multi_hypothesis_results,
    process_single_gene,
    process_all_genes,
)


__all__ = [
    "find_gene_alignment",
    "read_alignment",
    "configure_logging",
    "get_logger",
    "plan_concurrency",
    "init_worker",
    "codeml_slot",
    "init_db",
    "get_connection",
    "insert_gene_alignment",
    "insert_result",
    "fetch_alignment_for_gene",
    "iter_group_keys",
    "iter_results_for_group",
    "convert_biochem_result_to_dict",
    "merge_multi_hypothesis_results",
    "process_single_gene",
    "process_all_genes",
]
