"""Public utility exports for I/O, logging, DB helpers, and wrappers."""

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


def convert_convergence_result_to_dict(*args, **kwargs):
    """Proxy to :func:`src.utils.gene_wrapper.convert_convergence_result_to_dict`."""
    from .gene_wrapper import convert_convergence_result_to_dict as _impl

    return _impl(*args, **kwargs)


def convert_biochem_result_to_dict(*args, **kwargs):
    """Proxy to :func:`src.utils.gene_wrapper.convert_biochem_result_to_dict`."""
    from .gene_wrapper import convert_biochem_result_to_dict as _impl

    return _impl(*args, **kwargs)


def merge_multi_hypothesis_results(*args, **kwargs):
    """Proxy to :func:`src.utils.gene_wrapper.merge_multi_hypothesis_results`."""
    from .gene_wrapper import merge_multi_hypothesis_results as _impl

    return _impl(*args, **kwargs)


def process_single_gene(*args, **kwargs):
    """Proxy to :func:`src.utils.gene_wrapper.process_single_gene`."""
    from .gene_wrapper import process_single_gene as _impl

    return _impl(*args, **kwargs)


def process_all_genes(*args, **kwargs):
    """Proxy to :func:`src.utils.gene_wrapper.process_all_genes`."""
    from .gene_wrapper import process_all_genes as _impl

    return _impl(*args, **kwargs)


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
    "convert_convergence_result_to_dict",
    "convert_biochem_result_to_dict",
    "merge_multi_hypothesis_results",
    "process_single_gene",
    "process_all_genes",
]
