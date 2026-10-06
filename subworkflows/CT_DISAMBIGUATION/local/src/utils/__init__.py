# __init__.py — Package exports of the I/O, logging and concurrency utilities.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Public utility exports: alignment lookup and reading (io_utils), `configure_logging` (logger) and the worker
concurrency helpers (concurrency).

Imported by: the submodules are imported directly (`src.utils.<module>`) by contract_main.py, observed_b0_main.py,
explain_positions.py, disambiguation_perms_main.py, src/core/driver.py and src/utils/gene_wrapper.py
"""

from .io_utils import (
    find_gene_alignment,
    read_alignment,
)
from .logger import configure_logging
from .concurrency import (
    plan_concurrency,
    init_worker,
    codeml_slot,
)

__all__ = [
    "find_gene_alignment",
    "read_alignment",
    "configure_logging",
    "plan_concurrency",
    "init_worker",
    "codeml_slot",
]
