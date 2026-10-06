# concurrency.py — CPU concurrency planning and worker initialization for the per-gene pools.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Concurrency helpers: size a worker pool against the CPUs the process may use,
initialize each worker, and optionally gate concurrent codeml runs.

Usage example::

    # Plan for 4 threads per gene, with the number of workers chosen automatically
    workers, threads = plan_concurrency(None, 4, logger=my_logger)

Imported by: src/utils/gene_wrapper.py, src/core/driver.py, explain_positions.py, observed_b0_main.py
"""

import multiprocessing as mp
import os
from contextlib import contextmanager
from typing import Optional, Tuple

# Semaphore that gates concurrent codeml runs; set per worker by init_worker (None = no gate).
_CODEML_SEM = None


def _cpu_available() -> int:
    """Return the number of CPUs available to the current process."""
    try:
        return len(os.sched_getaffinity(0))  # respects taskset/cgroup limits
    except Exception:
        return mp.cpu_count()


def plan_concurrency(
    requested_workers: Optional[int], threads_per_gene: int, logger=None
) -> Tuple[int, int]:
    """Derive a (workers, threads) plan that does not oversubscribe the available CPUs.

    The number of workers is capped at available CPUs // threads, so that
    workers x threads never exceeds the CPUs the process may use.

    :param requested_workers: User-supplied worker count or None for auto.
    :type requested_workers: Optional[int]
    :param threads_per_gene: Requested threads per codeml/ASR run.
    :type threads_per_gene: int
    :param logger: Optional logger for diagnostics.
    :type logger: Optional[logging.Logger]
    :returns: Tuple containing (effective_workers, normalized_threads_per_gene).
    :rtype: Tuple[int, int]
    :example: ::

        workers, threads = plan_concurrency(None, 4)
    """
    threads = max(1, threads_per_gene or 1)
    avail = max(1, _cpu_available())
    max_workers = max(1, avail // threads)
    effective_workers = (
        max_workers
        if requested_workers is None
        else max(1, min(requested_workers, max_workers))
    )

    if logger:
        logger.info(
            "Concurrency plan: requested_workers=%s threads_per_gene=%s available_cpu=%s -> workers=%s",
            requested_workers,
            threads,
            avail,
            effective_workers,
        )
    return effective_workers, threads


def init_worker(threads_per_gene: int, codeml_sem=None) -> None:
    """Initialize a worker process: thread environment, logging and the optional codeml gate.

    Side effects: sets the OMP_NUM_THREADS environment variable, configures
    root logging, and stores the module-global semaphore used by :func:`codeml_slot`.
    Logging must be configured here because a worker started with the spawn or
    forkserver method inherits no logging configuration from the parent; without it
    every `logger.info(...)` of a worker is dropped.

    :param threads_per_gene: Number of threads each worker should use.
    :type threads_per_gene: int
    :param codeml_sem: Optional multiprocessing Semaphore used for gating.
    :type codeml_sem: Optional[multiprocessing.Semaphore]
    :returns: None
    :rtype: None
    """
    os.environ["OMP_NUM_THREADS"] = str(max(1, threads_per_gene))
    global _CODEML_SEM
    _CODEML_SEM = codeml_sem

    from src.utils.logger import configure_logging
    configure_logging()


@contextmanager
def codeml_slot():
    """Limit concurrent codeml runs with the shared semaphore of the worker.

    A no-op when init_worker received no semaphore (no gating).

    :returns: Yields control to a block guarded by the semaphore if present; otherwise, a no-op.
    :rtype: contextlib.AbstractContextManager
    :example: ::

        with codeml_slot():
            run_codeml()
    """
    if _CODEML_SEM is None:
        yield
        return
    with _CODEML_SEM:
        yield
