"""Concurrency Helpers for ASR Pipeline
=======================================

Unified helpers for planning CPU concurrency, initializing worker processes, and managing
thread limits to prevent oversubscription in ASR computations (e.g., codeml runs).

Usage Example
-------------
::

    # Plan concurrency for 4 threads per gene, auto workers
    workers, threads = plan_concurrency(None, 4, logger=my_logger)
    print(f"Using {workers} workers with {threads} threads each")

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-09
"""

import multiprocessing as mp
import os
from contextlib import contextmanager
from typing import Optional, Tuple

_CODEML_SEM = None  # set in worker initializer


def _cpu_available() -> int:
    """Return the number of CPUs available to the current process."""
    try:
        return len(os.sched_getaffinity(0))  # respects taskset/cgroup limits
    except Exception:
        return mp.cpu_count()


def plan_concurrency(requested_workers: Optional[int], threads_per_gene: int, logger=None) -> Tuple[int, int]:
    """Derive a safe (workers, threads) plan to reduce oversubscription.

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
    effective_workers = max_workers if requested_workers is None else max(1, min(requested_workers, max_workers))

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
    """Initializer for worker processes to set thread env and optional gate.

    Side effects: sets the OMP_NUM_THREADS environment variable and may set
    the module-global semaphore used by :func:`codeml_slot`.

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


@contextmanager
def codeml_slot():
    """Context manager that limits concurrent codeml runs when a shared semaphore is provided.

    If no semaphore is set (e.g., running locally without gating), this is a no-op.

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

