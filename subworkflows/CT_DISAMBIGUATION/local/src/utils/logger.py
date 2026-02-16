"""
Centralized logging configuration for CAAS pipelines
===============================================

Provides a single entry point to configure logging for both aggregation and
single-gene workflows. Uses synchronous stream/file handlers with a uniform
format to keep logs traceable across processes.

Usage Example
-------------
::

    from src.utils.logger import configure_logging, get_logger

    configure_logging(verbose=True)
    log = get_logger("my.module")
    log.debug("Hello world")

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-09
"""

from __future__ import annotations

import logging
import logging.config
from pathlib import Path
from typing import Optional

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(processName)s | %(name)s | %(message)s"


def configure_logging(
    *,
    verbose: bool = False,
    log_file: Optional[Path] = None,
    quiet_matplotlib: bool = True,
) -> None:
    """Configure root logging with a uniform, synchronous formatter.

    :param verbose: When True, set root log level to DEBUG; otherwise INFO.
    :type verbose: bool
    :param log_file: Optional path to also tee logs to a file (overwrite).
    :type log_file: Optional[Path]
    :param quiet_matplotlib: If True, silence matplotlib loggers (WARNING level).
    :type quiet_matplotlib: bool
    :returns: None
    :rtype: None
    :example: ::

        configure_logging(verbose=True, log_file=Path('logs/app.log'))
    """

    level = logging.DEBUG if verbose else logging.INFO

    handlers = {
        "console": {
            "class": "logging.StreamHandler",
            "level": level,
            "formatter": "standard",
            "stream": "ext://sys.stdout",
        },
    }
    root_handlers = ["console"]

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers["file"] = {
            "class": "logging.FileHandler",
            "level": level,
            "formatter": "standard",
            "filename": str(log_path),
            "mode": "w",
            "encoding": "utf-8",
        }
        root_handlers.append("file")

    logging.config.dictConfig(
        {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {"standard": {"format": LOG_FORMAT}},
            "handlers": handlers,
            "root": {"level": level, "handlers": root_handlers},
            "loggers": {
                "matplotlib": {"level": "WARNING" if quiet_matplotlib else level},
                "matplotlib.font_manager": {"level": "WARNING" if quiet_matplotlib else level},
            },
        }
    )


def get_logger(name: str) -> logging.Logger:
    """Get a namespaced logger.

    :param name: Logger namespace to fetch (e.g., 'src.utils').
    :type name: str
    :returns: Logger instance with the requested name.
    :rtype: logging.Logger
    :example: ::

        log = get_logger('src.utils')
        log.info('message')
    """
    return logging.getLogger(name)
