# logger.py — Single entry point for the logging configuration of the disambiguation scripts.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Logging configuration: one uniform format on stdout and, optionally, a log file.

Handlers are synchronous, so the lines of the different worker processes stay
traceable through the processName field of the format.

Usage example::

    import logging
    from src.utils.logger import configure_logging

    configure_logging(verbose=True)
    log = logging.getLogger("my.module")
    log.debug("Hello world")

Imported by: contract_main.py, disambiguation_perms_main.py, explain_positions.py,
             observed_b0_main.py, src/utils/concurrency.py (init_worker)
"""

from __future__ import annotations

import logging
import logging.config
from pathlib import Path
from typing import Optional

# One line per record: time | level | process | logger name | message.
LOG_FORMAT = "%(asctime)s | %(levelname)s | %(processName)s | %(name)s | %(message)s"


def configure_logging(
    *,
    verbose: bool = False,
    log_file: Optional[Path] = None,
    quiet_matplotlib: bool = True,
) -> None:
    """Configure root logging with a uniform, synchronous formatter.

    Replaces the root handlers on every call, so it is safe to call again in a worker process.

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
                "matplotlib.font_manager": {
                    "level": "WARNING" if quiet_matplotlib else level
                },
            },
        }
    )
