#!/usr/bin/env python3
# ========================= src/utils/logging_conf.py =========================

from __future__ import annotations
import logging

def configure_logging(level: str = "INFO") -> None:
    level = (level or "INFO").upper()
    lv = getattr(logging, level, logging.INFO)
    logging.basicConfig(
        level=lv,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
