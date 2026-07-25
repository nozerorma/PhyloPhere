#!/usr/bin/env python3
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/PhyloPhere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: gui/main.py
#

"""
PhyloPhere Runner GUI: generates the SBATCH array-job wrapper and single-phenotype
runner scripts from a multi-tab desktop form, instead of hand-editing them.

Usage:
    python -m gui.main [--project PROJECT.json]

Run `python -m gui.main --help` for the full option list.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import sys
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import project_io
from gui.app import run as run_app


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="PhyloPhere Runner GUI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--project", type=Path, help="Project JSON file to open on startup")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.project is not None and not args.project.is_file():
        print(f"error: project file not found: {args.project}", file=sys.stderr)
        sys.exit(1)
    sys.exit(run_app(initial_project_path=args.project))


if __name__ == "__main__":
    main()
