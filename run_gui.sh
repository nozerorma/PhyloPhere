#!/bin/bash
# run_gui.sh — Launches the PhyloPhere Runner GUI.
#
# Double-click this file in a file manager (if it's configured to run
# executable .sh scripts) or run it from a terminal: ./run_gui.sh
#
# Requires the `phylophere` conda/mamba environment (environment/install_env.sh)
# — it already pins PySide6 and Jinja2, so no separate GUI install step is needed.
#
# Uses `<tool> run -n phylophere` rather than `conda activate` — activation
# depends on shell hooks that are only loaded for interactive shells, which
# breaks when this script is launched non-interactively (double-clicked from a
# file manager, or via the .desktop entry from install_gui_launcher.sh). `run -n`
# works regardless of whether the current shell has ever been hooked up.

set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

ENV_NAME="phylophere"

if command -v micromamba >/dev/null 2>&1; then
    RUN=(micromamba run -n "$ENV_NAME")
elif command -v mamba >/dev/null 2>&1; then
    RUN=(mamba run -n "$ENV_NAME")
elif command -v conda >/dev/null 2>&1; then
    RUN=(conda run -n "$ENV_NAME")
else
    echo "error: none of micromamba, mamba, or conda were found on PATH." >&2
    echo "Set up the environment first with: ./environment/install_env.sh" >&2
    exit 1
fi

if ! "${RUN[@]}" python --version >/dev/null 2>&1; then
    echo "error: the '$ENV_NAME' environment was not found." >&2
    echo "Set it up first with: ./environment/install_env.sh" >&2
    exit 1
fi

exec "${RUN[@]}" python -m gui.main "$@"
