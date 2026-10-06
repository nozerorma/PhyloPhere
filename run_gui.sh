#!/bin/bash
# run_gui.sh — Launch the PhyloPhere Runner GUI.
# PhyloPhere | ./
#
# Called by:  the user (terminal, file-manager double-click, or the .desktop entry
#             written by install_gui_launcher.sh)
# Usage:      ./run_gui.sh [args passed to gui.main]
#
# Requires the `phylophere` environment (environment/install_env.sh), which pins
# PySide6 and Jinja2 (environment/phylophere.yml), so the GUI needs no separate install.
#
# The GUI starts through `<tool> run -n phylophere` instead of `conda activate`:
# activation relies on shell hooks that only interactive shells load, so it fails when the
# script is launched from a file manager or from the .desktop entry.
#
# ── PATH for desktop launchers ────────────────────────────────────────────────
# Desktop launchers start this script with a bare session PATH (typically
# /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin), because they source
# neither ~/.bashrc nor ~/.profile. Installers usually put micromamba on PATH through
# an rc-file hook that non-interactive shells skip, so `command -v micromamba` would
# fail there although it works in a terminal. The well-known install locations are
# therefore appended to PATH when they exist.
for _dir in \
    "$HOME/.local/bin" \
    "$HOME/micromamba/bin" \
    "$HOME/micromamba/condabin" \
    "$HOME/miniforge3/bin" \
    "$HOME/miniforge3/condabin" \
    "$HOME/miniconda3/bin" \
    "$HOME/miniconda3/condabin" \
    "$HOME/anaconda3/bin" \
    "$HOME/anaconda3/condabin" \
    "$HOME/mambaforge/bin" \
    "/opt/conda/bin"; do
    case ":$PATH:" in
        *":$_dir:"*) ;;
        *) [ -d "$_dir" ] && PATH="$PATH:$_dir" ;;
    esac
done
export PATH

# ── Environment variables ─────────────────────────────────────────────────────
# MAMBA_ROOT_PREFIX is normally exported by the same rc-file hook. Without it
# micromamba guesses a default, warns on every launch and may guess wrong, so it is
# set explicitly when the standard install layout is present.
if [ -z "${MAMBA_ROOT_PREFIX:-}" ] && [ -d "$HOME/micromamba/envs" ]; then
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
fi

# Qt loads no platform-theme plugin by default, so the GUI would render with Qt's plain
# built-in look instead of matching Plasma/GNOME (dark mode, accent color, fonts). The
# qt6-main package of the environment ships the xdg-desktop-portal theme plugin, which
# reads the native theme over the portal. A theme already set by the user (e.g. qt6ct)
# is kept.
: "${QT_QPA_PLATFORMTHEME:=xdgdesktopportal}"
export QT_QPA_PLATFORMTHEME

# ── Launch ────────────────────────────────────────────────────────────────────
set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

ENV_NAME="phylophere"

# The installed .desktop entry runs without a terminal (Terminal=false), so stderr is
# normally invisible and a failure would look like a shortcut that does nothing.
# Failures are therefore also reported through whichever desktop notifier exists.
notify_failure() {
    echo "$1" >&2
    if command -v notify-send >/dev/null 2>&1; then
        notify-send -u critical "PhyloPhere Runner GUI" "$1" || true
    elif command -v zenity >/dev/null 2>&1; then
        zenity --error --title="PhyloPhere Runner GUI" --text="$1" || true
    fi
}

if command -v micromamba >/dev/null 2>&1; then
    RUN=(micromamba run -n "$ENV_NAME")
elif command -v mamba >/dev/null 2>&1; then
    RUN=(mamba run -n "$ENV_NAME")
elif command -v conda >/dev/null 2>&1; then
    RUN=(conda run -n "$ENV_NAME")
else
    notify_failure "none of micromamba, mamba, or conda were found. Set up the environment first with: ../environment/install_env.sh"
    exit 1
fi

if ! "${RUN[@]}" python --version >/dev/null 2>&1; then
    notify_failure "the '$ENV_NAME' environment was not found. Set it up first with: ../environment/install_env.sh"
    exit 1
fi

exec "${RUN[@]}" python -m gui.main "$@"
