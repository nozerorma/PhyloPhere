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
#
# Desktop launchers (GNOME/KDE .desktop Exec=, double-click from a file
# manager) start this script with a bare session PATH — typically just
# /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin — because they
# don't source ~/.bashrc or ~/.profile. Most micromamba/conda installers only
# put the binary on PATH via the shell rc hook (and this user's ~/.bashrc, like
# most, early-returns for non-interactive shells), so `command -v micromamba`
# silently fails there even though the exact same command works from a
# terminal. Explicitly add the well-known install locations before searching.
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

# Likewise MAMBA_ROOT_PREFIX is normally exported by the same rc-file hook;
# without it micromamba guesses a default (harmless here since it happens to
# match, but it warns on every launch and isn't guaranteed elsewhere), so pin
# it explicitly whenever the standard install layout is present.
if [ -z "${MAMBA_ROOT_PREFIX:-}" ] && [ -d "$HOME/micromamba/envs" ]; then
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
fi

# Qt has no platform-theme plugin loaded by default, so without this the GUI
# renders with Qt's plain built-in look instead of matching Plasma/GNOME (dark
# mode, accent color, fonts) — even though nothing in the app code forces a
# style. conf/phylophere.yml's qt6-main package already ships the
# xdg-desktop-portal theme plugin (reads native theme over the portal, no
# extra system packages needed on modern Plasma/GNOME); only opt in if the
# user hasn't already set their own (e.g. qt6ct).
: "${QT_QPA_PLATFORMTHEME:=xdgdesktopportal}"
export QT_QPA_PLATFORMTHEME

set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

ENV_NAME="phylophere"

# Terminal=false in the installed .desktop entry means stderr below is
# normally invisible ("the shortcut does nothing"). Report failures through
# whatever notifier is available so they're not silently swallowed, in
# addition to printing them (still useful when run from a terminal).
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
