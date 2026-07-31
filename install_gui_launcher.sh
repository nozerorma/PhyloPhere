#!/bin/bash
# install_gui_launcher.sh — Adds a "PhyloPhere Runner GUI" entry to your desktop's
# application menu, so you can launch it like any other installed app.
#
# Optional — run_gui.sh works fine on its own from a terminal or file manager.
# This just makes the GUI show up in menus/launchers (GNOME Activities, KDE
# app menu, etc.) without needing to know where the repo lives.
#
# Safe to re-run (e.g. after moving the repo) — it always regenerates the entry.

set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APPS_DIR="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
DESKTOP_FILE="$APPS_DIR/phylophere-gui.desktop"

mkdir -p "$APPS_DIR"

cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Type=Application
Name=PhyloPhere Runner GUI
Comment=Generate PhyloPhere SBATCH/runner scripts
Exec=$REPO_DIR/run_gui.sh
Icon=$REPO_DIR/res/icon.png
Terminal=false
Categories=Science;Biology;
EOF

chmod +x "$DESKTOP_FILE"

if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database "$APPS_DIR" 2>/dev/null || true
fi

echo "Installed: $DESKTOP_FILE"
echo "The PhyloPhere Runner GUI should now appear in your application menu."
