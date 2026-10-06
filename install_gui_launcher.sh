#!/bin/bash
# install_gui_launcher.sh — Add a "PhyloPhere Runner GUI" entry to the desktop application menu.
# PhyloPhere | ./
#
# Called by:  the user, once per checkout
# Usage:      ./install_gui_launcher.sh
#
# Optional: run_gui.sh works on its own from a terminal or file manager. This script
# writes a .desktop file (in $XDG_DATA_HOME/applications, default ~/.local/share/applications)
# that points at run_gui.sh and the icon of this checkout, so the GUI appears in menus
# and launchers (GNOME Activities, KDE application menu) without knowing where the repo is.
# The file is rewritten on every run, so re-running after moving the repo updates the paths.

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
