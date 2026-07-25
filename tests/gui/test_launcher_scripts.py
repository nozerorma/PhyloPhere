#!/usr/bin/env python3
# test_launcher_scripts.py — Syntax and structure checks for the distribution launcher scripts.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
No AppImage/PyInstaller packaging — distribution is `run_gui.sh` (runs the GUI via
the existing `phylophere` env's PySide6/Jinja2) + `install_gui_launcher.sh`
(optional .desktop entry). These just need to stay valid, executable bash.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _bash_syntax_ok(path: Path) -> bool:
    result = subprocess.run(["bash", "-n", str(path)], capture_output=True)
    return result.returncode == 0


def test_run_gui_script_exists_and_is_executable():
    path = REPO_ROOT / "run_gui.sh"
    assert path.is_file()
    assert path.stat().st_mode & 0o111  # at least one executable bit set
    assert _bash_syntax_ok(path)


def test_run_gui_script_falls_back_across_micromamba_mamba_conda():
    text = (REPO_ROOT / "run_gui.sh").read_text()
    assert "micromamba" in text
    assert "mamba" in text
    assert "conda" in text
    # Must actually invoke `<tool> run -n phylophere ...` for the real command
    # (not just `conda activate`, which depends on interactive shell hooks that
    # aren't loaded when this is double-clicked from a file manager).
    assert 'RUN=(micromamba run -n "$ENV_NAME")' in text
    assert 'RUN=(mamba run -n "$ENV_NAME")' in text
    assert 'RUN=(conda run -n "$ENV_NAME")' in text
    assert '"${RUN[@]}" python -m gui.main' in text


def test_install_gui_launcher_script_exists_and_is_executable():
    path = REPO_ROOT / "install_gui_launcher.sh"
    assert path.is_file()
    assert path.stat().st_mode & 0o111
    assert _bash_syntax_ok(path)


def test_install_gui_launcher_references_run_gui_and_logo():
    text = (REPO_ROOT / "install_gui_launcher.sh").read_text()
    assert "run_gui.sh" in text
    assert "res/logo.png" in text
    assert "Exec=" in text
    assert "Icon=" in text
