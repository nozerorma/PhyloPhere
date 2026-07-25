#!/usr/bin/env python3
# secrets_io.py — Seqera/Tower access-token handling.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Reuses the pipeline's existing convention: conf/common.config's `tower {}` block
already reads a gitignored `token.tk` at the repo root automatically
(`accessToken = new File("$baseDir/token.tk").exists() ? ... : System.getenv(...)`).
The GUI never needs to write the token into a generated script — it just writes
straight to that file. `token.tk` is confirmed present in .gitignore.

The token itself is intentionally NOT part of ProjectConfig (see gui/models/runtime.py)
so it can never end up in the human-diffable JSON project file.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import stat
from pathlib import Path

TOKEN_FILENAME = "token.tk"


def write_tower_token(repo_dir: Path, token: str) -> Path:
    """Write `token` to <repo_dir>/token.tk with owner-only permissions."""
    token_path = Path(repo_dir) / TOKEN_FILENAME
    token_path.write_text(token.strip() + "\n")
    token_path.chmod(stat.S_IRUSR | stat.S_IWUSR)  # 0o600
    return token_path


def clear_tower_token(repo_dir: Path) -> None:
    """Remove <repo_dir>/token.tk, if present (e.g. user unchecks 'use Tower')."""
    token_path = Path(repo_dir) / TOKEN_FILENAME
    token_path.unlink(missing_ok=True)
