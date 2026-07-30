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

`repo_dir` is "wherever the pipeline checkout that will actually execute lives" —
when General > Remote host is set, that's a path on the cluster reached over SSH,
not on the machine running the GUI (same local/remote split as gui/remote.py's
path validation and browsing). `remote_host` must be passed through explicitly so
this stays true regardless of where the GUI happens to be running.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import stat
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import remote

TOKEN_FILENAME = "token.tk"


def write_tower_token(repo_dir: Path | str, token: str, remote_host: str = "") -> str:
    """Write `token` to <repo_dir>/token.tk with owner-only permissions, on
    `remote_host` over SSH if given, else on the local filesystem."""
    token_content = token.strip() + "\n"
    if remote_host:
        remote_path = f"{str(repo_dir).rstrip('/')}/{TOKEN_FILENAME}"
        remote.write_remote_file(remote_host, remote_path, token_content)
        return remote_path
    token_path = Path(repo_dir) / TOKEN_FILENAME
    token_path.write_text(token_content)
    token_path.chmod(stat.S_IRUSR | stat.S_IWUSR)  # 0o600
    return str(token_path)


def clear_tower_token(repo_dir: Path | str, remote_host: str = "") -> None:
    """Remove <repo_dir>/token.tk, if present (e.g. user unchecks 'use Tower')."""
    if remote_host:
        remote_path = f"{str(repo_dir).rstrip('/')}/{TOKEN_FILENAME}"
        remote.remove_remote_file(remote_host, remote_path)
        return
    token_path = Path(repo_dir) / TOKEN_FILENAME
    token_path.unlink(missing_ok=True)
