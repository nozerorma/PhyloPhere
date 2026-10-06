#!/usr/bin/env python3
# secrets_io.py — Write and remove the Seqera/Tower access token file.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Secrets: write and remove `token.tk`, the Seqera/Tower access token file.

The tower{} block of conf/common.config reads a gitignored `token.tk` at the repo root
(or, without the file, the TOWER_ACCESS_TOKEN environment variable). The GUI writes the
token to that file and never into a generated script. The token is not part of
ProjectConfig (see gui/models/runtime.py), so it cannot reach the JSON project file.

`repo_dir` is the checkout that will run the pipeline. With General > Remote host set it
is a path on the cluster, not on the machine running the GUI, and `remote_host` must be
passed for the file to be written there over SSH (gui/remote.py).

Imported by: gui/widgets/tabs/runtime_tab.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import stat
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui import remote

TOKEN_FILENAME = "token.tk"


def write_tower_token(repo_dir: Path | str, token: str, remote_host: str = "") -> str:
    """Write `token` to <repo_dir>/token.tk with owner-only permissions (mode 600).

    The file goes to `remote_host` over SSH when given, else to the local filesystem.
    Returns the path written.
    """
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
    """Remove <repo_dir>/token.tk if present (used when the token field is empty)."""
    if remote_host:
        remote_path = f"{str(repo_dir).rstrip('/')}/{TOKEN_FILENAME}"
        remote.remove_remote_file(remote_host, remote_path)
        return
    token_path = Path(repo_dir) / TOKEN_FILENAME
    token_path.unlink(missing_ok=True)
