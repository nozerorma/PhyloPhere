#!/usr/bin/env python3
# remote.py — SSH-based path existence checks + directory listing for remote hosts.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Dataset paths often live on a remote HPC cluster (reached over SSH) while the GUI
itself runs on a laptop. This module lets "Validate Paths" and the PathField
"Browse..." button work against a remote host instead of the local filesystem.

Deliberately separate from gui/generation/ (which must stay network-call-free —
its tests are hermetic and fast) and has no PySide6 import, so it's independently
unit-testable by mocking subprocess.run.

Assumes passwordless SSH key-based auth is already configured for the target host
(the norm for HPC cluster access) — runs ssh in BatchMode so a misconfigured host
fails fast with a clear error instead of hanging on a password prompt.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import shlex
import subprocess

DEFAULT_TIMEOUT = 15


class RemoteCheckError(Exception):
    """SSH itself failed (unreachable host, auth failure, timeout) — distinct from
    a path simply not existing, which is a normal (non-exceptional) result."""


def _run_ssh(host: str, remote_command: str, *, stdin: str | None, timeout: int) -> str:
    try:
        result = subprocess.run(
            [
                "ssh",
                "-o", "BatchMode=yes",
                "-o", f"ConnectTimeout={min(timeout, 30)}",
                "-o", "StrictHostKeyChecking=accept-new",
                host,
                remote_command,
            ],
            input=stdin,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise RemoteCheckError("ssh is not installed or not on PATH.") from exc
    except subprocess.TimeoutExpired as exc:
        raise RemoteCheckError(f"SSH to {host!r} timed out after {timeout}s.") from exc

    if result.returncode != 0:
        stderr = result.stderr.strip()
        raise RemoteCheckError(
            f"SSH to {host!r} failed (exit {result.returncode})."
            + (f" {stderr}" if stderr else " Is passwordless key-based auth configured?")
        )
    return result.stdout


def check_remote_paths(
    host: str, entries: list[tuple[str, str, str]], timeout: int = DEFAULT_TIMEOUT
) -> list[str]:
    """Like gui.generation.validate.validate_paths, but checks existence on `host`
    over SSH in a single round-trip instead of the local filesystem."""
    filled = [(label, path, kind) for label, path, kind in entries if path.strip()]
    if not filled:
        return []

    script_lines = []
    for i, (_label, path, kind) in enumerate(filled):
        flag = "-f" if kind == "file" else "-d"
        script_lines.append(f"test {flag} {shlex.quote(path)} || echo MISSING:{i}")
    script = "\n".join(script_lines)

    stdout = _run_ssh(host, "bash -s", stdin=script, timeout=timeout)

    missing_indices = set()
    for line in stdout.splitlines():
        if line.startswith("MISSING:"):
            missing_indices.add(int(line.split(":", 1)[1]))

    problems = []
    for i, (label, path, kind) in enumerate(filled):
        if i in missing_indices:
            noun = "file" if kind == "file" else "directory"
            problems.append(f"{label}: {noun} not found on {host} — {path}")
    return problems


def list_remote_directory(
    host: str, path: str, timeout: int = DEFAULT_TIMEOUT
) -> list[tuple[str, bool]]:
    """Returns [(name, is_dir), ...] for `path` on `host`, sorted directories-first.

    Uses `find -L` (follows symlinks) so symlinked data mounts — common on cluster
    filesystems — are classified as navigable directories rather than opaque files.
    """
    remote_command = (
        f"find -L {shlex.quote(path)} -mindepth 1 -maxdepth 1 "
        f"-printf '%y\\t%f\\n' 2>/dev/null | sort"
    )
    stdout = _run_ssh(host, remote_command, stdin=None, timeout=timeout)

    entries = []
    for line in stdout.splitlines():
        kind, _, name = line.partition("\t")
        if name:
            entries.append((name, kind == "d"))
    entries.sort(key=lambda e: (not e[1], e[0].lower()))  # directories first, then A-Z
    return entries
