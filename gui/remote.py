#!/usr/bin/env python3
# remote.py — Path checks, directory listings and file writes on a remote host over SSH.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Remote: SSH helpers for "Validate Paths", the PathField "Browse..." dialog, the
Regenerate HTML Reports dialog, and writing run scripts and the Tower token on the
remote host.

Dataset paths often live on a remote HPC cluster while the GUI runs on a laptop. The
module is kept apart from gui/generation/, which makes no network calls, and imports no
PySide6, so subprocess.run can be mocked to test it alone.

Authentication must be key-based and already configured for the host. ssh runs in
BatchMode, so a misconfigured host fails at once with an error instead of waiting at a
password prompt. Every failure of ssh itself raises RemoteCheckError.

Imported by: gui/secrets_io.py, gui/widgets/main_window.py,
gui/widgets/common/remote_browse_dialog.py, gui/widgets/common/regenerate_dialog.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import shlex
import subprocess

DEFAULT_TIMEOUT = 15


class RemoteCheckError(Exception):
    """ssh itself failed (unreachable host, authentication failure, timeout).

    A path that does not exist is not an error: it is reported in the return value.
    """


def _run_ssh(host: str, remote_command: str, *, stdin: str | None, timeout: int) -> str:
    """Run `remote_command` on `host` and return its stdout; raise RemoteCheckError on failure."""
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
    """Check that each (label, path, kind) entry exists on `host`, in one SSH round trip.

    Counterpart of gui.generation.validate.validate_paths for a remote filesystem.
    Entries with a blank path are skipped. Returns one message per missing path.
    """
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


def write_remote_file(
    host: str, path: str, content: str, mode: str = "600", timeout: int = DEFAULT_TIMEOUT
) -> None:
    """Write `content` to `path` on `host` and set its permissions to `mode`.

    The default 600 suits secrets (gui/secrets_io.py); the generated run scripts
    use mode="755" (MainWindow._save_generated_scripts).
    """
    remote_command = f"cat > {shlex.quote(path)} && chmod {mode} {shlex.quote(path)}"
    _run_ssh(host, remote_command, stdin=content, timeout=timeout)


def remove_remote_file(host: str, path: str, timeout: int = DEFAULT_TIMEOUT) -> None:
    """Remove `path` on `host`; a file that is already absent is not an error."""
    remote_command = f"rm -f {shlex.quote(path)}"
    _run_ssh(host, remote_command, stdin=None, timeout=timeout)


def list_all_remote(host: str, root: str, timeout: int = DEFAULT_TIMEOUT) -> dict[str, bool]:
    """List everything under `root` on `host` recursively, as {relative_posix_path: is_dir}.

    One SSH round trip. The Regenerate HTML Reports dialog
    (gui/widgets/common/regenerate_dialog.py) uses the listing as a stand-in for the
    local glob search of gui.generation.report_registry when the output directory is
    remote, so that report_registry itself makes no network calls.

    Runs `find -L` (follows symlinks, as list_remote_directory does) with `%y` (type)
    and `%P` (path relative to root), so the results need no path handling here.
    """
    remote_command = (
        f"find -L {shlex.quote(root)} -mindepth 1 -printf '%y\\t%P\\n' 2>/dev/null"
    )
    stdout = _run_ssh(host, remote_command, stdin=None, timeout=timeout)

    entries: dict[str, bool] = {}
    for line in stdout.splitlines():
        kind, _, relpath = line.partition("\t")
        if relpath:
            entries[relpath] = kind == "d"
    return entries


def list_remote_directory(
    host: str, path: str, timeout: int = DEFAULT_TIMEOUT
) -> list[tuple[str, bool]]:
    """List the entries of `path` on `host` as [(name, is_dir), ...], directories first.

    Runs `find -L` (follows symlinks), so symlinked data mounts, common on cluster
    filesystems, are listed as directories that can be entered.
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
