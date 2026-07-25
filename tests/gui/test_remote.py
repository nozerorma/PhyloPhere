#!/usr/bin/env python3
# test_remote.py — SSH-based remote path check/listing tests (subprocess mocked, no network).
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
import subprocess
from unittest.mock import patch

# ── Local ─────────────────────────────────────────────────────────────────────
import pytest

from gui.remote import RemoteCheckError, check_remote_paths, list_remote_directory


def _fake_result(stdout: str = "", stderr: str = "", returncode: int = 0):
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr=stderr)


def test_check_remote_paths_reports_missing_only():
    entries = [
        ("Runtime: tree file", "/data/tree.nwk", "file"),
        ("Runtime: work dir", "/data/work", "dir"),
        ("Runtime: unused", "", "file"),  # empty -> must be skipped, no ssh call needed for it
    ]
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.return_value = _fake_result(stdout="MISSING:1\n")
        problems = check_remote_paths("user@host", entries)

    assert len(problems) == 1
    assert "Runtime: work dir" in problems[0]
    assert "/data/work" in problems[0]
    assert "user@host" in problems[0]

    # Confirm the script sent over stdin tests exactly the 2 non-empty entries.
    call_kwargs = mock_run.call_args.kwargs
    assert "test -f /data/tree.nwk" in call_kwargs["input"]
    assert "test -d /data/work" in call_kwargs["input"]


def test_check_remote_paths_all_exist_returns_empty():
    entries = [("Runtime: tree file", "/data/tree.nwk", "file")]
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.return_value = _fake_result(stdout="")
        assert check_remote_paths("host", entries) == []


def test_check_remote_paths_no_entries_skips_ssh_call():
    with patch("gui.remote.subprocess.run") as mock_run:
        assert check_remote_paths("host", [("x", "", "file")]) == []
        mock_run.assert_not_called()


def test_ssh_connection_failure_raises_remote_check_error():
    entries = [("x", "/data", "dir")]
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.return_value = _fake_result(returncode=255, stderr="Permission denied")
        with pytest.raises(RemoteCheckError, match="Permission denied"):
            check_remote_paths("bad-host", entries)


def test_ssh_timeout_raises_remote_check_error():
    entries = [("x", "/data", "dir")]
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.side_effect = subprocess.TimeoutExpired(cmd="ssh", timeout=15)
        with pytest.raises(RemoteCheckError, match="timed out"):
            check_remote_paths("slow-host", entries)


def test_list_remote_directory_parses_and_sorts_dirs_first():
    stdout = "f\tzeta.txt\nd\talpha\nf\tbeta.csv\nd\tomega\n"
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.return_value = _fake_result(stdout=stdout)
        entries = list_remote_directory("host", "/data")

    assert entries == [("alpha", True), ("omega", True), ("beta.csv", False), ("zeta.txt", False)]


def test_list_remote_directory_uses_find_dash_capital_l_to_follow_symlinks():
    with patch("gui.remote.subprocess.run") as mock_run:
        mock_run.return_value = _fake_result(stdout="")
        list_remote_directory("host", "/data")

    remote_command = mock_run.call_args.args[0][-1]
    assert "find -L" in remote_command
