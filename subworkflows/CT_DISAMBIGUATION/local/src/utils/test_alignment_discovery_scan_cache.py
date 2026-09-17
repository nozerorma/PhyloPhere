#!/usr/bin/env python3
"""Regression net for the NFS-scan memoization fix (docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md,
2026-09-17 rework): `find_gene_alignment` (io_utils.py) and `_perms_worker_replay`'s directory-mode
discovery lookup (gene_wrapper.py) used to re-run a full directory scan (`glob("**/*")` /
`iterdir()`) on EVERY call -- an O(N_files) NFS readdir+stat sweep repeated once per gene
pre-Stage-2, once per CHUNK of a gene after it. Both are now memoized per directory via
`functools.lru_cache`. This test proves: (a) lookup results are unchanged (first-match-in-scan-order
semantics preserved), (b) the expensive scan genuinely runs once per directory, not once per call.

Run:  python -m pytest test_alignment_discovery_scan_cache.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # .../local
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import pytest  # noqa: E402

import src.utils.io_utils as io_utils  # noqa: E402
import src.utils.gene_wrapper as gw  # noqa: E402


def _mk_alignment_dir(tmp_path):
    d = tmp_path / "alignments"
    d.mkdir()
    (d / "ABCB1.Homo_sapiens.fa").write_text(">x\nAA\n")
    (d / "ABCB10.Homo_sapiens.fa").write_text(">x\nAA\n")
    return d


def test_find_gene_alignment_correct_and_missing(tmp_path):
    d = _mk_alignment_dir(tmp_path)
    io_utils._scan_alignment_dir.cache_clear()
    got = io_utils.find_gene_alignment(d, "ABCB1")
    assert got.name == "ABCB1.Homo_sapiens.fa"
    with pytest.raises(FileNotFoundError):
        io_utils.find_gene_alignment(d, "NOT_A_GENE")


def test_find_gene_alignment_prefix_exactness_preserved(tmp_path):
    # ABCB1 must not match ABCB10's file -- exact prefix-before-first-dot only.
    d = _mk_alignment_dir(tmp_path)
    io_utils._scan_alignment_dir.cache_clear()
    got1 = io_utils.find_gene_alignment(d, "ABCB1")
    got10 = io_utils.find_gene_alignment(d, "ABCB10")
    assert got1.name != got10.name


def test_find_gene_alignment_scans_directory_once(tmp_path, monkeypatch):
    d = _mk_alignment_dir(tmp_path)
    io_utils._scan_alignment_dir.cache_clear()

    calls = {"n": 0}
    real_glob = Path.glob

    def counting_glob(self, pattern):
        calls["n"] += 1
        return real_glob(self, pattern)

    monkeypatch.setattr(Path, "glob", counting_glob)

    io_utils.find_gene_alignment(d, "ABCB1")
    io_utils.find_gene_alignment(d, "ABCB10")
    io_utils.find_gene_alignment(d, "ABCB1")  # repeat, still cached

    assert calls["n"] == 1, "directory should be scanned exactly once across 3 lookups"


def test_scan_perm_discovery_dir_correct_and_cached(tmp_path, monkeypatch):
    d = tmp_path / "perm_disc"
    d.mkdir()
    (d / "ABCB1.perm_replay.discovery.output").write_text("cycle\tposition\nb_1~H1\t10\n")
    (d / "ABCB10.perm_replay.discovery.output").write_text("cycle\tposition\nb_1~H1\t20\n")
    gw._scan_perm_discovery_dir.cache_clear()

    got = gw._scan_perm_discovery_dir(d)
    assert set(got.keys()) == {"ABCB1", "ABCB10"}
    assert got["ABCB1"].name.startswith("ABCB1.")

    calls = {"n": 0}
    real_iterdir = Path.iterdir

    def counting_iterdir(self):
        calls["n"] += 1
        return real_iterdir(self)

    monkeypatch.setattr(Path, "iterdir", counting_iterdir)
    gw._scan_perm_discovery_dir(d)  # same path -> cache hit, no new iterdir() call
    assert calls["n"] == 0
