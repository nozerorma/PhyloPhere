#!/usr/bin/env python3
"""Regression net for the Tier 1 fix in
docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md: `_perms_worker` used to write each
cycle's (fg, bg) labeling to a trait file and immediately re-read+re-parse it via
`parse_trait_pairs` before every call to `analyze_gene_disambiguation`. That
round-trip is now skipped: `_perms_worker` builds the equivalent in-memory
`trait_pairs = {1: list(zip(fg, bg))}` dict directly and passes it via the new
`trait_pairs=` kwarg on `analyze_gene_disambiguation`.

This test proves the in-memory shape is exactly what `parse_trait_pairs` would
have derived from the file it replaces, across the shapes that actually occur in
production: equal-length fg/bg, and the FOP-mirror case where the cycle tag (and
therefore the trait filename) embeds an "H<n>" hypothesis token that used to
change parse_trait_pairs' contrast *key* (never the pairs themselves).

Run:  python -m pytest test_perms_worker_trait_pairs.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

_LOCAL = Path(__file__).resolve().parents[1]  # .../src
if str(_LOCAL) not in sys.path:
    sys.path.insert(0, str(_LOCAL))

from data.loaders import parse_trait_pairs  # noqa: E402


def _write_legacy_trait_file(fg, bg, out_path: Path) -> None:
    """Reproduces the exact on-disk format `_write_cycle_trait_file` used to
    write (gene_wrapper.py, pre-Tier-1): species, trait (1=fg/high, 0=bg/low),
    pair index = position in the resample lists (fg[k] <-> bg[k] is pair k+1)."""
    with open(out_path, "w") as f:
        for k, sp in enumerate(fg, start=1):
            f.write(f"{sp}\t1\t{k}\n")
        for k, sp in enumerate(bg, start=1):
            f.write(f"{sp}\t0\t{k}\n")


def _in_memory_trait_pairs(fg, bg):
    """Exactly what _perms_worker now builds (gene_wrapper.py)."""
    return {1: list(zip(fg, bg))}


def test_equal_length_fg_bg_matches_file_based_parse():
    fg = ["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla"]
    bg = ["Mus_musculus", "Rattus_norvegicus", "Cavia_porcellus"]
    with tempfile.TemporaryDirectory() as td:
        trait_path = Path(td) / "trait_cyc1.tab"
        _write_legacy_trait_file(fg, bg, trait_path)
        expected = parse_trait_pairs(trait_path)
    got = _in_memory_trait_pairs(fg, bg)
    assert expected == got


def test_mismatched_length_drops_the_same_unpaired_tail():
    # fg longer than bg: parse_trait_pairs only keeps pair_ids present on BOTH
    # sides (disambiguate_single.py's `if high_species and low_species`), which
    # is exactly what zip()'s shorter-iterable truncation reproduces.
    fg = ["A", "B", "C", "D", "E"]
    bg = ["X", "Y", "Z"]
    with tempfile.TemporaryDirectory() as td:
        trait_path = Path(td) / "trait_cyc2.tab"
        _write_legacy_trait_file(fg, bg, trait_path)
        expected = parse_trait_pairs(trait_path)
    got = _in_memory_trait_pairs(fg, bg)
    assert expected == got
    assert len(got[1]) == 3


def test_fop_mirror_filename_changes_only_the_contrast_key_not_the_pairs():
    """FOP-mirror cycle tags look like "<base>~H<m>"; the trait filename used to
    be f"trait_{cycle_tag}.tab", so parse_trait_pairs' `re.search(r"H(\\d+)", ...)`
    picks up that "H<m>" and uses contrast key `m` instead of the plain-file
    default of 1. The in-memory dict always uses key 1 regardless. This is safe
    because analyze_gene_disambiguation's `_resolve_contrast` only reads a CAAS
    metadata row's OWN `trait` field for an "H<n>" tag (never the dict's key) and
    falls back to "the lone contrast" whenever there's exactly one — which is
    always true for a single replayed cycle. This test locks in the one thing
    that must stay true for that safety argument to hold: the pairs list itself
    is identical no matter which key wraps it.
    """
    fg = ["Homo_sapiens", "Pan_troglodytes"]
    bg = ["Mus_musculus", "Rattus_norvegicus"]
    with tempfile.TemporaryDirectory() as td:
        trait_path = Path(td) / "trait_cyc3~H2.tab"
        _write_legacy_trait_file(fg, bg, trait_path)
        expected = parse_trait_pairs(trait_path)
    assert list(expected.keys()) == [2]  # confirms the H-token IS picked up
    got = _in_memory_trait_pairs(fg, bg)
    assert list(expected.values()) == list(got.values())
