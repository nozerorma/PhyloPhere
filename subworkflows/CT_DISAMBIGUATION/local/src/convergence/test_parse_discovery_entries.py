#!/usr/bin/env python3
"""`_parse_discovery_entries` — the one parser shared by `_perms_worker`'s two
bootstrap-discovery layouts (V3-6 consolidation of two near-identical copies).

Run: python -m pytest test_parse_discovery_entries.py
"""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # .../local
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import src.utils.gene_wrapper as gw  # noqa: E402


CYCLE_TAGS = {"c1", "c2"}


def _concat_stream():
    # single concatenated perm_discovery_file layout: has a `gene` column
    return io.StringIO(
        "gene\tcycle\tposition\tcaas\tamino_encoded\tis_conserved_meta\tconserved_pair\tcaap_group\n"
        "GENEA\tc1\t10\tV/A\tV\tFALSE\tsp1:sp2\tUS\n"
        "GENEB\tc1\t99\tL/M\tL\tTRUE\tx:y\tGS3\n"          # other gene -> filtered out
        "GENEA\tc2\t11\tI/T\tI\tFALSE\t\tGS2\n"
        "GENEA\tc9\t12\tI/T\tI\tFALSE\t\tUS\n"              # cycle not in tags -> dropped
    )


def _shard_stream():
    # per-gene shard layout: no `gene` column, every row is this gene
    return io.StringIO(
        "cycle\tposition\tcaas\tamino_encoded\tis_conserved_meta\tconserved_pair\tcaap_group\n"
        "c1\t10\tV/A\tV\tFALSE\tsp1:sp2\tUS\n"
        "c2\t11\tI/T\tI\tTRUE\t\tGS2\n"
    )


def test_concat_layout_filters_by_gene_and_cycle():
    out = gw._parse_discovery_entries(_concat_stream(), CYCLE_TAGS, gene_filter="GENEA")
    assert set(out) == {"c1", "c2"}
    assert [e.position for e in out["c1"]] == [10]
    e = out["c1"][0]
    assert e.position_one_based == 11
    assert e.caap_group == "US"
    assert e.conserved_pair == "sp2"          # ":"-prefix stripped
    assert e.trait1_aa and e.trait0_aa        # caas "V/A" split


def test_shard_layout_no_gene_column():
    out = gw._parse_discovery_entries(_shard_stream(), CYCLE_TAGS)
    assert set(out) == {"c1", "c2"}
    assert out["c2"][0].is_conserved_meta is True
    assert out["c2"][0].caap_group == "GS2"


def test_concat_layout_requires_gene_column_when_filtering():
    # a stream with no `gene` column but a gene_filter given -> nothing parsed
    out = gw._parse_discovery_entries(_shard_stream(), CYCLE_TAGS, gene_filter="GENEA")
    assert out == {}


def test_empty_stream():
    assert gw._parse_discovery_entries(io.StringIO(""), CYCLE_TAGS) == {}


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
