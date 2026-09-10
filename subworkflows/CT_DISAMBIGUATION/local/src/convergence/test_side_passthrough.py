#!/usr/bin/env python3
"""`side` is the first-class direction key (top / bottom / none). T4b retired the
change_top/change_bottom/change_side triplet; this locks in that:
  * ConvergenceResult / PositionAxes carry a `side` field;
  * convert_convergence_result_to_dict copies it verbatim;
  * the flat-CSV header (_generate_dynamic_fields) lists `side`.

Run:  python -m pytest test_side_passthrough.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.data.models import ConvergenceResult  # noqa: E402
from src.convergence.disambiguate_single import PositionAxes  # noqa: E402
from src.utils.gene_wrapper import convert_convergence_result_to_dict  # noqa: E402
from src.reporting.disambiguation_writers import _generate_dynamic_fields  # noqa: E402


def _mk(**kw) -> ConvergenceResult:
    base = dict(gene="G", position=10, tag="t", caas="A/B", ancestral="A",
                derived="B", convergence_type="convergent")
    base.update(kw)
    return ConvergenceResult(**base)


def test_model_has_side_field():
    r = _mk()
    assert r.side == "none"
    pa = PositionAxes(position=1, caap_group="US", asr_path_score=0.0)
    assert pa.side == "none"
    assert PositionAxes(1, "US", 0.0, "top").side == "top"


def test_dict_copies_side_verbatim():
    for s in ("top", "bottom", "none"):
        d = convert_convergence_result_to_dict(_mk(side=s), None)
        assert d["side"] == s


def test_dict_missing_side_defaults_none():
    """An attribute bag that predates the `side` field still gets a sane value
    (the DB round-trip re-hydrates results as SimpleNamespace)."""
    from types import SimpleNamespace
    r = SimpleNamespace(gene="G", position=10, tag="t", caas="A/B")
    d = convert_convergence_result_to_dict(r, None)
    assert d["side"] == "none"


def test_flat_header_lists_side():
    fields = _generate_dynamic_fields(max_pairs=2)
    assert "side" in fields
    assert "change_side" not in fields
    assert "change_top" not in fields
    assert "change_bottom" not in fields


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
