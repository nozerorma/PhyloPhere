#!/usr/bin/env python3
"""T2a: `side` is a first-class column that passes `change_side` through
unchanged (no "both" split, no cardinality change until T3b).

Covers the three emission points touched in T2a:
  * ConvergenceResult / PositionAxes carry a `side` field;
  * convert_convergence_result_to_dict copies it (falling back to change_side);
  * the flat-CSV header (_generate_dynamic_fields) lists it right after
    change_side.

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
    pa = PositionAxes(position=1, caap_group="US", asr_path_score=0.0,
                      change_top="no_change", change_bottom="no_change")
    assert pa.side == "none"
    assert PositionAxes(1, "US", 0.0, "convergent", "no_change",
                        "top").side == "top"


def test_dict_copies_side_verbatim():
    for cs in ("top", "bottom", "both", "none"):
        d = convert_convergence_result_to_dict(_mk(change_side=cs, side=cs), None)
        assert d["side"] == cs == d["change_side"]


def test_dict_falls_back_to_change_side():
    """An attribute bag that predates the `side` field still gets a sane value
    (the DB round-trip re-hydrates results as SimpleNamespace)."""
    from types import SimpleNamespace
    r = SimpleNamespace(gene="G", position=10, tag="t", caas="A/B",
                        change_side="bottom", change_top="no_change",
                        change_bottom="no_change")
    d = convert_convergence_result_to_dict(r, None)
    assert d["side"] == "bottom"


def test_flat_header_lists_side_after_change_side():
    fields = _generate_dynamic_fields(max_pairs=2)
    assert "side" in fields
    assert fields.index("side") == fields.index("change_side") + 1


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
