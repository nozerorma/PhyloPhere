"""validation.harness.cli — thin command line over the harness.

    python -m validation.harness.cli emit <spec.json> --out <dir> [--quantile 0.25]
    python -m validation.harness.cli score --truth <t.tsv> --run <dir> [--ref-row ID]   (STUB)
    python -m validation.harness.cli calibrate --run <dir> --pcol <col>                 (STUB)

``emit`` is wired. ``score`` / ``calibrate`` need the PhyloPhere-output adapter,
which is deliberately not guessed here — it should be written against a real
``runs/`` tree so the column mapping matches what the pipeline actually emits.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .phenotype import load_dataset_spec
from .trait_matrix import emit_trait_matrix


def _cmd_emit(args: argparse.Namespace) -> int:
    spec = load_dataset_spec(args.spec)
    manifest = emit_trait_matrix(spec, args.out, quantile=args.quantile)
    print(json.dumps(manifest, indent=2))
    if manifest["warnings"]:
        print("\nWARNINGS:", *manifest["warnings"], sep="\n  ", file=sys.stderr)
        return 2
    return 0


def _cmd_stub(name: str) -> int:
    print(
        f"`{name}` is not implemented yet. It needs the adapter that maps a "
        f"PhyloPhere runs/ output table to (id, score, call) / (p-value) tuples. "
        f"Write it against a real run — see validation/tier1/README.md.",
        file=sys.stderr,
    )
    return 64


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validation.harness.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="expand a DatasetSpec into the trait matrix")
    e.add_argument("spec", type=Path)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--quantile", type=float, default=0.25)
    e.set_defaults(func=_cmd_emit)

    s = sub.add_parser("score", help="(stub) score a run against a truth set")
    s.add_argument("--truth", type=Path, required=True)
    s.add_argument("--run", type=Path, required=True)
    s.add_argument("--ref-row", default=None)
    s.set_defaults(func=lambda a: _cmd_stub("score"))

    c = sub.add_parser("calibrate", help="(stub) null-calibration report from a run")
    c.add_argument("--run", type=Path, required=True)
    c.add_argument("--pcol", required=True)
    c.set_defaults(func=lambda a: _cmd_stub("calibrate"))

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
