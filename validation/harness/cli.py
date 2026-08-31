"""validation.harness.cli — thin command line over the harness.

    python -m validation.harness.cli emit  <spec.json> --out <dir> [--quantile 0.25]

``emit`` expands a DatasetSpec into the trait matrix (Tier 1+).

The ``score`` sub-command was Tier 0's gate report; Tier 0 is demoted
(``validation/.demoted/``, D-DIR-01). Tier 1 will re-add a ``score`` command
with a site-truth adapter — see ``validation/.demoted/tier0_adapter.py`` for the
old shape and ``validation/harness/truthset.py`` for the truth-set loaders.
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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validation.harness.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="expand a DatasetSpec into the trait matrix")
    e.add_argument("spec", type=Path)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--quantile", type=float, default=0.25)
    e.set_defaults(func=_cmd_emit)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
