"""validation.harness.cli — thin command line over the harness.

    python -m validation.harness.cli emit  <spec.json> --out <dir> [--quantile 0.25]
    python -m validation.harness.cli score --truth <sites.tsv> --run <run_complete_dir>

``emit``  expands a DatasetSpec into the trait matrix (Tier 1+).
``score`` scores a completed Tier 1 pipeline run against a site-level truth set
          (validated-convergent-site recovery + rank; ``tier1_adapter.py``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .phenotype import load_dataset_spec
from .trait_matrix import emit_trait_matrix
from .tier1_adapter import score_tier1


def _cmd_emit(args: argparse.Namespace) -> int:
    spec = load_dataset_spec(args.spec)
    manifest = emit_trait_matrix(spec, args.out, quantile=args.quantile)
    print(json.dumps(manifest, indent=2))
    if manifest["warnings"]:
        print("\nWARNINGS:", *manifest["warnings"], sep="\n  ", file=sys.stderr)
        return 2
    return 0


def _cmd_score(args: argparse.Namespace) -> int:
    ref_row = None
    if args.ref_row_file:
        # a fasta with the single gapped reference sequence
        txt = Path(args.ref_row_file).read_text()
        ref_row = "".join(
            ln.strip() for ln in txt.splitlines() if not ln.startswith(">")
        )
    rep = score_tier1(
        args.truth, args.run,
        ref_row=ref_row, slop=args.slop, dataset=args.dataset,
    )
    print(json.dumps(rep.to_dict(), indent=2))

    r = rep.recall_by_tier
    print("\n" + "=" * 60, file=sys.stderr)
    print(f"{rep.n_truth} truth sites | {rep.n_scored_positions} scored positions",
          file=sys.stderr)
    for tier, d in r.items():
        rp = rep.median_rank_percentile.get(tier)
        rp_s = f"{rp:.2f}" if rp == rp else "-"       # nan check
        print(f"  {tier:<14} recall {d['recovered']}/{d['n']} "
              f"= {d['recall']:.2f}  CI{d['wilson_95ci']}  median rank-pct {rp_s}",
          file=sys.stderr)
    for m, c in rep.comparator.items():
        print(f"  vs {m}: we recovered {c['we_also_recovered']}; "
              f"missed {c['they_called_we_missed']}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validation.harness.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="expand a DatasetSpec into the trait matrix")
    e.add_argument("spec", type=Path)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--quantile", type=float, default=0.25)
    e.set_defaults(func=_cmd_emit)

    s = sub.add_parser("score", help="score a Tier 1 run vs a site-level truth set")
    s.add_argument("--truth", type=Path, required=True,
                   help="truthsets/tier1/<name>.sites.tsv")
    s.add_argument("--run", type=Path, required=True,
                   help="the <trait>_complete output dir (or its parent)")
    s.add_argument("--ref-row-file", type=Path, default=None,
                   help="fasta of the gapped reference row, if the fixture carries one "
                        "(omit when alignment column == reference position, e.g. PEPC)")
    s.add_argument("--slop", type=int, default=0,
                   help="search +/- N columns for a residue-matching CAAS "
                        "(guards truth-table transcription offsets; residues must still match)")
    s.add_argument("--dataset", default=None,
                   help="key into DEFAULT_COMPARATORS, e.g. 'pepc_c4' for the ConDor calls")
    s.set_defaults(func=_cmd_score)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
