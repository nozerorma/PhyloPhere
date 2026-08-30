"""validation.harness.cli — thin command line over the harness.

    python -m validation.harness.cli emit <spec.json> --out <dir> [--quantile 0.25]
    python -m validation.harness.cli calibrate --run <runs/tier0> [--pcol meta_caas_boot]
    python -m validation.harness.cli score     --run <runs/tier0> [--json out.json]

``emit`` expands a DatasetSpec into the trait matrix. ``calibrate`` / ``score``
read a staged Tier 0 ``runs/`` tree via ``tier0_adapter`` — see that module for
the file/column mapping.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from . import tier0_adapter as t0
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


def _fmt(x: float) -> str:
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.3f}"


def _cmd_calibrate(args: argparse.Namespace) -> int:
    res = t0.calibrate(args.run, args.pcol, alpha=args.alpha)
    d = res.to_dict()
    if args.json:
        Path(args.json).write_text(json.dumps(d, indent=2))
    print(json.dumps(d, indent=2))
    print(
        f"\n{res.pcol}: n={res.n_pvalues} over {res.n_replicates} null context(s) | "
        f"KS D={_fmt(res.report.ks_stat)} p={_fmt(res.report.ks_pvalue)} | "
        f"mean={_fmt(res.report.mean)} | "
        f"type-I@{args.alpha}={_fmt(res.report.type1_error.get(args.alpha))}",
        file=sys.stderr,
    )
    print(f"VERDICT: {res.verdict}", file=sys.stderr)
    return 0 if res.verdict == "OK" else 3


def _cmd_score(args: argparse.Namespace) -> int:
    res = t0.score(args.run)
    d = res.to_dict()
    if args.json:
        Path(args.json).write_text(json.dumps(d, indent=2))
    print(json.dumps(d, indent=2))
    s = res.summary
    print(
        f"\n{res.n_power_replicates} power replicate(s)\n"
        f"  gene   : ROC-AUC {_fmt(s['gene_roc_auc'])}  PR-AUC {_fmt(s['gene_pr_auc'])}  "
        f"precision@k {_fmt(s['gene_precision_at_k'])}\n"
        f"  sites  : recall {_fmt(s['site_recall'])}  precision {_fmt(s['site_precision'])}  "
        f"by-mechanism {{{', '.join(f'{m}:{_fmt(v)}' for m, v in s['site_recall_by_mechanism'].items())}}}\n"
        f"  schemes: identical_aa->US {_fmt(s['identical_aa_us_recall'])}  "
        f"grouped_caap->GS {_fmt(s['grouped_caap_gs_recall'])}  "
        f"(US leakage {_fmt(s['grouped_caap_us_leakage'])})\n"
        f"  contrast: jaccard {_fmt(s['contrast_jaccard'])}  fg-precision {_fmt(s['contrast_fg_precision'])}  "
        f"pairs-recovered {_fmt(s['contrast_pairs_recovered_frac'])}",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validation.harness.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="expand a DatasetSpec into the trait matrix")
    e.add_argument("spec", type=Path)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--quantile", type=float, default=0.25)
    e.set_defaults(func=_cmd_emit)

    c = sub.add_parser("calibrate", help="null-calibration (KS uniformity) from a Tier 0 run")
    c.add_argument("--run", type=Path, required=True)
    c.add_argument("--pcol", default="meta_caas_boot", choices=sorted(t0.PCOLS))
    c.add_argument("--alpha", type=float, default=0.05)
    c.add_argument("--json", type=Path, default=None)
    c.set_defaults(func=_cmd_calibrate)

    s = sub.add_parser("score", help="planted recovery + contrast recovery from a Tier 0 run")
    s.add_argument("--run", type=Path, required=True)
    s.add_argument("--json", type=Path, default=None)
    s.set_defaults(func=_cmd_score)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
