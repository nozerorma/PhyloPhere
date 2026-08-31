"""validation.harness.cli — thin command line over the harness.

    python -m validation.harness.cli emit  <spec.json> --out <dir> [--quantile 0.25]
    python -m validation.harness.cli score --run <runs/tier0> [--json out.json]

``emit`` expands a DatasetSpec into the trait matrix (Tier 1+). ``score`` reads a
staged Tier 0 ``runs/`` tree via ``tier0_adapter`` and prints the gate report:
prioritisation, null-vs-power separation, site recovery by mechanism/scheme, and
contrast recovery. There is no ``calibrate`` — the pipeline's permulation
p-values are foreground-specificity scores, not Uniform(0,1) p-values, so a KS
uniformity test is the wrong instrument (see project memory).
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


def _fmt(x) -> str:
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.3f}"


def _cmd_emit(args: argparse.Namespace) -> int:
    spec = load_dataset_spec(args.spec)
    manifest = emit_trait_matrix(spec, args.out, quantile=args.quantile)
    print(json.dumps(manifest, indent=2))
    if manifest["warnings"]:
        print("\nWARNINGS:", *manifest["warnings"], sep="\n  ", file=sys.stderr)
        return 2
    return 0


def _cmd_score(args: argparse.Namespace) -> int:
    res = t0.score(args.run)
    d = res.to_dict()
    if args.json:
        Path(args.json).write_text(json.dumps(d, indent=2))
    print(json.dumps(d, indent=2))

    out = [f"\nVERDICT: {res.verdict}   ({res.n_power_replicates} power replicate(s))"]
    for sep in res.separation:
        out.append(
            f"  [{sep.archetype}] null vs power (detected-CAAS count): AUC {_fmt(sep.auc_npos)}  "
            f"planted p50 {_fmt(sep.planted_npos_p50)}  vs  null p95 {_fmt(sep.null_npos_p95)} "
            f"(max {_fmt(sep.null_npos_max)})  -> {'separated' if sep.separated else 'NOT separated'}"
            f"   [caas_score AUC {_fmt(sep.auc_caas_score)}]"
        )
    for a, v in res.summary.items():
        out.append(
            f"  [{a}] prioritisation: precision@k(n_pos) {_fmt(v['gene_precision_at_k_npos'])}  "
            f"precision@k(caas_score) {_fmt(v['gene_precision_at_k'])}  "
            f"rank-pctile {_fmt(v['gene_planted_rank_pctile_p50'])}  "
            f"slice_global25 {_fmt(v['planted_in_slice_global25'])}  "
            f"slice_global5 {_fmt(v['planted_in_slice_global5'])}\n"
            f"  [{a}] sites: recall {_fmt(v['site_recall'])}  "
            f"identical_aa->US {_fmt(v['identical_aa_us_recall'])}  "
            f"grouped_caap->GS {_fmt(v['grouped_caap_gs_recall'])} (US also {_fmt(v['grouped_caap_us_leakage'])})\n"
            f"  [{a}] contrast: jaccard {_fmt(v['contrast_jaccard'])}  "
            f"pairs recovered {_fmt(v['contrast_pairs_recovered_frac'])}"
        )
    print("\n".join(out), file=sys.stderr)
    return 0 if res.verdict == "PASS" else 3


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validation.harness.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="expand a DatasetSpec into the trait matrix")
    e.add_argument("spec", type=Path)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--quantile", type=float, default=0.25)
    e.set_defaults(func=_cmd_emit)

    s = sub.add_parser("score", help="Tier 0 gate report from a staged runs/ tree")
    s.add_argument("--run", type=Path, required=True)
    s.add_argument("--json", type=Path, default=None)
    s.set_defaults(func=_cmd_score)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
