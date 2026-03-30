#!/usr/bin/env python3
# ================================ src/caas_main.py ================================

from __future__ import annotations

import argparse
import logging
import sys

from utils.logging_conf import configure_logging


def _read_tax_map(path: str) -> tuple[dict, dict]:
    import pandas as pd

    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = df.columns.str.strip()
    need = {"tax_id", "species"}
    if not need.issubset(df.columns):
        raise SystemExit("Tax map must have columns: tax_id, species")
    df["species"] = df["species"].astype(str).str.strip()
    df["tax_id"] = df["tax_id"].astype(str).str.strip()
    sp2tax = dict(zip(df["species"], df["tax_id"]))
    tax2sp = dict(zip(df["tax_id"], df["species"]))
    return sp2tax, tax2sp


def _detect_tree_format(path: str, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    p = path.lower()
    return "nexus" if (p.endswith(".nex") or p.endswith(".nexus")) else "newick"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="CAAS validation long table + internal-species site-level phylogenetic PGLS"
    )
    ap.add_argument("--trait-file", required=True)
    ap.add_argument("--species-col", default="species")
    ap.add_argument("--trait-name", required=True, help="Single trait column to analyze.")
    ap.add_argument("--caas-file", required=True)
    ap.add_argument("--alignments-dir", required=True)
    ap.add_argument("--alignments-format", default="auto")

    ap.add_argument("--q-lower", type=float, default=0.25)
    ap.add_argument("--q-upper", type=float, default=0.75)
    ap.add_argument("--extremes-method", choices=["quantile", "1sd", "2sd"], default="quantile",
                    help="Method for defining trait extremes (default: quantile)")
    ap.add_argument("--use-median-sep", action=argparse.BooleanOptionalAction, default=True,
                    help="Use median separation for quantile method (default: True)")
    ap.add_argument("--keep-ambiguous", action="store_true")
    ap.add_argument("--include-other", action="store_true")

    ap.add_argument("--tree-file", required=True)
    ap.add_argument("--tree-format", default="auto", choices=["auto", "newick", "nexus"])
    ap.add_argument("--tax-map-file", help="CSV/TSV with columns: tax_id, species")

    ap.add_argument("--pgls-out", required=True, help="Output TSV for internal-species phylogenetic tests")
    ap.add_argument("--extremes-DEBUG-out", help="Optional TSV of LOW/HIGH membership per trait/species")
    ap.add_argument("--intval-out", help="Optional TSV: per (Gene,Position,species) AA snapshot")
    ap.add_argument("--site-diag-out", help="Optional TSV with per-site diagnostics/log rows")
    ap.add_argument("--threads", type=int, default=1, help="Number of worker threads for site PGLS (default: 1)")
    ap.add_argument("--site-fdr-alpha", type=float, default=0.05,
                    help="Benjamini-Hochberg FDR alpha for per-site tests, applied per trait (default 0.05)")
    ap.add_argument("--min-per-class", type=int, default=6,
                    help="Minimum number of internal species carrying each of the TOP and BOTTOM amino-acid groups (default: 6)")

    ap.add_argument("--phylo-select", choices=["auto", "fixed"], default="fixed")
    ap.add_argument("--phylo-candidates", default="lambda")
    ap.add_argument("--phylo-bounds", choices=["wide", "phylolm"], default="wide")
    ap.add_argument("--phylo-model", choices=["bm", "lambda", "ou", "eb", "delta", "kappa"], default="lambda")
    ap.add_argument("--lambda-value", type=float, default=None,
                    help="Optional fixed lambda value. If omitted with --phylo-model lambda, lambda is fitted per site.")
    ap.add_argument("--ou-alpha", type=float, default=0.0)

    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    configure_logging(args.log_level)

    trait_name = args.trait_name.strip()
    if not trait_name:
        raise SystemExit("--trait-name must be a single non-empty trait column")
    if "," in trait_name:
        raise SystemExit("--trait-name accepts exactly one trait column; comma-separated values are not supported")

    import pandas as pd
    from caas.pgls_runner import DIAG_COLUMNS, RESULT_COLUMNS, run_caas_pgls
    from caas.validate import build_validation_table, extremes_DEBUG_table

    try:
        from Bio import Phylo
    except ImportError as e:
        raise SystemExit(f"Biopython is required to read trees but is not available: {e}")

    traits_df = pd.read_csv(args.trait_file, sep=None, engine="python")
    traits_df.columns = traits_df.columns.str.strip()
    traits_df[args.species_col] = traits_df[args.species_col].astype(str).str.strip()
    if trait_name not in traits_df.columns:
        raise SystemExit(f"Trait file missing column '{trait_name}'")
    traits_df[trait_name] = pd.to_numeric(traits_df[trait_name], errors="coerce")

    if args.extremes_DEBUG_out:
        debug_df = extremes_DEBUG_table(
            traits_df,
            args.species_col,
            trait_name,
            args.q_lower,
            args.q_upper,
            extremes_method=args.extremes_method,
            use_median_separation=args.use_median_sep,
        )
        debug_df.to_csv(args.extremes_DEBUG_out, sep="\t", index=False)
        logging.info(f"Wrote extremes membership to {args.extremes_DEBUG_out}")

    caas_df = pd.read_csv(args.caas_file, sep=None, engine="python")
    caas_df.columns = caas_df.columns.str.strip()

    sp2tax, tax2sp = (None, None)
    if args.tax_map_file:
        sp2tax, tax2sp = _read_tax_map(args.tax_map_file)

    valid_df = build_validation_table(
        caas_df=caas_df,
        traits_df=traits_df,
        species_col=args.species_col,
        trait_name=trait_name,
        align_dir=args.alignments_dir,
        align_fmt=args.alignments_format,
        q_lower=args.q_lower,
        q_upper=args.q_upper,
        keep_ambiguous=args.keep_ambiguous,
        include_other=args.include_other,
        extremes_method=args.extremes_method,
        use_median_separation=args.use_median_sep,
        DEBUG=(logging.getLogger().level <= logging.DEBUG),
        sp2tax=sp2tax,
        allowed_tax_ids=None,
    )

    if args.intval_out:
        cols = [c for c in ["Gene", "Position", "Substitution", "species", "observed_AA", "observed_dir",
                            "trait", "trait_value", "extreme_side", "q1", "q3"] if c in valid_df.columns]
        valid_df[cols].drop_duplicates().to_csv(args.intval_out, sep="\t", index=False)
        logging.info(f"Wrote intval snapshot to {args.intval_out}")

    if valid_df.empty:
        logging.warning("Validation table is empty; writing empty PGLS output.")
        pd.DataFrame(columns=RESULT_COLUMNS).to_csv(args.pgls_out, sep="\t", index=False)
        if args.site_diag_out:
            pd.DataFrame(columns=DIAG_COLUMNS).to_csv(args.site_diag_out, sep="\t", index=False)
        sys.exit(0)

    fmt = _detect_tree_format(args.tree_file, args.tree_format)
    tree = Phylo.read(args.tree_file, fmt)
    tip_set = {str(t.name) for t in tree.get_terminals()}

    select_cfg = {
        "mode": args.phylo_select,
        "candidates": [m.strip() for m in args.phylo_candidates.split(",") if m.strip()],
        "bounds": args.phylo_bounds,
        "fixed_model": args.phylo_model,
        "fixed_theta": (args.lambda_value if args.phylo_model == "lambda"
                        else (args.ou_alpha if args.phylo_model == "ou" else None)),
    }

    res_df, site_diag = run_caas_pgls(
        valid_df=valid_df,
        tree=tree,
        tip_set=tip_set,
        sp2tax=sp2tax,
        min_per_class=args.min_per_class,
        select_cfg=select_cfg,
        threads=max(1, int(args.threads)),
        fdr_alpha=float(args.site_fdr_alpha),
    )

    res_df.to_csv(args.pgls_out, sep="\t", index=False)
    logging.info(f"Wrote phylogenetic test summary to {args.pgls_out}")
    if args.site_diag_out:
        site_diag.to_csv(args.site_diag_out, sep="\t", index=False)
        logging.info(f"Wrote site diagnostics to {args.site_diag_out}")


if __name__ == "__main__":
    main()
