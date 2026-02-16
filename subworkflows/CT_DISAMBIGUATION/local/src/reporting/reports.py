    #!/usr/bin/env python3

"""Configuration for report generation."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any, cast
import pandas as pd
from dataclasses import dataclass


from src.reporting.disambiguation_writers import (
    write_caas_convergence_csvs,
    export_from_db,
)
from src.reporting.disambiguation_json import (
    export_aggregated_convergence_json,
    export_gene_summaries_json,
)
from src.reporting.gene_lists import export_gene_lists

# Optional: if you want plots as part of the report run
from src.plots.bulk_plots import generate_bulk_plots

logger = logging.getLogger(__name__)


DEBUG_JSON_FIELDS = {
    "node_mapping",
    "node_state_details",
    "node_posteriors",
    "site_posteriors",
    "tree_plot",
    "tip_details",
    "diagnostics",
    "node_species",
    "node_mapping_details",
}

@dataclass
class ReportConfig:
    gene_name: str
    output_dir: Path
    include_plots: bool = True
    embed_images: bool = False
    show_detailed_stats: bool = True
    include_non_significant: bool = False


def load_master_csv(path: Path) -> pd.DataFrame:
    if not path or not path.exists():
        return pd.DataFrame()
    # Your load_df helper might already exist elsewhere;
    # keep this simple here to reduce cross-module coupling.
    return pd.read_csv(path)


def orchestrate(
    *,
    results_dir: Path,
    caas_csv: Optional[Path] = None,
    ensembl_csv: Optional[Path] = None,
    db_path: Optional[Path] = None,
    trait_pairs_json: Optional[Path] = None,
    include_non_significant: bool = False,
    n_gene_trees: int = 5,
    run_plots: bool = True,
    run_json: bool = True,
) -> Dict[str, Any]:
    """
    High-level orchestrator.

    - If db_path is provided, exports master CSV from DB-backed aggregation output.
    - Otherwise uses provided caas_csv.
    - Then emits gene lists + JSON summaries.
    - Optionally generates bulk plots.
    """

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------
    # 1) Source of truth handling
    # ---------------------------
    if db_path:
        # trait_pairs is required by your DB export helper. trait_pairs_json is
        # an optional Path to JSON mapping; load it and normalize keys to ints.
        trait_pairs = {}
        if trait_pairs_json and Path(trait_pairs_json).exists():
            import json
            try:
                raw = json.loads(Path(trait_pairs_json).read_text(encoding="utf-8"))
                # Convert string keys (from JSON) to int keys expected by exporters
                trait_pairs = {
                    int(k): v for k, v in raw.items()
                }
            except Exception:
                # Fallback: attempt to use raw JSON as is
                try:
                    trait_pairs = json.loads(Path(trait_pairs_json).read_text(encoding="utf-8"))
                except Exception:
                    trait_pairs = {}

        caas_files, summary_json = export_from_db(
            db_path=db_path,
            output_dir=results_dir,
            trait_pairs=trait_pairs,
        )
        caas_csv = results_dir / "caas_convergence_master.csv"
        logger.info("DB export complete: %s", summary_json)

    if not caas_csv:
        raise FileNotFoundError("No CAAS master CSV provided and no DB export requested.")

    # ---------------------------
    # 2) Load master results
    # ---------------------------
    df = load_master_csv(caas_csv)
    if df.empty:
        logger.warning("Master CAAS CSV is empty: %s", caas_csv)
        return {"status": "empty"}

    # ---------------------------
    # 3) JSON exports
    # ---------------------------
    # Convert df rows to dicts for JSON exporters
    results = df.to_dict(orient="records")

    # -----------------------------------------------------------------------
    # Ensure canonical CSV export for consumers: CAAS master
    # -----------------------------------------------------------------------
    try:
        # Write CAAS master CSV (and no-change debug CSV) from results
        write_caas_convergence_csvs(results, results_dir)
        caas_csv = results_dir / "caas_convergence_master.csv"
        logger.info("Wrote CAAS convergence CSVs to %s", results_dir)
    except Exception as exc:  # pragma: no cover - best-effort writer
        logger.warning("Failed to write CAAS convergence CSVs: %s", exc)

    # ---------------------------
    # JSON summaries
    # ---------------------------
    # Optional JSON exports for downstream consumers
    json_out = None
    aggregated_json = None
    per_gene_jsons = []
    if run_json:
        json_out = results_dir / "json_summaries"
        json_out.mkdir(exist_ok=True)

        aggregated_json = export_aggregated_convergence_json(
            cast(List[Dict[str, Any]], results), results_dir / "caas_convergence_summary.json"
        )
        per_gene_jsons = export_gene_summaries_json(
            cast(List[Dict[str, Any]], results), json_out
        )

    # ---------------------------
    # 4) Gene list exports
    # ---------------------------
    gene_lists_out = export_gene_lists(
        results=cast(List[Dict[str, Any]], results),
        output_root=results_dir,
    )

    # ---------------------------
    # 5) Optional plots
    # ---------------------------
    plots_out = None
    if run_plots:
        # Delegate plotting to your already-refactored bulk module
        plots_out = results_dir / "plots_bulk"
        generate_bulk_plots(
            caas_csv=caas_csv,
            ensembl_csv=ensembl_csv,
            output_dir=plots_out,
            include_non_significant=include_non_significant,
            n_gene_trees=n_gene_trees,
        )

    return {
        "status": "ok",
        "caas_csv": str(caas_csv),
        "aggregated_json": str(aggregated_json),
        "per_gene_json_count": len(per_gene_jsons),
        "gene_lists": gene_lists_out,
        "plots_dir": str(plots_out) if plots_out else None,
    }


def main():
    parser = argparse.ArgumentParser(description="CAAS reporting orchestrator")
    parser.add_argument("--results-dir", type=Path, required=True)

    # Choose ONE source:
    parser.add_argument("--db", type=Path, help="Aggregation SQLite DB (preferred)")
    parser.add_argument("--trait-pairs-json", type=Path, help="Trait pairs mapping for DB export")

    parser.add_argument("--caas-csv", type=Path, help="Master CAAS CSV if not using --db")

    parser.add_argument("--ensembl-csv", type=Path, help="Ensembl annotations CSV for plots")

    parser.add_argument("--include-non-significant", action="store_true")
    parser.add_argument("--n-gene-trees", type=int, default=5)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--debug-json-fields", nargs="*", default=[], help="Additional debug JSON fields to include")

    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()


    orchestrate(
        results_dir=args.results_dir,
        caas_csv=args.caas_csv,
        ensembl_csv=args.ensembl_csv,
        db_path=args.db,
        trait_pairs_json=args.trait_pairs_json,
        include_non_significant=args.include_non_significant,
        n_gene_trees=args.n_gene_trees,
        run_plots=not args.no_plots,
        run_json=args.debug_json_fields,
    )


if __name__ == "__main__":
    main()
