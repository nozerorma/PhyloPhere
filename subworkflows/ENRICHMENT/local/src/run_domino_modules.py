#!/usr/bin/env python3
# run_domino_modules.py — Active module identification with DOMINO's Python API, one gene list at a time.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
RunDominoModules: finds the active modules of each gene list on the STRING network and
reports, for every module, its genes and its Bonferroni-corrected hypergeometric p-value.

src.core.domino.main() is called directly instead of the `domino` CLI because the
per-module score is otherwise lost: get_final_modules() computes it, but returns only
the modules, and the CLI writes gene names only. install_capturing_get_final_modules()
replaces that function with one that applies the same test (hypergeom.sf + hypergeom.pmf,
tail inclusive; module_threshold divided by the number of putative modules; ascending
score order) and also stores the scores. Every other step (network caching, slice
pruning, PCST optimization) is upstream DOMINO.

Called by:  DOMINO_RUN_MODULES Nextflow process (domino.nf → run_domino_modules.py)
Inputs:     --network          network.sif (build_domino_network.py)
            --slices           slices.txt (output of the `slicer` CLI)
            --gene-lists-dir   directory of *.txt active-gene files, one gene per line
Outputs:    <list>_domino_modules.tsv       node, cluster (cluster = rank by p-value)
            <list>_domino_module_stats.tsv  cluster, n_genes, genes (comma-separated), p_value, p_adj
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import glob
import os
import sys

# ── Third-party ───────────────────────────────────────────────────────────────
import pandas as pd
from scipy.stats import hypergeom


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Run DOMINO active-module detection per gene list, with p-values.")
    p.add_argument("--network", required=True, help="network.sif (built by build_domino_network.py)")
    p.add_argument("--slices", required=True, help="slices.txt (output of the `slicer` CLI)")
    p.add_argument("--gene-lists-dir", required=True,
                   help="directory of *.txt active-gene files, one per percentile-tier slice")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--slice-threshold", type=float, default=0.3,
                   help="DOMINO's own default (src/runner.py) for retaining a slice as relevant")
    p.add_argument("--module-threshold", type=float, default=0.05,
                   help="DOMINO's own default (src/runner.py) for accepting a putative module as final")
    p.add_argument("--threads", type=int, default=1,
                   help="constants.N_OF_THREADS override (upstream default is a hardcoded 40, which "
                        "oversubscribes a shared cluster node); pass task.cpus from Nextflow")
    return p.parse_args()


# ── DOMINO wrapper ────────────────────────────────────────────────────────────


def install_capturing_get_final_modules(domino_core):
    """Replace domino_core.get_final_modules with a copy that also records the scores.

    Returns the dict the copy fills on each call: sig_scores (Bonferroni-corrected score of
    every accepted module, ascending) and n_putative (number of putative modules tested).
    """
    captured = {}

    def _get_final_modules_capturing(G, G_putative_modules, module_threshold):
        module_sigs = []
        n_putative = len(G_putative_modules)
        for cur_G_module in G_putative_modules:
            pertubed_nodes_in_cc = [n for n in cur_G_module.nodes if G.nodes[n]["pertubed_node"]]
            pertubed_nodes = [n for n in G.nodes if G.nodes[n]["pertubed_node"]]
            sig_score = (
                hypergeom.sf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes), len(cur_G_module.nodes))
                + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes), len(cur_G_module.nodes))
            )
            final_module_threshold = module_threshold / n_putative if n_putative else module_threshold
            if sig_score <= final_module_threshold:
                module_sigs.append((cur_G_module, sig_score / n_putative if n_putative else sig_score))

        module_sigs = sorted(module_sigs, key=lambda a: a[1])
        captured["sig_scores"] = [s for _, s in module_sigs]
        captured["n_putative"] = n_putative
        return [m for m, _ in module_sigs]

    domino_core.get_final_modules = _get_final_modules_capturing
    return captured


def run_one_list(domino_core, captured, list_path, network_file, slices_file, slice_threshold, module_threshold):
    """Run DOMINO on one active-gene list.

    Returns (rows_modules, rows_stats): one row per gene (node, cluster) and one per module
    (cluster, n_genes, genes, p_value, p_adj). p_adj is p_value times n_putative, capped at 1.
    """
    captured.clear()
    try:
        final_modules = domino_core.main(
            active_genes_file=list_path,
            network_file=network_file,
            slices_file=slices_file,
            slice_threshold=slice_threshold,
            module_threshold=module_threshold,
        )
    except ValueError as e:
        # Some DOMINO builds raise a ValueError mentioning union_all when zero modules
        # survive modularity slicing. That is a legitimate outcome for a small or sparse
        # network, not malformed input, so it yields zero modules instead of failing the
        # whole batch of gene lists.
        if "union_all" not in str(e):
            raise
        print(f"[run_domino_modules] {os.path.basename(list_path)}: 0 modules "
              f"survived modularity slicing (network too sparse for this "
              f"gene list) -- treating as zero final modules.", file=sys.stderr)
        final_modules = []
    sig_scores = captured.get("sig_scores", [None] * len(final_modules))
    n_putative = captured.get("n_putative", len(final_modules))

    rows_modules = []
    rows_stats = []
    for i, (mod, p) in enumerate(zip(final_modules, sig_scores), start=1):
        genes = sorted(mod.nodes)
        for g in genes:
            rows_modules.append({"node": g, "cluster": i})
        p_adj = min(1.0, p * n_putative) if p is not None else None
        rows_stats.append({
            "cluster": i,
            "n_genes": len(genes),
            "genes": ",".join(genes),
            "p_value": p,
            "p_adj": p_adj,
        })
    return rows_modules, rows_stats


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    import src.constants as constants
    constants.USE_CACHE = True
    constants.N_OF_THREADS = max(1, args.threads)

    import src.core.domino as domino_core
    captured = install_capturing_get_final_modules(domino_core)

    list_paths = sorted(glob.glob(os.path.join(args.gene_lists_dir, "*.txt")))
    if not list_paths:
        raise RuntimeError(f"no *.txt gene list files found in {args.gene_lists_dir}")

    for list_path in list_paths:
        list_name = os.path.splitext(os.path.basename(list_path))[0]
        print(f"[run_domino_modules] {list_name}: running DOMINO...", file=sys.stderr)
        rows_modules, rows_stats = run_one_list(
            domino_core, captured, list_path, args.network, args.slices,
            args.slice_threshold, args.module_threshold,
        )
        print(f"[run_domino_modules] {list_name}: {len(rows_stats)} final modules "
              f"({sum(r['n_genes'] for r in rows_stats)} genes total)", file=sys.stderr)

        pd.DataFrame(rows_modules, columns=["node", "cluster"]).to_csv(
            os.path.join(args.output_dir, f"{list_name}_domino_modules.tsv"), sep="\t", index=False
        )
        pd.DataFrame(rows_stats, columns=["cluster", "n_genes", "genes", "p_value", "p_adj"]).to_csv(
            os.path.join(args.output_dir, f"{list_name}_domino_module_stats.tsv"), sep="\t", index=False
        )


if __name__ == "__main__":
    main()
