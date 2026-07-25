#!/usr/bin/env python3
# =============================================================================
# run_domino_modules.py — active module identification via DOMINO's Python API
# =============================================================================
# Calls src.core.domino.main() directly instead of shelling out to the `domino`
# CLI. This matters for exactly one reason: DOMINO's own module_threshold gate
# inside get_final_modules() computes a real Bonferroni-corrected hypergeometric
# sig_score per module (scipy.stats.hypergeom), but the reference return
# statement is `return [a[0] for a in module_sigs]` — it keeps only the module,
# throws the score away. The CLI (src/runner.py's main_domino) then only ever
# writes gene names to modules.out. There is no flag or API path that returns
# the score; the only way to recover it is to intercept get_final_modules()
# before it discards its own second tuple element.
#
# _get_final_modules_capturing() below is byte-for-byte the same statistical
# logic as the upstream get_final_modules() (verified against
# github.com/Shamir-Lab/DOMINO/blob/master/src/core/domino.py) — same
# hypergeom.sf + hypergeom.pmf tail-inclusive test, same
# module_threshold/len(putative_modules) Bonferroni correction, same sort order
# — it just also stashes the per-module sig_score instead of dropping it. It is
# monkeypatched onto the domino module before calling main(), so every other
# step (network build/caching, slice pruning, PCST optimization) runs exactly
# as upstream.
# =============================================================================

import argparse
import glob
import os
import sys

import pandas as pd
from scipy.stats import hypergeom


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


def install_capturing_get_final_modules(domino_core):
    """Monkeypatch domino_core.get_final_modules with a copy that also records
    each accepted module's sig_score (see module docstring for why)."""
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
    captured.clear()
    final_modules = domino_core.main(
        active_genes_file=list_path,
        network_file=network_file,
        slices_file=slices_file,
        slice_threshold=slice_threshold,
        module_threshold=module_threshold,
    )
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
