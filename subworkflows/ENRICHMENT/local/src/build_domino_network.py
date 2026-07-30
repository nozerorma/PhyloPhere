#!/usr/bin/env python3
# =============================================================================
# build_domino_network.py — builds DOMINO's input network (network.sif)
# =============================================================================
# DOMINO has no background/universe concept (confirmed from src/core/domino.py's
# main()): whatever is in the network file IS the implicit background. So the
# network itself has to be built and filtered here, not left to DOMINO.
#
# STRING v12.0's raw links file is cached once, unfiltered, exactly like
# 13.AMI_analysis.Rmd's ensure_string_cache() does for the R/STRINGdb path —
# score_threshold is applied here at use-time so the same cache serves any
# threshold. Node IDs are mapped ENSP -> gene symbol in a single hop via
# protein.info.v12.0.txt.gz; do NOT round-trip through Ensembl gene IDs, that
# loses information STRING already gives directly.
#
# Output is a plain 3-column SIF (geneA\tpp\tgeneB) — this is exactly what
# DOMINO's own build_network() expects (src/core/network_builder.py reads
# columns 0 and 2 of a tab-separated file, ignoring column 1). combined_score
# is read only to filter edges here and used to be discarded once an edge
# passed threshold, leaving 13.AMI_analysis.Rmd's plot_cluster_network() with
# no score to vary edge width/opacity by (every edge rendered identically
# thick). Also written to a sidecar network_edge_scores.tsv (gene1, gene2,
# combined_score) instead — network.sif itself stays untouched/3-column since
# DOMINO's own tools (slicer, run_domino_modules.py) expect exactly that shape.
# =============================================================================

import argparse
import gzip
import os
import shutil
import sys
import urllib.request

STRING_BASE = "https://stringdb-downloads.org/download"


def parse_args():
    p = argparse.ArgumentParser(description="Build DOMINO's network.sif from STRING v12.0 links.")
    p.add_argument("--species", type=int, default=9606)
    p.add_argument("--version", default="12.0")
    p.add_argument("--string-db-dir", default=None,
                   help="pre-cached STRING files directory, checked before downloading (same convention as string_db_dir elsewhere in this pipeline)")
    p.add_argument("--cache-dir", default="string_cache",
                   help="where raw STRING files are cached/downloaded to")
    p.add_argument("--cleaned-background", required=True,
                   help="gene symbol list (one per line) restricting both edge endpoints")
    p.add_argument("--score-threshold", type=int, default=700,
                   help="minimum STRING combined_score to keep an edge (0-1000 scale)")
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def ensure_cached(kind, species, version, string_db_dir, cache_dir):
    """kind is 'links' or 'info'. Returns the local path to the (possibly just-downloaded) file."""
    fname = f"{species}.protein.{kind}.v{version}.txt.gz"
    dest = os.path.join(cache_dir, fname)
    os.makedirs(cache_dir, exist_ok=True)

    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest

    if string_db_dir:
        local_src = os.path.join(string_db_dir, fname)
        if os.path.exists(local_src) and os.path.getsize(local_src) > 0:
            print(f"[build_domino_network] using pre-cached {local_src}", file=sys.stderr)
            shutil.copy(local_src, dest)
            return dest

    rel = f"protein.{kind}.v{version}/{fname}"
    url = f"{STRING_BASE}/{rel}"
    print(f"[build_domino_network] downloading {url}", file=sys.stderr)
    urllib.request.urlretrieve(url, dest)
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        raise RuntimeError(f"failed to obtain {fname} (from {string_db_dir} or {url})")
    return dest


def load_id2symbol(info_path):
    """protein.info.v12.0.txt.gz: #string_protein_id  preferred_name  protein_size  annotation (tab-separated)."""
    id2sym = {}
    with gzip.open(info_path, "rt") as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            id2sym[parts[0]] = parts[1]
    return id2sym


def load_background(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def build_edges(links_path, id2sym, background, score_threshold):
    """STRING links file: 'protein1 protein2 combined_score', space-separated,
    each undirected edge listed once per direction — dedup via a sorted-pair
    dict keyed on the gene pair, keeping the max score seen for that pair
    (the two directions can carry slightly different values in STRING's
    links file; max is the conservative/informative choice for display)."""
    edge_scores = {}
    n_lines = 0
    with gzip.open(links_path, "rt") as f:
        next(f)  # header
        for line in f:
            n_lines += 1
            parts = line.split()
            if len(parts) != 3:
                continue
            p1, p2, score_str = parts
            try:
                score = int(score_str)
            except ValueError:
                continue
            if score < score_threshold:
                continue
            s1 = id2sym.get(p1)
            s2 = id2sym.get(p2)
            if s1 is None or s2 is None or s1 == s2:
                continue
            if s1 not in background or s2 not in background:
                continue
            key = (s1, s2) if s1 < s2 else (s2, s1)
            edge_scores[key] = max(score, edge_scores.get(key, 0))
    print(f"[build_domino_network] scanned {n_lines} links, kept {len(edge_scores)} edges "
          f"(score >= {score_threshold}, both endpoints in background)", file=sys.stderr)
    return edge_scores


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    links_path = ensure_cached("links", args.species, args.version, args.string_db_dir, args.cache_dir)
    info_path = ensure_cached("info", args.species, args.version, args.string_db_dir, args.cache_dir)

    id2sym = load_id2symbol(info_path)
    background = load_background(args.cleaned_background)
    print(f"[build_domino_network] background genes: {len(background)}, "
          f"STRING IDs with a symbol: {len(id2sym)}", file=sys.stderr)

    edge_scores = build_edges(links_path, id2sym, background, args.score_threshold)
    if len(edge_scores) == 0:
        raise RuntimeError("no edges survived filtering — check score_threshold and cleaned_background overlap with STRING")

    nodes_covered = {g for pair in edge_scores for g in pair}
    print(f"[build_domino_network] network covers {len(nodes_covered)}/{len(background)} "
          f"background genes ({100 * len(nodes_covered) / len(background):.1f}%)", file=sys.stderr)

    sif_path = os.path.join(args.output_dir, "network.sif")
    with open(sif_path, "w") as f:
        for a, b in edge_scores:
            f.write(f"{a}\tpp\t{b}\n")
    print(f"[build_domino_network] wrote {sif_path} ({len(edge_scores)} edges, {len(nodes_covered)} nodes)", file=sys.stderr)

    scores_path = os.path.join(args.output_dir, "network_edge_scores.tsv")
    with open(scores_path, "w") as f:
        f.write("gene1\tgene2\tcombined_score\n")
        for (a, b), score in edge_scores.items():
            f.write(f"{a}\t{b}\t{score}\n")
    print(f"[build_domino_network] wrote {scores_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
