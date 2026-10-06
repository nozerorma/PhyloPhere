#!/usr/bin/env python3
# build_domino_network.py — Builds DOMINO's input network (network.sif) from STRING links.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
BuildDominoNetwork: filters the STRING protein-links file of one species to the edges
whose two endpoints are both in the gene background, and writes them as the network
DOMINO reads.

DOMINO has no background argument: the genes of the network file are its implicit
background. The network is therefore restricted to the analysis background here, not
left to DOMINO.

The raw STRING links and protein-info files are cached unfiltered (the same cache
convention as ensure_string_cache() in 13.AMI_analysis.Rmd), so one cache serves any
--score-threshold, which is applied at use time. Node IDs are mapped from ENSP to gene
symbol in a single hop with the protein-info file; no round trip through Ensembl gene IDs.

network.sif keeps the three columns DOMINO's own tools (slicer, run_domino_modules.py)
expect (gene A, "pp", gene B; only columns 0 and 2 are read). The combined score, which
the SIF cannot carry, goes to a sidecar table that 13.AMI_analysis.Rmd uses to scale the
width of the network edges.

Called by:  DOMINO_BUILD_NETWORK Nextflow process (domino.nf → build_domino_network.py)
Inputs:     --cleaned-background  gene symbols, one per line; both endpoints of an edge must be in it
            --string-db-dir       optional directory of pre-cached STRING files, checked before downloading
            --species, --version  STRING species taxon ID and release (files <species>.protein.<kind>.v<version>.txt.gz)
Outputs:    network.sif               geneA<TAB>pp<TAB>geneB, one line per undirected edge
            network_edge_scores.tsv   gene1, gene2, combined_score (0-1000, max over the two directions)
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import gzip
import os
import shutil
import sys
import urllib.request


# ── Constants ─────────────────────────────────────────────────────────────────

STRING_BASE = "https://stringdb-downloads.org/download"


# ── CLI ───────────────────────────────────────────────────────────────────────


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


# ── STRING files ──────────────────────────────────────────────────────────────


def ensure_cached(kind, species, version, string_db_dir, cache_dir):
    """Return the path of the cached STRING file of `kind` ('links' or 'info').

    Lookup order: the cache directory, then string_db_dir (copied into the cache), then a
    download from STRING. Raises RuntimeError when no source yields a non-empty file.
    """
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
    """Map STRING protein ID (ENSP) to preferred gene symbol.

    Reads the protein-info file: tab-separated, one header line, columns
    #string_protein_id, preferred_name, protein_size, annotation.
    """
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
    """Set of gene symbols of the background file (one per line, blank lines ignored)."""
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def build_edges(links_path, id2sym, background, score_threshold):
    """Edges of the links file as {(gene_a, gene_b): combined_score}, gene_a < gene_b.

    The links file is space-separated (protein1 protein2 combined_score) and lists each
    undirected edge once per direction. Edges are deduplicated on the sorted gene pair,
    keeping the higher score, since the two directions can differ slightly. An edge is
    kept when its score reaches score_threshold, both proteins map to different symbols
    and both symbols are in the background.
    """
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


# ── Main ──────────────────────────────────────────────────────────────────────


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
