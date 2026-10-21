#!/usr/bin/env python3
# accumulation_test.py — Per-gene exact test of CAAS accumulation.
# PhyloPhere | subworkflows/CT_ACCUMULATION/local/

"""
Accumulation test: does a gene hold more CAAS positions than expected from its share of the position universe?

Universe U
----------
U is the set of positions an observed CAAS could come from, the same universe SCORING and POSENRICH use:
the positions CAAStools tested (background.output) in the genes that survive post-processing
(cleaned_background_main.txt), kept when they score caas_score > 0 in at least one cycle and side of the
filtered permulation null (perm_pos_cycle_caas.tsv.gz), plus every observed CAAS position. A position the
null never scores cannot produce a CAAS under the permulation labelings, so it is not part of the
opportunity of a gene.

Test
----
Per direction (top, bottom, both) let M be the number of distinct positions of U, K the number of distinct
observed CAAS positions of the direction (a position detected on both sides counts once in 'both'), n_g the
number of positions of gene g in U and obs_g its observed CAAS positions. Under the null the K observed
positions are a uniform draw without replacement from U, so the count of gene g is hypergeometric:
    p_g = P(X >= obs_g),  X ~ Hypergeometric(M, n_g, K)
    MeanNull_g = K n_g / M,  SD_g = sqrt(K (n_g/M) (1 - n_g/M) (M - K) / (M - 1))
The gene length enters through n_g. Genes without observed CAAS have p = 1 and are not tested; the
Benjamini-Hochberg adjustment over the tested genes is done downstream in accum_gene_lists.nf and
scoring_compute.R. Positions are treated as independent draws, which understates the variance of genes with
many positions in U.

Called by:  CT_ACCUMULATION Nextflow process (ctacc_run.nf CT_ACCUMULATION_TEST)
Inputs:     --caas-csv (filtered_discovery.tsv: one row per gene, position and side), --background-positions
            (caastools background.output), --bg-caas (cleaned_background_main.txt), --null-cycle-caas
            (perm_pos_cycle_caas.tsv.gz)
Outputs:    --output CSV with columns Gene, ActualCount, PoolPositions, MeanNull, SD, PValue
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import os

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

POS_BITS = 32                      # position keys: gene id << 32 | position
COLUMNS = ['Gene', 'ActualCount', 'PoolPositions', 'MeanNull', 'SD', 'PValue']


# ── Universe U ────────────────────────────────────────────────────────────────


def pair_keys(gene_ids, positions):
    """One int64 key per (gene id, position)."""
    return (np.asarray(gene_ids, dtype=np.int64) << POS_BITS) + np.asarray(positions, dtype=np.int64)


def load_universe_genes(path):
    """Genes that survive post-processing (cleaned_background_main.txt: one gene per line), in file order."""
    with open(path) as f:
        genes = [line.strip() for line in f if line.strip() and line.strip() != "Gene"]
    return list(dict.fromkeys(genes))


def load_tested_keys(background_file, gene_to_id):
    """Positions CAAStools tested, as int64 keys (gene id, position) of the genes in `gene_to_id`.

    Format mirrors posenrich's build_background(): `gene<TAB>comma-separated positions`, one line per gene,
    in the coordinate space of filtered_discovery's `Position`. Genes outside `gene_to_id` (removed by
    post-processing) can never contribute an observed CAAS and are dropped.
    """
    keys = []
    with open(background_file) as f:
        for line in f:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) < 2 or parts[0] not in gene_to_id or parts[1] in ("", "NULL", "Position"):
                continue
            pos = np.array([int(p) for p in parts[1].split(",") if p.strip()], dtype=np.int64)
            keys.append(pair_keys(np.full(pos.size, gene_to_id[parts[0]]), pos))
    return np.unique(np.concatenate(keys)) if keys else np.empty(0, dtype=np.int64)


def load_null_positive_keys(null_path, gene_to_id):
    """Keys (gene id, position) with caas_score > 0 in any cycle and side of the filtered null table.

    The table is read in chunks and reduced to distinct keys chunk by chunk, so memory follows the number of
    distinct null-positive positions, not the number of (position, cycle) rows. Genes of the null outside
    `gene_to_id` are ignored.
    """
    cols = ['Gene', 'Position', 'caas_score']
    if os.path.getsize(null_path) == 0:
        return np.empty(0, dtype=np.int64)
    header = pd.read_csv(null_path, sep='\t', nrows=0).columns
    missing = [c for c in cols if c not in header]
    if missing:
        raise ValueError(f"{null_path} lacks column(s) {missing}: regenerate the CAAS permulation null (CAAS_CORE_MERGE).")
    parts = []
    for ch in pd.read_csv(null_path, sep='\t', usecols=cols, chunksize=2_000_000,
                          dtype={'Gene': str, 'Position': str, 'caas_score': str}):
        ch = ch[pd.to_numeric(ch['caas_score'], errors='coerce') > 0]
        gene = ch['Gene'].map(gene_to_id)
        pos = pd.to_numeric(ch['Position'], errors='coerce')
        ok = (gene.notna() & pos.notna()).to_numpy()
        parts.append(np.unique(pair_keys(gene[ok].to_numpy(), pos[ok].to_numpy())))
    return np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.int64)


def load_observed_keys(caas_csv, gene_to_id):
    """Observed CAAS positions by side: (keys of the top rows, keys of the bottom rows).

    filtered_discovery.tsv has one row per (gene, position, side); a position detected on both sides appears
    in both arrays. Rows of genes outside `gene_to_id` are dropped with a warning.
    """
    empty = np.empty(0, dtype=np.int64)
    if os.path.getsize(caas_csv) == 0:
        logging.warning("CAAS input file is empty: no positions passed the filter")
        return empty, empty
    df = pd.read_csv(caas_csv, sep=None, engine='python', usecols=lambda c: c in ('Gene', 'Position', 'side'))
    missing = [c for c in ('Gene', 'Position', 'side') if c not in df.columns]
    if missing:
        raise ValueError(f"{caas_csv} lacks column(s) {missing}: expected filtered_discovery.tsv.")
    gene = df['Gene'].map(gene_to_id)
    if gene.isna().any():
        logging.warning(f"{int(gene.isna().sum())} observed rows belong to genes outside the cleaned background and are dropped")
    df = df[gene.notna().to_numpy()]
    keys = pair_keys(gene[gene.notna()].to_numpy(), df['Position'].to_numpy())
    side = df['side'].astype(str).str.strip().str.lower().to_numpy()
    return np.unique(keys[side == 'top']), np.unique(keys[side == 'bottom'])


# ── Test ──────────────────────────────────────────────────────────────────────


def hypergeometric_test(obs, n_pool, M, K):
    """Per-gene p-value, null mean and null SD of the hypergeometric count.

    obs and n_pool are per-gene arrays (observed CAAS positions and positions of the gene in U), M the
    size of U and K the number of observed CAAS positions. Genes with obs = 0 get p = 1.
    """
    obs, n_pool = np.asarray(obs, dtype=np.int64), np.asarray(n_pool, dtype=np.int64)
    p = np.ones(len(obs))
    tested = obs > 0
    p[tested] = hypergeom.sf(obs[tested] - 1, M, n_pool[tested], K)
    share = n_pool / M
    mean = K * share
    sd = np.sqrt(K * share * (1 - share) * (M - K) / (M - 1)) if M > 1 else np.zeros(len(obs))
    return p, mean, sd


def run(args):
    """Run the exact accumulation test of one direction and write the results CSV.

    Reads from args: caas_csv, background_positions, bg_caas, null_cycle_caas, change_side, output.
    """
    genes = load_universe_genes(args.bg_caas)
    gene_to_id = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)
    logging.info(f"Universe (cleaned background) genes: {n_genes}")

    tested = load_tested_keys(args.background_positions, gene_to_id)
    null_positive = load_null_positive_keys(args.null_cycle_caas, gene_to_id)
    if null_positive.size == 0:
        raise ValueError(f"{args.null_cycle_caas} has no position with caas_score > 0 (a null without cycles, for "
                         f"example --caas_full_perms 0): the universe U cannot be defined.")
    top, bottom = load_observed_keys(args.caas_csv, gene_to_id)

    # U: tested and null-positive positions, plus every observed position (either side)
    observed_any = np.union1d(top, bottom)
    universe = np.union1d(np.intersect1d(tested, null_positive, assume_unique=True), observed_any)
    n_outside_tested = int(np.setdiff1d(observed_any, tested, assume_unique=True).size)
    logging.info(f"Tested positions: {tested.size}; null-positive: {null_positive.size}; |U| = {universe.size}; "
                 f"{n_outside_tested} observed positions lie outside the tested set and stay in U")
    M = int(universe.size)
    n_pool = np.bincount(universe >> POS_BITS, minlength=n_genes)

    # Observed CAAS positions of the direction: distinct positions, so a position detected on both sides
    # is one observed CAAS position in 'both'
    observed = {'top': top, 'bottom': bottom, 'both': observed_any}[args.change_side]
    K = int(observed.size)
    actual = np.bincount(observed >> POS_BITS, minlength=n_genes)
    logging.info(f"Direction '{args.change_side}': {K} distinct observed CAAS positions; "
                 f"fraction of U observed = {K / M:.4f}")

    p, mean, sd = hypergeometric_test(actual, n_pool, M, K)
    logging.info(f"Genes tested (observed CAAS > 0): {int((actual > 0).sum())} of {n_genes}")

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    pd.DataFrame({'Gene': genes, 'ActualCount': actual, 'PoolPositions': n_pool,
                  'MeanNull': mean, 'SD': sd, 'PValue': p}, columns=COLUMNS).to_csv(args.output, index=False)
    logging.info(f"Results: {args.output}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="CT_ACCUMULATION: exact per-gene accumulation test",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--caas-csv", required=True, help="filtered_discovery.tsv (one row per gene, position and side)")
    parser.add_argument("--background-positions", required=True,
                        help="caastools background.output (gene<TAB>tested msa positions)")
    parser.add_argument("--bg-caas", required=True, help="cleaned_background_main.txt (genes that survive post-processing)")
    parser.add_argument("--null-cycle-caas", required=True,
                        help="perm_pos_cycle_caas.tsv.gz from CAAS_CORE_MERGE, the filtered null that SCORING also uses")
    parser.add_argument("--change-side", default="both", choices=["top", "bottom", "both"],
                        help="Phenotype direction tested (default: both, each position counted once)")
    parser.add_argument("--output", required=True, help="Results CSV")
    parser.add_argument("--log-level", default="INFO", choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format='%(asctime)s - %(levelname)s - %(message)s')
    run(args)


if __name__ == "__main__":
    main()
