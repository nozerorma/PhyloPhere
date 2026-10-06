#!/usr/bin/env python3
# posenrich_enrich.py — Position-level gene-set enrichment by path-sum permulation.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
PosenrichEnrich: tests whether the CAAS scores of the alignment positions of a gene set
exceed what the CAAS permulation null produces, for every gene set of the GMT files and
of the characterization layers, in three directions (global, top, bottom).

Statistic. The score of a term is the sum of the observed CAAS scores of its positions,
T_obs = sum_{p in term} s_p. Most background positions score 0 and add nothing, so the sum
is magnitude-weighted and needs no cutoff on the score. A fixed-cutoff test (Fisher) over
a background of this size would call biologically trivial deviations significant and would
discard the magnitude.

Null. The null is the real permulation cycles of the CAAS null (perm_pos_cycle_caas.tsv.gz),
the same preference fcs_enrich.R's fcs_run_permulation gives the FCS path-sum test: each
cycle gives one null sum per term, and p = (1 + number of null sums >= T_obs) / (n_cycles + 1).
perm_nes is (T_obs - mean) / sd of the null sums. A label shuffle of the nonzero scores
over the background is run only with --allow-label-shuffle, because it treats the scores
as independent draws and ignores the phylogenetic dependence the CAAS null accounts for
(anticonservative, exploratory). Without a null, p_value, p_adj, perm_nes, null_mean and
null_sd are NA and nothing is significant; the observed sums are still written.

p_adj is the BH adjustment within each (ranking, database). A term is significant (sig)
when p_adj < --padj-thr and perm_nes > 0.

Background. The positions tested by caastools (--background) restricted to the genes of
--universe, plus any scored position. For the cosmic_orthogroups and pai3d_orthogroups
databases it is further restricted to the genes the source database can annotate
(--cosmic-coverage, --pai3d-coverage from build_position_gmt.py).

Called by:  POSENRICH_RUN, POSENRICH_RUN_BATCHED Nextflow processes (posenrich.nf → posenrich_enrich.py)
Inputs:     --obs-scores        position_scores.tsv (Gene, Position, CAAS_score, side)
            --gmt-dir           directory of *.gmt files (position IDs Gene:Position)
            --characterization  characterization_layers.tsv, optional (name, description, members)
            --universe, --background   gene universe and caastools background.output (gene, tested positions)
            --caas-null-prepped | --caas-cycle-null   the CAAS null (prepped pickle takes precedence)
Outputs:    posenrich_characterization.tsv  ranking, database, pathway, description, layer_size,
                n_pos_with_score, obs_sum, null_mean, null_sd, perm_nes, p_value, p_adj,
                direction (enriched or depleted by the sign of perm_nes; empty without a null),
                background_n, pct_<flag> (one per flag_* column of --annot-file), n_scored, sig
            posenrich_leading_edge.tsv      gene:position driver members of the significant terms
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os
import sys
import glob
import pickle
import argparse

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import scipy.sparse as sp


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Position-level Path Sum Permulation enrichment.")
    p.add_argument("--obs-scores", required=True,
                   help="position_scores.tsv (Gene, Position, CAAS_score, side)")
    p.add_argument("--gmt-dir", required=True)
    p.add_argument("--characterization", default=None,
                   help="characterization_layers.tsv (broad functional layers)")
    p.add_argument("--annot-file", default=None,
                   help="SCORING's fcs_stats.tsv (gene + flag_* columns) - cross-module "
                        "corroboration flags (gate_sig/fade/rer/accum), reported "
                        "as the %% of distinct genes in each overlap carrying each flag")
    p.add_argument("--universe", required=True,
                   help="cleaned_background_main.txt gene list (postproc-surviving)")
    p.add_argument("--background", required=True,
                   help="caastools background.output (gene<TAB>tested positions); "
                        "restricted to --universe genes = the position background.")
    p.add_argument("--cosmic-coverage", default=None,
                   help="cosmic_coverage_genes.txt from build_position_gmt.py - genes "
                        "COSMIC itself could annotate; restricts the background used "
                        "for cosmic_orthogroups")
    p.add_argument("--pai3d-coverage", default=None,
                   help="pai3d_coverage_genes.txt from build_position_gmt.py - genes "
                        "PAI3D itself could annotate; restricts the background used "
                        "for pai3d_orthogroups")
    p.add_argument("--caas-cycle-null", default=None,
                   help="perm_pos_cycle_caas.tsv.gz from CAAS_CORE_MERGE (Gene, Position, "
                        "side, cycle, caas_score, n_schemes) - the CAAS permulation null. When "
                        "supplied, its real cycles are the null for p_value/p_adj/perm_nes, "
                        "mirroring fcs_enrich.R's fcs_run_permulation null_mat preference. "
                        "Omit or pass a NO_FILE* sentinel for no null (those values are NA "
                        "unless --allow-label-shuffle). Ignored when "
                        "--caas-null-prepped is given.")
    p.add_argument("--caas-null-prepped", default=None,
                   help="caas_null_prepped.pkl from POSENRICH_PREP_NULL (posenrich_prep_caas_null.py) "
                        "- the same CAAS permulation null as --caas-cycle-null, already parsed and "
                        "split by direction once for the whole batched run instead of once per batch "
                        "task. Takes precedence over --caas-cycle-null when both are given. Omit or "
                        "pass a NO_FILE* sentinel for no null (those values are NA unless "
                        "--allow-label-shuffle).")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-size", type=int, default=5,
                   help="min positions per set in background (GMT sources only)")
    p.add_argument("--max-size", type=int, default=0,
                   help="max positions per set in background (0 = no cap; GMT sources only)")
    p.add_argument("--n-perms", type=int, default=100000,
                   help="number of label permutations of the label-shuffle test, used only with --allow-label-shuffle "
                        "(the permulation test reads its draws from the CAAS null)")
    p.add_argument("--allow-label-shuffle", action="store_true",
                   help="when no CAAS permulation null is supplied, run a private label shuffle instead of leaving the "
                        "null-based values NA. It ignores the phylogeny (anticonservative) and is exploratory.")
    p.add_argument("--perm-chunk-size", type=int, default=1000,
                   help="permutations materialized at once as a dense (n_terms x chunk) "
                        "array before being folded into running sum/sumsq/count accumulators "
                        "(default 1000). Peak memory scales with this, not with --n-perms -- "
                        "see run_permulation_for_terms' docstring.")
    p.add_argument("--seed", type=int, default=1998,
                   help="random seed for permulations (default 1998)")
    p.add_argument("--padj-thr", type=float, default=0.15,
                   help="BH-adjusted p-value significance threshold")
    p.add_argument("--position-lists-dir", required=False, default=None,
                   help="Accepted for backward compatibility (unused in continuous permulation)")
    p.add_argument("--char-fracs", default=None, help="Accepted for backward compatibility (unused)")
    p.add_argument("--fold-thr", type=float, default=1.5, help="Accepted for backward compatibility (unused)")
    return p.parse_args()


# ── Universe and background ───────────────────────────────────────────────────


def load_universe_genes(path):
    """Set of gene symbols of a one-per-line file, ignoring blank lines and a Gene header."""
    genes = set()
    with open(path) as f:
        for line in f:
            g = line.strip()
            if g and g != "Gene":
                genes.add(g)
    return genes


def build_background(background_file, universe_genes):
    """Position IDs (Gene:Position) that caastools tested, restricted to universe genes.

    background_file is the caastools background.output (gene, comma-separated positions,
    tab-separated). Exits with an error when it is absent, because without it the
    background would be undefined.
    """
    if not background_file or background_file.startswith("NO_FILE") or not os.path.exists(background_file):
        sys.exit(
            f"[posenrich] ERROR: background file is required but was not supplied or "
            f"does not exist: {background_file!r}."
        )
    bg = []
    with open(background_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene, positions = parts[0], parts[1]
            if gene == "Gene" or gene not in universe_genes:
                continue
            if positions in ("", "NULL", "Position"):
                continue
            for p in positions.split(","):
                p = p.strip()
                if p:
                    bg.append(f"{gene}:{p}")
    return bg


# ── Gene sets ─────────────────────────────────────────────────────────────────


def load_gmts(gmt_dir):
    """Read every non-empty *.gmt of gmt_dir into {database: (terms, descriptions)}.

    terms maps the set name to its member position IDs (GMT columns 3 onward); the
    description is GMT column 2, or the set name when empty. The database is the file name.
    """
    gmts = {}
    for f in sorted(glob.glob(os.path.join(gmt_dir, "*.gmt"))):
        if os.path.getsize(f) == 0:
            continue
        db = os.path.basename(f)[:-4]
        terms = {}
        descs = {}
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                terms[parts[0]] = parts[2:]
                descs[parts[0]] = parts[1] if parts[1] else parts[0]
        if terms:
            gmts[db] = (terms, descs)
    return gmts


def read_charset(path):
    """Characterization layers as {name: (description, set of members)}; empty when absent."""
    layers = {}
    if not path or not os.path.exists(path):
        return layers
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            layers[parts[0]] = (parts[1], set(parts[2:]))
    return layers


def load_annot(path):
    """Cross-module corroboration flags per gene, from SCORING's fcs_stats.tsv.

    Returns ({gene: {flag_*: bool}}, flag column names); empty when the file is absent or
    has no gene or flag_* column.
    """
    if not path or not os.path.exists(path):
        return {}, []
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        return {}, []
    flag_cols = [c for c in df.columns if c.startswith("flag_")]
    if not flag_cols:
        return {}, []
    annot = {
        row["gene"]: {c: bool(row[c]) for c in flag_cols}
        for _, row in df.iterrows()
    }
    return annot, flag_cols


# ── Observed scores and CAAS null ─────────────────────────────────────────────


def direction_rows(df, direction):
    """One row per pos_id for a direction (global, top or bottom).

    "top" and "bottom" keep the rows of that side; "global" keeps the better side (highest
    CAAS_score), the same max-over-sides reduction SCORING applies to its undirected axis.
    """
    if direction == "global":
        sub = df
    elif direction == "top":
        sub = df[df["side"] == "top"]
    else:
        sub = df[df["side"] == "bottom"]
    if sub.empty:
        return sub
    return sub.loc[sub.groupby("pos_id", sort=False)["CAAS_score"].idxmax()]


def collapse_null_sides(sub):
    """Reduce a (pos_id, side, cycle, score) null subset to one score per (pos_id, cycle).

    The score is the maximum over sides, as direction_rows() does for the observed scores.
    """
    return (sub.groupby(["pos_id", "cycle"], observed=True, sort=False)["score"]
               .max().reset_index())


def load_caas_cycle_null(path):
    """Load perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, caas_score, n_schemes).

    Returns (long_df, all_cycle_levels). long_df has columns pos_id, side, cycle and score,
    where score is the caas_score of the position on that side in that null cycle (the
    units of CAAS_score in position_scores.tsv; 0 when unscored). all_cycle_levels is every
    distinct cycle of the file regardless of side, so a cycle without hits on one side is
    still a null draw that contributes 0 to the term sums of that side.

    Returns (None, None) when the path is absent, a NO_FILE* sentinel, missing or empty:
    callers treat it as "no CAAS null available".
    """
    if not path or os.path.basename(path).startswith("NO_FILE") or not os.path.exists(path):
        return None, None
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    if "caas_score" not in header:
        raise ValueError(f"{path} has no caas_score column: it predates the shared position score. "
                         "Regenerate the CAAS permulation null.")
    df = pd.read_csv(path, sep="\t", usecols=["Gene", "Position", "side", "cycle", "caas_score"],
                     float_precision="round_trip")
    if df.empty:
        return None, None
    df["pos_id"] = df["Gene"].astype(str) + ":" + df["Position"].astype(str)
    df["score"] = df["caas_score"].fillna(0.0)
    all_cycle_levels = np.sort(df["cycle"].unique())
    return df[["pos_id", "side", "cycle", "score"]], all_cycle_levels


def load_prepped_caas_null(path):
    """Load caas_null_prepped.pkl (posenrich_prep_caas_null.py).

    It holds the null of load_caas_cycle_null() already split by direction and reduced to
    one score per (pos_id, cycle). Returns (by_direction, all_cycle_levels), where
    by_direction is {"global"/"top"/"bottom": DataFrame[pos_id, cycle, score]}, equal to
    null_direction_subset(long_df, direction) for each direction.

    Returns (None, None) when the path is absent, a NO_FILE* sentinel, missing, or the
    artifact is empty (same contract as load_caas_cycle_null).
    """
    if not path or os.path.basename(path).startswith("NO_FILE") or not os.path.exists(path):
        return None, None
    with open(path, "rb") as fh:
        prepped = pickle.load(fh)
    if prepped.get("cycle_levels") is None:
        return None, None
    by_direction = {d: prepped[d] for d in ("global", "top", "bottom")}
    return by_direction, prepped["cycle_levels"]


def null_direction_subset(long_df, direction):
    """Direction-filtered long_df of load_caas_cycle_null, one score per (pos_id, cycle).

    Sides are handled as in direction_rows(): "global" takes the maximum over sides.
    """
    if direction == "global":
        sub = long_df
    elif direction == "top":
        sub = long_df[long_df["side"] == "top"]
    else:
        sub = long_df[long_df["side"] == "bottom"]
    return collapse_null_sides(sub)


def caas_null_term_sums(M_mat, bg_idx_map, N, null_sub, all_cycle_levels):
    """Term sums under the CAAS permulation null, as an (n_terms x n_cycles) array.

    M_mat is the term indicator matrix (n_terms x N) shared with the label-shuffle null.
    The columns span all_cycle_levels, every cycle of the file, not only cycles with a
    nonzero row in this direction and background, so a cycle without hits is a null draw
    contributing 0 instead of being dropped. Returns None when no CAAS null was supplied;
    the caller then leaves the null-based values undefined (or, on request, runs a label
    shuffle).
    """
    if null_sub is None or len(all_cycle_levels) == 0:
        return None
    sub = null_sub[null_sub["pos_id"].isin(bg_idx_map)]
    n_cycles = len(all_cycle_levels)
    cycle_to_col = {c: i for i, c in enumerate(all_cycle_levels)}
    row_idx = sub["pos_id"].map(bg_idx_map).to_numpy(dtype=np.int64)
    col_idx = sub["cycle"].map(cycle_to_col).to_numpy(dtype=np.int64)
    null_mat = sp.csr_matrix(
        (sub["score"].to_numpy(dtype=np.float32), (row_idx, col_idx)),
        shape=(N, n_cycles)
    )
    return M_mat.dot(null_mat).toarray()  # (n_terms x n_cycles)


def annotate_overlap(overlap, annot, flag_names):
    """Percentage of the distinct genes of the overlap carrying each cross-module flag.

    Returns {pct_<flag>: percent}, or {} when there are no flags.
    """
    if not flag_names:
        return {}
    genes = {p.rsplit(":", 1)[0] for p in overlap}
    if not genes:
        return {f"pct_{f[len('flag_'):]}": np.nan for f in flag_names}
    return {
        f"pct_{f[len('flag_'):]}":
            100.0 * sum(annot.get(g, {}).get(f, False) for g in genes) / len(genes)
        for f in flag_names
    }


def restrict_background_to_coverage(bg_set, coverage_genes):
    """Positions of bg_set whose gene is in coverage_genes."""
    return {p for p in bg_set if p.rsplit(":", 1)[0] in coverage_genes}


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values, in the order of the input."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


# ── Path-sum permulation ──────────────────────────────────────────────────────


def run_permulation_for_terms(terms, descs, obs_scores_dict, background, min_size, max_size,
                             n_perms=10000, seed=1998, annot=None, flag_names=None,
                             perm_chunk_size=1000, caas_null_sub=None, caas_null_cycles=None, allow_label_shuffle=False):
    """Path-sum permulation test of every term of one database in one direction.

    Observed and null term sums are sparse matrix products of the term indicator matrix
    (n_terms x N background positions) with the score vector or the null matrix. Returns
    one dict per term (rows of the characterization table before p_adj, n_scored and sig;
    "_overlap" holds the scored member positions) or [] when no term passes the size filters.

    With the label shuffle, permutations are generated and multiplied in chunks of
    perm_chunk_size, so peak memory scales with the chunk and not with n_perms or N. The
    null mean, sd and exceedance count are accumulated across chunks (sum, sum of squares,
    count), which gives the same values as one unchunked permutation matrix.
    """
    N = len(background)
    if N == 0:
        return []
    bg_idx_map = {p: i for i, p in enumerate(background)}
    bg_set = set(background)

    # 1. Filter valid terms and construct sparse indicator matrix M (n_terms x N)
    valid_terms = []
    rows, cols = [], []
    term_idx = 0

    for term, members in terms.items():
        m_bg = set(members) & bg_set
        mm = len(m_bg)
        if min_size and mm < min_size:
            continue
        if max_size and mm > max_size:
            continue
        
        term_idx_list = [bg_idx_map[p] for p in m_bg]
        if not term_idx_list:
            continue

        rows.extend([term_idx] * len(term_idx_list))
        cols.extend(term_idx_list)
        valid_terms.append((term, descs.get(term, term), m_bg))
        term_idx += 1

    n_terms = len(valid_terms)
    if n_terms == 0:
        return []

    M_mat = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n_terms, N))

    # 2. Score vector V_obs
    V_obs = np.zeros(N, dtype=np.float64)
    for pos_id, score in obs_scores_dict.items():
        if pos_id in bg_idx_map:
            V_obs[bg_idx_map[pos_id]] = float(score)

    obs_sums = M_mat.dot(V_obs)

    # 3. Null source. The CAAS permulation null is preferred over a label shuffle, as
    # fcs_run_permulation does for the FCS path-sum test: a shuffle treats every position
    # score as an independent draw and ignores the phylogenetic dependence the CAAS null
    # accounts for. When the CAAS null is available it replaces the shuffle instead of
    # being combined with it, because an anticonservative test adds no independent
    # evidence. p_value and p_adj reflect whichever null was used.
    null_cycle_sums = caas_null_term_sums(M_mat, bg_idx_map, N, caas_null_sub, caas_null_cycles)

    if null_cycle_sums is None and not allow_label_shuffle:
        # No phylogenetic null: the observed sums are reported and every value that needs a null is undefined. A label
        # shuffle would put a weaker, anticonservative test under the names of the permulation test (it is opt-in).
        out = []
        for i, (term, desc, m_bg) in enumerate(valid_terms):
            driver_positions = {p for p in m_bg if obs_scores_dict.get(p, 0.0) > 0}
            row = dict(
                pathway=term, description=desc, layer_size=len(m_bg), n_pos_with_score=len(driver_positions),
                obs_sum=float(obs_sums[i]), null_mean=np.nan, null_sd=np.nan, perm_nes=np.nan, p_value=np.nan,
                direction="", background_n=N, _overlap=driver_positions,
            )
            row.update(annotate_overlap(driver_positions, annot or {}, flag_names or []))
            out.append(row)
        return out

    if null_cycle_sums is not None:
        n_draws = null_cycle_sums.shape[1]
        null_mu = null_cycle_sums.mean(axis=1)
        null_var = np.maximum(np.square(null_cycle_sums).mean(axis=1) - null_mu**2, 0.0)
        null_sd = np.sqrt(null_var)
        null_sd_safe = np.where(null_sd == 0, 1.0, null_sd)
        perm_nes = (obs_sums - null_mu) / null_sd_safe
        counts = np.sum(null_cycle_sums >= obs_sums[:, None], axis=1)
        pvals = (counts + 1.0) / (n_draws + 1.0)
    else:
        # Label shuffle: the nonzero scores of this direction are placed at random
        # positions of the background (chunked, see the docstring).
        nz_idx = np.where(V_obs > 0)[0]
        nz_vals = V_obs[nz_idx]
        K = len(nz_idx)

        if K == 0:
            # All scores 0 -> no signal
            out = []
            for i, (term, desc, m_bg) in enumerate(valid_terms):
                row = dict(
                    pathway=term, description=desc, layer_size=len(m_bg),
                    n_pos_with_score=0, obs_sum=0.0, null_mean=0.0, null_sd=0.0,
                    perm_nes=0.0, p_value=1.0,
                    direction="depleted", background_n=N,
                    _overlap=set()
                )
                row.update(annotate_overlap(set(), annot or {}, flag_names or []))
                out.append(row)
            return out

        rng = np.random.default_rng(seed)

        # Each chunk builds its own sparse permutation matrix and is folded into the
        # running accumulators, so peak memory is bounded by perm_chunk_size.
        sum_null = np.zeros(n_terms, dtype=np.float64)
        sumsq_null = np.zeros(n_terms, dtype=np.float64)
        counts = np.zeros(n_terms, dtype=np.int64)

        done = 0
        while done < n_perms:
            chunk = min(perm_chunk_size, n_perms - done)

            p_rows = np.empty(K * chunk, dtype=np.int32)
            p_cols = np.empty(K * chunk, dtype=np.int32)
            p_vals = np.tile(nz_vals, chunk)
            for j in range(chunk):
                rnd_idx = rng.choice(N, size=K, replace=False)
                p_rows[j*K : (j+1)*K] = rnd_idx
                p_cols[j*K : (j+1)*K] = j

            P_chunk = sp.csr_matrix((p_vals, (p_rows, p_cols)), shape=(N, chunk))
            null_chunk = M_mat.dot(P_chunk).toarray()  # (n_terms x chunk)

            sum_null += null_chunk.sum(axis=1)
            sumsq_null += np.square(null_chunk).sum(axis=1)
            counts += np.sum(null_chunk >= obs_sums[:, None], axis=1)

            done += chunk

        null_mu = sum_null / n_perms
        # Population variance from sum-of-squares (matches np.std's default ddof=0);
        # clip at 0 to guard float round-off pushing a near-zero variance negative.
        null_var = np.maximum(sumsq_null / n_perms - null_mu**2, 0.0)
        null_sd = np.sqrt(null_var)
        null_sd_safe = np.where(null_sd == 0, 1.0, null_sd)

        perm_nes = (obs_sums - null_mu) / null_sd_safe
        pvals = (counts + 1.0) / (n_perms + 1.0)

    # 4. Format results
    out = []
    for i, (term, desc, m_bg) in enumerate(valid_terms):
        driver_positions = {p for p in m_bg if obs_scores_dict.get(p, 0.0) > 0}
        row = dict(
            pathway=term,
            description=desc,
            layer_size=len(m_bg),
            n_pos_with_score=len(driver_positions),
            obs_sum=float(obs_sums[i]),
            null_mean=float(null_mu[i]),
            null_sd=float(null_sd[i]),
            perm_nes=float(perm_nes[i]),
            p_value=float(pvals[i]),
            direction=("enriched" if perm_nes[i] > 0 else "depleted"),
            background_n=N,
            _overlap=driver_positions
        )
        row.update(annotate_overlap(driver_positions, annot or {}, flag_names or []))
        out.append(row)

    return out


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    universe_genes = load_universe_genes(args.universe)
    print(f"[posenrich] universe (cleaned_background) genes: {len(universe_genes)}", flush=True)

    background = build_background(args.background, universe_genes)
    obs = pd.read_csv(args.obs_scores, sep="\t")
    for col in ("Gene", "Position", "CAAS_score", "side"):
        if col not in obs.columns:
            sys.exit(f"[posenrich] obs-scores missing column: {col}")
    obs = obs[obs["Gene"].isin(universe_genes)].copy()
    obs["CAAS_score"] = pd.to_numeric(obs["CAAS_score"], errors="coerce").fillna(0.0)
    obs["pos_id"] = obs["Gene"].astype(str) + ":" + obs["Position"].astype(str)

    background = sorted(set(background) | set(obs["pos_id"]))
    bg_set = set(background)
    N = len(background)
    print(f"[posenrich] background positions: {N} | scored: {(obs['CAAS_score']>0).sum()}", flush=True)

    coverage_restricted_bg = {}
    for db_name, coverage_path in (("cosmic_orthogroups", args.cosmic_coverage),
                                    ("pai3d_orthogroups", args.pai3d_coverage)):
        if coverage_path and os.path.exists(coverage_path):
            coverage_genes = load_universe_genes(coverage_path)
            restricted_bg = restrict_background_to_coverage(bg_set, coverage_genes)
            coverage_restricted_bg[db_name] = sorted(restricted_bg)
            print(f"[posenrich] {db_name}: background restricted to "
                  f"{len(coverage_genes)} coverage genes "
                  f"({len(restricted_bg)}/{N} positions)", flush=True)

    gmts = load_gmts(args.gmt_dir)
    char_layers = read_charset(args.characterization)
    annot, flag_names = load_annot(args.annot_file)
    if flag_names:
        print(f"[posenrich] cross-module flags: {', '.join(flag_names)} "
              f"({len(annot)} genes annotated)", flush=True)

    sources = {db: (terms, descs, True) for db, (terms, descs) in gmts.items()}
    if char_layers:
        char_terms = {name: members for name, (desc, members) in char_layers.items()}
        char_descs = {name: desc for name, (desc, members) in char_layers.items()}
        sources["characterization"] = (char_terms, char_descs, False)
    print(f"[posenrich] sources: {len(sources)} ({', '.join(sources)})", flush=True)

    # Hypothesis descriptors are read per direction from the row direction_rows()
    # keeps, so a leading-edge entry describes the same side that set its score.
    # position_scores.tsv carries the hypothesis list as participating_hypotheses;
    # it is written out here under the leading-edge schema's supporting_hypotheses.
    has_hyp = "n_hypotheses" in obs.columns
    supp_col = next((c for c in ("participating_hypotheses", "supporting_hypotheses") if c in obs.columns), None)

    caas_null_by_direction = None
    if args.caas_null_prepped:
        caas_null_by_direction, caas_null_cycles = load_prepped_caas_null(args.caas_null_prepped)
        caas_null_long = None
        if caas_null_by_direction is not None:
            print(f"[posenrich] CAAS permulation null: {len(caas_null_cycles)} cycles "
                  f"loaded (prepped) from {args.caas_null_prepped} -> used as the primary null", flush=True)
        else:
            print("[posenrich] no CAAS permulation null supplied", flush=True)
    else:
        caas_null_long, caas_null_cycles = load_caas_cycle_null(args.caas_cycle_null)
        if caas_null_long is not None:
            print(f"[posenrich] CAAS permulation null: {len(caas_null_cycles)} cycles "
                  f"loaded from {args.caas_cycle_null} -> used as the primary null", flush=True)
        else:
            print("[posenrich] no CAAS permulation null supplied", flush=True)
    if caas_null_by_direction is None and caas_null_long is None:
        print("[posenrich] without a null, p_value, p_adj, perm_nes and the null columns are NA and nothing is significant"
              + ("" if args.allow_label_shuffle else " (--allow-label-shuffle runs the weaker label-shuffle test instead)"), flush=True)

    directions = ["global", "top", "bottom"]
    rows = []
    leading_edge_rows = []
    for direction in directions:
        dir_rows = direction_rows(obs, direction)
        obs_scores = dict(zip(dir_rows["pos_id"], dir_rows["CAAS_score"]))
        hyp_dict = dict(zip(dir_rows["pos_id"], dir_rows["n_hypotheses"])) if has_hyp else {}
        supp_dict = (dict(zip(dir_rows["pos_id"], dir_rows[supp_col].fillna("").astype(str)))
                     if supp_col else {})
        n_scored = sum(1 for s in obs_scores.values() if s > 0)
        if n_scored == 0:
            continue
        print(f"[posenrich] {direction}: {n_scored} scored positions | running Path Sum Permulation (N_perms={args.n_perms})...", flush=True)

        if caas_null_by_direction is not None:
            null_sub = caas_null_by_direction.get(direction)
        elif caas_null_long is not None:
            null_sub = null_direction_subset(caas_null_long, direction)
        else:
            null_sub = None

        for db, (terms, descs, apply_size_filter) in sources.items():
            db_bg = coverage_restricted_bg.get(db, background)
            res = run_permulation_for_terms(
                terms, descs, obs_scores, db_bg,
                args.min_size if apply_size_filter else 0,
                args.max_size if apply_size_filter else 0,
                n_perms=args.n_perms, seed=args.seed,
                annot=annot, flag_names=flag_names,
                perm_chunk_size=args.perm_chunk_size,
                caas_null_sub=null_sub, caas_null_cycles=caas_null_cycles,
                allow_label_shuffle=args.allow_label_shuffle
            )
            if not res:
                continue

            no_null = null_sub is None and not args.allow_label_shuffle
            padj = np.full(len(res), np.nan) if no_null else bh_adjust([r["p_value"] for r in res])
            for i, r in enumerate(res):
                r["p_adj"] = padj[i]
                r["n_scored"] = n_scored
                r["sig"] = bool(r["p_adj"] < args.padj_thr and r["perm_nes"] > 0)
                overlap = r.pop("_overlap")
                if r["sig"]:
                    for pos_id in sorted(overlap):
                        leading_edge_rows.append(dict(
                            ranking=direction, database=db, pathway=r["pathway"],
                            gene=pos_id.rsplit(":", 1)[0],
                            gene_position=pos_id,
                            CAAS_score=obs_scores.get(pos_id, 0.0),
                            n_hypotheses=hyp_dict.get(pos_id, 1),
                            supporting_hypotheses=supp_dict.get(pos_id, "")
                        ))
                rows.append(dict(ranking=direction, database=db, **r))

    result_cols = (
        ["ranking", "database", "pathway", "description", "layer_size", "n_pos_with_score",
         "obs_sum", "null_mean", "null_sd", "perm_nes", "p_value", "p_adj",
         "direction", "background_n"]
        + [f"pct_{f[len('flag_'):]}" for f in flag_names]
        + ["n_scored", "sig"]
    )
    results = pd.DataFrame(rows, columns=result_cols) if not rows else pd.DataFrame(rows, columns=result_cols)
    if not results.empty:
        results = results.sort_values(["p_adj", "p_value"], na_position="last")
    out_path = os.path.join(args.output_dir, "posenrich_characterization.tsv")
    results.to_csv(out_path, sep="\t", index=False)
    print(f"[posenrich] wrote {out_path} ({len(results)} rows)", flush=True)

    leading_edge_cols = ["ranking", "database", "pathway", "gene", "gene_position", "CAAS_score"]
    if has_hyp:
        leading_edge_cols.extend(["n_hypotheses", "supporting_hypotheses"])
    leading_edge = pd.DataFrame(
        leading_edge_rows,
        columns=leading_edge_cols,
    )
    le_path = os.path.join(args.output_dir, "posenrich_leading_edge.tsv")
    leading_edge.to_csv(le_path, sep="\t", index=False)
    print(f"[posenrich] wrote {le_path} ({len(leading_edge)} rows)", flush=True)


if __name__ == "__main__":
    main()
