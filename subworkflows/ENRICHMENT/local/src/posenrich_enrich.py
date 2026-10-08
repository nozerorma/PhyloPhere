#!/usr/bin/env python3
# posenrich_enrich.py — Position-level gene-set enrichment by the Lachenbruch two-part test.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
PosenrichEnrich: tests whether the CAAS scores of the alignment positions of a gene set stand out among the scored
positions of the background and exceed what the CAAS permulation null produces, for every gene set of the GMT files and
of the characterization layers, in three directions (global, top, bottom).

Test. The Lachenbruch two-part test, the same as the FCS (fcs_enrich.R, fcs_lach_prepare / fcs_lach_stats), by the same
formulas. Part 1 (prevalence) is the upper tail of the hypergeometric distribution of the scored positions of the term
among the m scored positions of the direction over the background. Part 2 (magnitude) is the rank-sum of the CAAS scores
of the scored positions of the term against those of the other scored positions, by the normal approximation with
continuity correction and the tie term of the scored positions; it needs at least 2 scored positions in the term and 2
outside it. Each p is floored at 1e-15, turned into chi-square with 1 df, and the sum is chi-square with 2 df (lach_pval).
A fixed-cutoff test (Fisher) over a background of this size would call biologically trivial deviations significant and
would discard the magnitude.

Null. The null is the real permulation cycles of the CAAS null (perm_pos_cycle_caas.tsv.gz). The observed scores and every
cycle go through the same function, so lach_p.perm = (1 + number of cycles with chi-square >= observed) / (n_cycles + 1)
compares like with like. Without a CAAS null lach_p.perm is NA, nothing is significant, and the analytic columns are still
written.

lach_p.adj is the BH adjustment of lach_pval within each (ranking, database). A term is significant (sig) when
lach_p.adj < --fdr-lachenbruch and lach_p.perm < --pperm-thr, the gates of the Lachenbruch vote of the FCS.

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
                n_pos_with_score, lach_chi_binary, lach_chi_nonzero, lach_chi_total, lach_frac_magnitude,
                lach_pval, lach_p.adj, lach_p.perm, background_n,
                pct_<flag> (one per flag_* column of --annot-file), n_scored, sig
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
from scipy import stats as st


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Position-level Lachenbruch two-part enrichment.")
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
                        "supplied, its real cycles are the null of the permulation p (lach_p.perm). "
                        "Omit or pass a NO_FILE* sentinel for no null (the permulation p is then NA). "
                        "Ignored when --caas-null-prepped is given.")
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
    p.add_argument("--fdr-lachenbruch", type=float, default=0.15,
                   help="BH FDR gate of the Lachenbruch test (analytic p), as the FCS")
    p.add_argument("--pperm-thr", type=float, default=0.025,
                   help="permulation p gate of the Lachenbruch test (lach_p.perm), as the FCS")
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


def caas_null_matrix(bg_idx_map, N, null_sub, all_cycle_levels):
    """The CAAS permulation null as a sparse (N background positions x n_cycles) matrix of scores, or None without a null.

    The columns span all_cycle_levels, every cycle of the file, so a cycle without hits is a null draw of zeros.
    """
    if null_sub is None or len(all_cycle_levels) == 0:
        return None
    sub = null_sub[null_sub["pos_id"].isin(bg_idx_map)]
    cycle_to_col = {c: i for i, c in enumerate(all_cycle_levels)}
    row_idx = sub["pos_id"].map(bg_idx_map).to_numpy(dtype=np.int64)
    col_idx = sub["cycle"].map(cycle_to_col).to_numpy(dtype=np.int64)
    return sp.csr_matrix((sub["score"].to_numpy(dtype=np.float32), (row_idx, col_idx)),
                         shape=(N, len(all_cycle_levels)))


# ── Lachenbruch two-part test (the definition of fcs_enrich.R: fcs_lach_prepare, fcs_lach_stats) ──


def rank_sum_sd(n1, n2, tie):
    """Standard deviation of the rank sum of n1 values against n2, with the tie term (sum of t^3 - t over the tied groups)."""
    n = n1 + n2
    return np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / np.maximum(n * (n - 1), 1)))


def lachenbruch_columns(X, M, n1, n_genes):
    """Lachenbruch chi-squares of the terms M (n_terms x n_genes, 0/1) in every column of X (n_genes x n_cols, sparse, >= 0).

    A value of 0 is no signal. Part 1 is the upper tail of the hypergeometric distribution of the scored positions of the
    term (n1 positions) among the m scored positions of the column over n_genes; part 2 is the rank-sum of the scored
    positions of the term against the other scored positions of the column, normal approximation with continuity correction
    and the tie term of the scored positions, defined when the term has at least 2 scored positions and so do the others.
    Each p is floored at 1e-15 and turned into chi-square with 1 df. Returns (chi1, chi2) as (n_terms x n_cols) arrays.
    """
    X = sp.csc_matrix(X, dtype=np.float64)
    X.data[~(X.data > 0)] = 0.0
    X.eliminate_zeros()
    ranks = X.copy()
    tie = np.zeros(X.shape[1])
    m = np.diff(X.indptr).astype(np.float64)
    for j in range(X.shape[1]):
        a, b = X.indptr[j], X.indptr[j + 1]
        if b > a:
            v = X.data[a:b]
            ranks.data[a:b] = st.rankdata(v)
            _, c = np.unique(v, return_counts=True)
            tie[j] = np.sum(c.astype(np.float64) ** 3 - c)
    pos = X.copy()
    pos.data[:] = 1.0
    n1 = np.asarray(n1, dtype=np.float64)
    k1 = M.dot(pos).toarray()
    mm = m[None, :]
    p1 = st.hypergeom.sf(k1 - 1, n_genes, mm, n1[:, None])
    chi1 = st.chi2.isf(np.maximum(p1, 1e-15), 1)
    n2p = mm - k1
    u2 = M.dot(ranks).toarray() - k1 * (k1 + 1) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (u2 - k1 * n2p / 2 - 0.5) / rank_sum_sd(k1, n2p, tie[None, :])
    chi2 = st.chi2.isf(np.maximum(st.norm.sf(z), 1e-15), 1)
    chi2[~((k1 >= 2) & (n2p >= 2)) | np.isnan(chi2)] = 0.0
    return chi1, chi2


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


# ── Terms ─────────────────────────────────────────────────────────────────────


def run_enrichment_for_terms(terms, descs, obs_scores_dict, background, min_size, max_size,
                             annot=None, flag_names=None, caas_null_sub=None, caas_null_cycles=None):
    """Lachenbruch test of every term of one database in one direction.

    The term indicator matrix (n_terms x N background positions) is multiplied with the observed score vector and with
    every cycle of the CAAS null (lachenbruch_columns). Returns one dict per term (the rows of the characterization table
    before lach_p.adj, n_scored and sig; "_overlap" holds the scored member positions) or [] when no term passes the size
    filters. Without a null, lach_p.perm is NaN and the analytic columns are still defined.
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

    # 3. The test. Column 0 is the observed ranking and the other columns are the cycles of the CAAS null.
    null_mat = caas_null_matrix(bg_idx_map, N, caas_null_sub, caas_null_cycles)
    X = sp.csc_matrix(V_obs.reshape(-1, 1)) if null_mat is None else sp.hstack([sp.csc_matrix(V_obs.reshape(-1, 1)), null_mat], format="csc")
    sizes = np.asarray(M_mat.sum(axis=1)).ravel()
    chi1, chi2 = lachenbruch_columns(X, M_mat, sizes, N)
    chi = chi1 + chi2
    pval = st.chi2.sf(chi[:, 0], 2)
    if null_mat is None:
        perm = np.full(n_terms, np.nan)
    else:
        perm = (np.sum(chi[:, 1:] >= chi[:, [0]] - 1e-12, axis=1) + 1.0) / (chi.shape[1] - 1 + 1.0)

    # 4. Format results
    out = []
    for i, (term, desc, m_bg) in enumerate(valid_terms):
        driver_positions = {p for p in m_bg if obs_scores_dict.get(p, 0.0) > 0}
        row = dict(
            pathway=term,
            description=desc,
            layer_size=len(m_bg),
            n_pos_with_score=len(driver_positions),
            lach_chi_binary=float(chi1[i, 0]),
            lach_chi_nonzero=float(chi2[i, 0]),
            lach_chi_total=float(chi[i, 0]),
            lach_frac_magnitude=float(chi2[i, 0] / chi[i, 0]) if chi[i, 0] > 0 else np.nan,
            lach_pval=float(pval[i]),
            background_n=N,
            _overlap=driver_positions
        )
        row["lach_p.perm"] = float(perm[i])
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
        print("[posenrich] without a null, lach_p.perm is NA and nothing is significant", flush=True)

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
        print(f"[posenrich] {direction}: {n_scored} scored positions | running the Lachenbruch two-part test...", flush=True)

        if caas_null_by_direction is not None:
            null_sub = caas_null_by_direction.get(direction)
        elif caas_null_long is not None:
            null_sub = null_direction_subset(caas_null_long, direction)
        else:
            null_sub = None

        for db, (terms, descs, apply_size_filter) in sources.items():
            db_bg = coverage_restricted_bg.get(db, background)
            res = run_enrichment_for_terms(
                terms, descs, obs_scores, db_bg,
                args.min_size if apply_size_filter else 0,
                args.max_size if apply_size_filter else 0,
                annot=annot, flag_names=flag_names,
                caas_null_sub=null_sub, caas_null_cycles=caas_null_cycles
            )
            if not res:
                continue

            lach_adj = bh_adjust([r["lach_pval"] for r in res])
            for i, r in enumerate(res):
                r["lach_p.adj"] = lach_adj[i]
                r["n_scored"] = n_scored
                r["sig"] = bool(lach_adj[i] < args.fdr_lachenbruch and r["lach_p.perm"] < args.pperm_thr)
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
         "lach_chi_binary", "lach_chi_nonzero", "lach_chi_total", "lach_frac_magnitude",
         "lach_pval", "lach_p.adj", "lach_p.perm", "background_n"]
        + [f"pct_{f[len('flag_'):]}" for f in flag_names]
        + ["n_scored", "sig"]
    )
    results = pd.DataFrame(rows, columns=result_cols) if not rows else pd.DataFrame(rows, columns=result_cols)
    if not results.empty:
        results = results.sort_values(["lach_p.adj", "lach_pval"], na_position="last")
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
