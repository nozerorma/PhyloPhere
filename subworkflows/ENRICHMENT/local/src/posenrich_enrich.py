#!/usr/bin/env python3
# =============================================================================
# posenrich_enrich.py - Position-wise enrichment via Path Sum Permulation Tests
# =============================================================================
# Rank-free, magnitude-aware position-level gene-set enrichment via sparse-matrix
# vectorized label permutations over the full ~1.47M position background.
#
# Rationale:
#   At the position level (honest testable background of ~1.47M alignment columns),
#   <0.2% of positions carry a non-zero CAAS score. Fixed-cutoff Fisher tests
#   suffer from extreme sample-size sensitivity (N=1.47M), assigning p < 0.001 to
#   biologically meaningless deviations (e.g. fold = 1.01) while discarding continuous
#   score magnitude.
#
#   Position-Level Path Sum Permulation solves this by summing raw CAAS score
#   magnitudes for every position in a gene set:
#     T_obs = sum_{p in term} s_p
#   The ~99.8% zero-scoring positions contribute 0 to the pathway sum.
#   Permuting the non-zero scores randomly across the ~1.47M background pool
#   N_perms times yields an exact, magnitude-weighted empirical null distribution.
#
# Output:
#   posenrich_characterization.tsv : ranking, database, pathway, description,
#                                    layer_size, n_pos_with_score, obs_sum,
#                                    null_mean, null_sd, perm_nes, p_value,
#                                    p_adj, direction, background_n, n_scored, sig
#   posenrich_leading_edge.tsv     : gene:position driver members for significant terms
# =============================================================================

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp


# ── args ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Position-level Path Sum Permulation enrichment.")
    p.add_argument("--obs-scores", required=True,
                   help="position_scores.tsv (Gene, Position, CAAS_score, change_side)")
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
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-size", type=int, default=5,
                   help="min positions per set in background (GMT sources only)")
    p.add_argument("--max-size", type=int, default=0,
                   help="max positions per set in background (0 = no cap; GMT sources only)")
    p.add_argument("--n-perms", type=int, default=100000,
                   help="number of label permutations for Path Sum Permulation (default 10000)")
    p.add_argument("--seed", type=int, default=42,
                   help="random seed for permulations (default 42)")
    p.add_argument("--padj-thr", type=float, default=0.15,
                   help="BH-adjusted p-value significance threshold")
    p.add_argument("--position-lists-dir", required=False, default=None,
                   help="Accepted for backward compatibility (unused in continuous permulation)")
    p.add_argument("--char-fracs", default=None, help="Accepted for backward compatibility (unused)")
    p.add_argument("--fold-thr", type=float, default=1.5, help="Accepted for backward compatibility (unused)")
    return p.parse_args()


# ── universe / background from MAP files ─────────────────────────────────────
def load_universe_genes(path):
    genes = set()
    with open(path) as f:
        for line in f:
            g = line.strip()
            if g and g != "Gene":
                genes.add(g)
    return genes


def build_background(background_file, universe_genes):
    """Honest testable position background = positions CAAStools tested,
    RESTRICTED to cleaned-background genes."""
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


# ── GMT / characterization-layer loading ─────────────────────────────────────
def load_gmts(gmt_dir):
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
    """Cross-module corroboration flags, per gene (SCORING's fcs_stats.tsv)."""
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


# ── observed scores ──────────────────────────────────────────────────────────
def direction_scores(df, direction):
    """Return {pos_id: CAAS_score} for a direction (global/top/bottom)."""
    if direction == "global":
        sub = df
    elif direction == "top":
        sub = df[df["change_side"].isin(["top", "both"])]
    else:
        sub = df[df["change_side"].isin(["bottom", "both"])]
    return dict(zip(sub["pos_id"], sub["CAAS_score"]))


def annotate_overlap(overlap, annot, flag_names):
    """Cross-module corroboration for distinct genes in overlap."""
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
    """Restrict background positions to genes in coverage_genes."""
    return {p for p in bg_set if p.rsplit(":", 1)[0] in coverage_genes}


def bh_adjust(pvals):
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


# ── Sparse Path Sum Permulation Engine ───────────────────────────────────────
def run_permulation_for_terms(terms, descs, obs_scores_dict, background, min_size, max_size,
                             n_perms=10000, seed=42, annot=None, flag_names=None):
    """
    Position-Level Path Sum Permulation test.
    Vectorized sparse matrix multiplication over background pool (N positions).
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

    # 3. Sparse permutation matrix P (N x n_perms)
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
                perm_nes=0.0, p_value=1.0, direction="depleted", background_n=N,
                _overlap=set()
            )
            row.update(annotate_overlap(set(), annot or {}, flag_names or []))
            out.append(row)
        return out

    rng = np.random.default_rng(seed)
    p_rows = np.empty(K * n_perms, dtype=np.int32)
    p_cols = np.empty(K * n_perms, dtype=np.int32)
    p_vals = np.tile(nz_vals, n_perms)

    for j in range(n_perms):
        rnd_idx = rng.choice(N, size=K, replace=False)
        p_rows[j*K : (j+1)*K] = rnd_idx
        p_cols[j*K : (j+1)*K] = j

    P_mat = sp.csr_matrix((p_vals, (p_rows, p_cols)), shape=(N, n_perms))

    # 4. Multiply sparse matrices to get null sums
    null_sums = M_mat.dot(P_mat).toarray()  # (n_terms x n_perms)

    null_mu = np.mean(null_sums, axis=1)
    null_sd = np.std(null_sums, axis=1)
    null_sd_safe = np.where(null_sd == 0, 1.0, null_sd)

    perm_nes = (obs_sums - null_mu) / null_sd_safe
    counts = np.sum(null_sums >= obs_sums[:, None], axis=1)
    pvals = (counts + 1.0) / (n_perms + 1.0)

    # 5. Format results
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


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    universe_genes = load_universe_genes(args.universe)
    print(f"[posenrich] universe (cleaned_background) genes: {len(universe_genes)}", flush=True)

    background = build_background(args.background, universe_genes)
    obs = pd.read_csv(args.obs_scores, sep="\t")
    for col in ("Gene", "Position", "CAAS_score", "change_side"):
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

    hyp_dict = dict(zip(obs["pos_id"], obs["n_hypotheses"])) if "n_hypotheses" in obs.columns else {}
    supp_dict = dict(zip(obs["pos_id"], obs["supporting_hypotheses"])) if "supporting_hypotheses" in obs.columns else {}

    directions = ["global", "top", "bottom"]
    rows = []
    leading_edge_rows = []
    for direction in directions:
        obs_scores = direction_scores(obs, direction)
        n_scored = sum(1 for s in obs_scores.values() if s > 0)
        if n_scored == 0:
            continue
        print(f"[posenrich] {direction}: {n_scored} scored positions | running Path Sum Permulation (N_perms={args.n_perms})...", flush=True)

        for db, (terms, descs, apply_size_filter) in sources.items():
            db_bg = coverage_restricted_bg.get(db, background)
            res = run_permulation_for_terms(
                terms, descs, obs_scores, db_bg,
                args.min_size if apply_size_filter else 0,
                args.max_size if apply_size_filter else 0,
                n_perms=args.n_perms, seed=args.seed,
                annot=annot, flag_names=flag_names
            )
            if not res:
                continue

            padj = bh_adjust([r["p_value"] for r in res])
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
    if hyp_dict:
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
