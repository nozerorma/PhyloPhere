#!/usr/bin/env python3
# =============================================================================
# posenrich_enrich.py — Position-wise enrichment via the XL-mHG test
# =============================================================================
# This is NOT the gene-level FCS: FCS is a Wilcoxon rank-shift over gene scores;
# this is a position-wise over-representation test that accounts for per-position
# CAAS evidence via the minimum-hypergeometric statistic. Rationale: the honest
# testable position background (all non-gapped "selected" columns of the
# cleaned-background genes) is <1% non-zero in CAAS_score, which pins the
# Mann-Whitney AUC at ~0.5. The XL-mHG (minimum-hypergeometric) test is a
# top-concentration test: it finds the score cutoff where a position set is most
# over-represented at the head of the ranking, so it is robust to the huge tied
# zero-mass and remains powerful. See the gene-level FCS (fcs_enrich.R) which
# keeps Wilcoxon — that background is ~40%+ non-zero, so it has resolution there.
#
# Observed statistic  : XL-mHG on positions ranked by CAAS_score (from the
#                       scoring module's position_scores.tsv), per direction
#                       (global / top / bottom via change_side).
# Analytic p (p.adj)  : the exact XL-mHG p-value, BH-adjusted per (direction, db).
#                       X>=2 neutralises one-hit-wonders; L confines the peak to
#                       the top of the ranking. This is the ONLY significance here:
#                       the exact XL-mHG p already corrects the cutoff-optimisation,
#                       so no permulation (p.perm) is used. (The CAAS permulation
#                       null is also a NAIVE random relabeling — its pairs are
#                       random cross-tree, not the observed matched sister pairs —
#                       so it would not isolate the phenotype at position level.
#                       n_le_genes flags within-gene clumping instead.)
# =============================================================================

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from xlmhglite import xlmhg_test


# ── args ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Position-level XL-mHG enrichment.")
    p.add_argument("--obs-scores", required=True,
                   help="position_scores.tsv (Gene, Position, CAAS_score, change_side)")
    p.add_argument("--gmt-dir", required=True)
    p.add_argument("--characterization", default=None,
                   help="characterization_layers.tsv (global overlap layers)")
    p.add_argument("--universe", required=True,
                   help="cleaned_background_main.txt gene list (postproc-surviving)")
    p.add_argument("--background", required=True,
                   help="caastools background.output (gene<TAB>tested positions); "
                        "restricted to --universe genes = the position background")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-size", type=int, default=5,
                   help="min positions per set (in background)")
    p.add_argument("--max-size", type=int, default=0,
                   help="max positions per set (0 = no cap)")
    p.add_argument("--min-overlap", type=int, default=2,
                   help="XL-mHG X: min set positions above the cutoff (>=2 kills one-hit-wonders)")
    p.add_argument("--l-frac", type=float, default=1.0,
                   help="XL-mHG L as a fraction of the scored region (peak confined to top)")
    p.add_argument("--char-fracs", default="0.10,0.05,0.01",
                   help="characterization foreground top-fractions among scored positions")
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
    """Honest testable position background = the positions CAAStools actually
    tested (caastools `background.output`: `gene<TAB>comma-separated positions`,
    or NULL), RESTRICTED to the cleaned-background genes (i.e. background.output
    minus the negative overlap with cleaned_background — the postproc-dropped
    genes). Positions are in the same prot_ali_col space as the GMTs and scores."""
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


# ── GMT loading ──────────────────────────────────────────────────────────────
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


# ── XL-mHG on one ranking ────────────────────────────────────────────────────
def rank_index(scores_by_pos, background):
    """Order background positions by score desc (stable tie-break by pos id).
    Returns (ordered_positions, index_of_pos, n_scored)."""
    scores = np.array([scores_by_pos.get(p, 0.0) for p in background], dtype=float)
    # stable, deterministic: sort by (-score, position string)
    order = sorted(range(len(background)), key=lambda i: (-scores[i], background[i]))
    ordered = [background[i] for i in order]
    idx = {p: r for r, p in enumerate(ordered)}
    n_scored = int(np.count_nonzero(scores > 0))
    return ordered, idx, n_scored


def percentile_thresholds(scores_by_pos, fracs=(0.10, 0.05, 0.01)):
    """Score cutoffs for the top-10/5/1% of the SCORED positions in a direction
    (the scoring module's is_top_Xpct convention)."""
    vals = sorted((s for s in scores_by_pos.values() if s > 0), reverse=True)
    thr = {}
    for f in fracs:
        if vals:
            k = max(1, int(round(f * len(vals))))
            thr[f] = vals[k - 1]
        else:
            thr[f] = float("inf")
    return thr


def run_xlmhg_for_sets(terms, idx, N, n_scored, x, l_frac, num_g, num_g_max,
                       scores_by_pos, ordered, tiers):
    """Run XL-mHG for every term in one database against one ranking."""
    L = max(x, min(N, int(round(l_frac * max(n_scored, 1)))))
    v = np.zeros(N, dtype=np.uint8)
    out = []
    for term, members in terms.items():
        member_ranks = [idx[m] for m in members if m in idx]
        K = len(member_ranks)
        if K < num_g or (num_g_max > 0 and K > num_g_max):
            continue
        n_scored_members = sum(1 for m in members
                               if scores_by_pos.get(m, 0.0) > 0 and m in idx)
        # composition: set members in the top 10/5/1% of the score distribution
        n_top10 = sum(1 for m in members if scores_by_pos.get(m, 0.0) >= tiers[0.10])
        n_top5 = sum(1 for m in members if scores_by_pos.get(m, 0.0) >= tiers[0.05])
        n_top1 = sum(1 for m in members if scores_by_pos.get(m, 0.0) >= tiers[0.01])
        v[member_ranks] = 1
        try:
            stat, cutoff, pval = xlmhg_test(v, X=x, L=L)
        except Exception:
            stat, cutoff, pval = 1.0, 0, 1.0
        finally:
            v[member_ranks] = 0
        # leading edge = members at or above the mHG cutoff
        k_cut = sum(1 for r in member_ranks if r < cutoff)
        fold = ((k_cut / cutoff) / (K / N)) if cutoff > 0 and K > 0 else np.nan
        le = [ordered[r] for r in sorted(member_ranks) if r < cutoff]
        # distinct genes in the leading edge — a cheap, analytic guard against
        # within-gene clumping (an enrichment driven by 1-2 genes is the
        # non-independence artifact the permulation null would have flagged).
        n_le_genes = len({m.rsplit(":", 1)[0] for m in le})
        out.append(dict(pathway=term, stat=float(stat), cutoff=int(cutoff),
                        pval=float(pval), set_size=K,
                        n_scored_members=n_scored_members,
                        n_top10=n_top10, n_top5=n_top5, n_top1=n_top1,
                        leading_edge_size=k_cut, n_le_genes=n_le_genes,
                        fold=fold, leading_edge=",".join(le)))
    return out


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


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    char_fracs = [float(x) for x in args.char_fracs.split(",")]

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
    # Union scored positions defensively (a scored position is by definition
    # testable, so never drop it) — a no-op when background.output is complete.
    n_bg_only = len(set(background))
    background = sorted(set(background) | set(obs["pos_id"]))
    N = len(background)
    print(f"[posenrich] background positions: {N} (background.output∩universe={n_bg_only}) "
          f"| scored: {(obs['CAAS_score']>0).sum()}", flush=True)

    gmts = load_gmts(args.gmt_dir)
    print(f"[posenrich] GMT databases: {len(gmts)}", flush=True)

    directions = ["global", "top", "bottom"]
    rows = []
    for direction in directions:
        obs_scores = direction_scores(obs, direction)
        ordered, idx, n_scored = rank_index(obs_scores, background)
        if n_scored == 0:
            continue
        tiers = percentile_thresholds(obs_scores)
        print(f"[posenrich] {direction}: {n_scored} scored positions ranked", flush=True)

        for db, (terms, descs) in gmts.items():
            res = run_xlmhg_for_sets(terms, idx, N, n_scored, args.min_overlap,
                                     args.l_frac, args.min_size, args.max_size,
                                     obs_scores, ordered, tiers)
            if not res:
                continue
            padj = bh_adjust([r["pval"] for r in res])
            for i, r in enumerate(res):
                rows.append(dict(ranking=direction, database=db,
                                 pathway=r["pathway"],
                                 description=descs.get(r["pathway"], r["pathway"]),
                                 stat=r["stat"], cutoff=r["cutoff"],
                                 pval=r["pval"], p_adj=padj[i],
                                 set_size=r["set_size"],
                                 n_scored_members=r["n_scored_members"],
                                 n_top10=r["n_top10"], n_top5=r["n_top5"],
                                 n_top1=r["n_top1"],
                                 leading_edge_size=r["leading_edge_size"],
                                 n_le_genes=r["n_le_genes"],
                                 fold=r["fold"], leading_edge=r["leading_edge"]))

    results = pd.DataFrame(rows)
    if not results.empty:
        results = results.sort_values(["p_adj", "pval"], na_position="last")
    out_path = os.path.join(args.output_dir, "posenrich_results.tsv")
    results.to_csv(out_path, sep="\t", index=False)
    print(f"[posenrich] wrote {out_path} ({len(results)} rows)", flush=True)

    # ── Characterization (Fisher overlap at top 10/5/1% of scored) ────────────
    layers = read_charset(args.characterization)
    char_rows = []
    if layers:
        gscore = direction_scores(obs, "global")
        scored_sorted = sorted([p for p in background if gscore.get(p, 0.0) > 0],
                               key=lambda p: -gscore[p])
        bg_set = set(background)
        for frac in char_fracs:
            n = max(1, int(round(frac * len(scored_sorted))))
            foreground = set(scored_sorted[:n])
            for name, (desc, members) in layers.items():
                M = members & bg_set
                k = len(foreground & M)
                mm, nn, NN = len(M), n, N
                expected = nn * mm / NN if NN else np.nan
                fold = (k / expected) if expected and expected > 0 else np.nan
                table = [[k, nn - k], [mm - k, NN - nn - (mm - k)]]
                try:
                    _, p = fisher_exact(table, alternative="two-sided")
                except Exception:
                    p = np.nan
                char_rows.append(dict(top_frac=frac, layer=name, description=desc,
                                      layer_size=mm, foreground_size=nn, observed=k,
                                      expected=expected, fold_enrichment=fold,
                                      p_value=p,
                                      direction=("enriched" if fold and fold >= 1 else "depleted")))
    char_df = pd.DataFrame(char_rows)
    if not char_df.empty:
        char_df["p_adj"] = bh_adjust(char_df["p_value"].fillna(1.0).values)
    char_path = os.path.join(args.output_dir, "posenrich_characterization.tsv")
    char_df.to_csv(char_path, sep="\t", index=False)
    print(f"[posenrich] wrote {char_path} ({len(char_df)} rows)", flush=True)


if __name__ == "__main__":
    main()
