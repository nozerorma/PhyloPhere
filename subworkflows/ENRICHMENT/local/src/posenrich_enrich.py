#!/usr/bin/env python3
# =============================================================================
# posenrich_enrich.py — Position-wise enrichment via fixed-cutoff Fisher tests
# =============================================================================
# This is NOT the gene-level FCS: FCS is a Wilcoxon rank-shift over gene scores;
# this is a position-wise over-representation test. Rationale for testing at the
# position level at all: the honest testable position background (all
# non-gapped "selected" columns of the cleaned-background genes) is <1%
# non-zero in CAAS_score, which pins a rank-shift test's AUC at ~0.5. See the
# gene-level FCS (fcs_enrich.R) which keeps Wilcoxon — that background is ~40%+
# non-zero, so it has resolution there.
#
# Observed statistic  : Fisher's exact test, one 2x2 table per
#                       (ranking, database, pathway, top_frac): foreground =
#                       the top 10/5/1% of scored positions in that ranking;
#                       background = the full honest tested position pool
#                       (never restricted to "scored", so negative evidence is
#                       never silently discarded).
# Significance         : dual-gated on p_adj (BH per ranking x database) AND
#                       fold_enrichment clearing a minimum effect size. This
#                       matters because Fisher/hypergeometric tests gain
#                       enormous power at this N (~1.47M) -- a barely-there
#                       deviation from the null ratio (e.g. fold=0.95) can
#                       reach a "significant" p-value on its own. That is the
#                       general fact that any correctly-calibrated test's
#                       power grows with N, so p-value alone is not an
#                       adequate decision criterion at this scale. Pairing it
#                       with a fold-enrichment floor (posenrich_fold_thr,
#                       ~DAVID/PANTHER "moderate enrichment" convention)
#                       mirrors the dual-threshold design already used by the
#                       gene-level FCS report (p.adj AND p.perm).
# =============================================================================

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


# ── args ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Position-level Fisher-exact enrichment.")
    p.add_argument("--obs-scores", required=True,
                   help="position_scores.tsv (Gene, Position, CAAS_score, change_side)")
    p.add_argument("--gmt-dir", required=True)
    p.add_argument("--characterization", default=None,
                   help="characterization_layers.tsv (broad functional layers)")
    p.add_argument("--annot-file", default=None,
                   help="SCORING's fcs_stats.tsv (gene + flag_* columns) — cross-module "
                        "corroboration flags (gate_sig/gate_fdr/fade/rer/accum), reported "
                        "as the %% of distinct genes in each overlap carrying each flag")
    p.add_argument("--universe", required=True,
                   help="cleaned_background_main.txt gene list (postproc-surviving)")
    p.add_argument("--background", required=True,
                   help="caastools background.output (gene<TAB>tested positions); "
                        "restricted to --universe genes = the position background. "
                        "Always required: without it the Fisher-test background collapses "
                        "to scored positions only, invalidating the test. Use "
                        "--posenrich_background_file (params) to supply it in standalone runs.")
    p.add_argument("--cosmic-coverage", default=None,
                   help="cosmic_coverage_genes.txt from build_position_gmt.py — genes "
                        "COSMIC itself could annotate; restricts the background used "
                        "for the cosmic_orthogroups source only (see main() sources loop)")
    p.add_argument("--pai3d-coverage", default=None,
                   help="pai3d_coverage_genes.txt from build_position_gmt.py — genes "
                        "PAI3D itself could annotate; restricts the background used "
                        "for the pai3d_orthogroups source only")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-size", type=int, default=5,
                   help="min positions per set in background (GMT sources only)")
    p.add_argument("--max-size", type=int, default=0,
                   help="max positions per set in background (0 = no cap; GMT sources only)")
    p.add_argument("--char-fracs", default="0.25,0.10,0.05,0.01",
                   help="foreground top-fractions among scored positions, per ranking")
    p.add_argument("--position-lists-dir", required=True,
                   help="SCORING's published position_lists/slice_{top,bottom,global}{25,10,5,1}.tsv "
                        "dir (scoring_compute.R) -- the sole source of posenrich's foreground. SCORING "
                        "is a mandatory upstream dependency of posenrich (never optional), so there is "
                        "no local re-ranking fallback: a missing directory or slice file is a hard error.")
    p.add_argument("--padj-thr", type=float, default=0.15,
                   help="BH-adjusted p-value significance threshold")
    p.add_argument("--fold-thr", type=float, default=1.5,
                   help="minimum fold-enrichment magnitude for significance "
                        "(fold >= thr for enrichment, or <= 1/thr for depletion)")
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
    genes). Positions are in the same prot_ali_col space as the GMTs and scores.

    This file is always mandatory. Without it the Fisher-test background collapses
    to scored positions only (via the union in main()), making every foreground
    hit appear enriched — statistically invalid. A missing or sentinel value is
    therefore a hard error, not a silent fallback.
    In the _complete run scenario (CT skipped, reusing _exploratory outputs) the
    caller must supply the exploratory run's caastools/background.output via
    --posenrich_background_file / POSENRICH_BACKGROUND_FILE; the run_single.sh.j2
    template does this in its reuse_exploratory block.
    """
    if not background_file or background_file.startswith("NO_FILE") or not os.path.exists(background_file):
        sys.exit(
            f"[posenrich] ERROR: background file is required but was not supplied or "
            f"does not exist: {background_file!r}.\n"
            f"In a _complete run reusing a prior _exploratory pass, set "
            f"POSENRICH_BACKGROUND_FILE to the exploratory run's "
            f"caastools/background.output (the run_single.sh.j2 template handles "
            f"this automatically via its reuse_exploratory block).\n"
            f"In a standalone enrichment-only run, pass --posenrich_background_file "
            f"pointing to the relevant caastools/background.output."
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
    """Cross-module corroboration flags, per gene (SCORING's fcs_stats.tsv).
    Returns ({gene: {flag_name: bool}}, [flag_name, ...]); flag columns are
    picked up dynamically (any column starting with 'flag_'), mirroring
    fcs_enrich.R's grep('^flag_', names(stats)). Empty/missing file -> no flags,
    so callers simply skip the pct_* columns."""
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
    """Return {pos_id: CAAS_score} for a direction (global/top/bottom).
    global: every position, no change_side filter (kept as-is — no
    SCORING-canonical "global" convention exists to migrate to). top/bottom:
    change_side-filtered first (matching SCORING's 11.Scoring_report.Rmd
    convention: pos_top/pos_bottom = filter(change_side %in% ...)), so
    "none"-labeled positions are structurally excluded before any percentile
    is ever computed."""
    if direction == "global":
        sub = df
    elif direction == "top":
        sub = df[df["change_side"].isin(["top", "both"])]
    else:
        sub = df[df["change_side"].isin(["bottom", "both"])]
    return dict(zip(sub["pos_id"], sub["CAAS_score"]))


def load_published_foreground(position_lists_dir, direction, char_fracs):
    """Read SCORING's published position_lists/slice_<direction><pct>.tsv files
    (scoring_compute.R) as the foreground_by_frac for one direction. SCORING is
    a mandatory upstream dependency -- there is no local re-ranking fallback,
    so a missing directory or slice file is a hard error, not a silent
    degrade. Returns {frac: (n, set-of-pos_id)}, one entry per requested frac
    (an empty set for a frac whose published slice is legitimately empty, e.g.
    no scored positions in that direction -- distinct from the file being
    absent, which is an error)."""
    if not position_lists_dir or not os.path.isdir(position_lists_dir):
        sys.exit(f"[posenrich] --position-lists-dir not found or not a directory: "
                  f"{position_lists_dir!r}. SCORING's position_lists/ output is mandatory.")
    out = {}
    for frac in char_fracs:
        pct = int(round(frac * 100))
        path = os.path.join(position_lists_dir, f"slice_{direction}{pct}.tsv")
        if not os.path.exists(path):
            sys.exit(f"[posenrich] missing published slice: {path}. "
                      f"SCORING must publish position_lists/slice_{direction}{pct}.tsv.")
        df = pd.read_csv(path, sep="\t")
        pos_ids = set() if df.empty else set(df["Gene"].astype(str) + ":" + df["Position"].astype(str))
        out[frac] = (len(pos_ids), pos_ids)
    return out


# ── Fisher overlap test ───────────────────────────────────────────────────────
def annotate_overlap(overlap, annot, flag_names):
    """Cross-module corroboration: of the DISTINCT genes contributing to this
    overlap (not the raw position count — a single gene contributing many
    overlapping positions must not dominate the percentage), what fraction
    carry each flag_* column from SCORING's fcs_stats.tsv? Returns
    {pct_<flag w/o 'flag_' prefix>: value}, NaN when the overlap is empty."""
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


def run_fisher_for_terms(terms, descs, foreground_by_frac, bg_set, N, min_size, max_size,
                         annot=None, flag_names=None):
    """Fisher-exact overlap test for every term in one database against one
    ranking, at each fixed top-fraction cutoff already precomputed for this
    ranking (foreground_by_frac: frac -> (foreground_size, foreground_set)).
    Background (N) is always the full honest tested position pool — never
    restricted to the scored subset, so tested-negative positions are never
    silently dropped from the null. min_size/max_size (set size within the
    background) apply only to GMT-style term collections; pass 0/0 to skip
    that filter for the broad characterization layers. annot/flag_names (from
    load_annot) add pct_<flag> cross-module corroboration columns when given."""
    out = []
    for term, members in terms.items():
        M = set(members) & bg_set
        mm = len(M)
        if min_size and mm < min_size:
            continue
        if max_size and mm > max_size:
            continue
        desc = descs.get(term, term)
        for frac, (nn, foreground) in foreground_by_frac.items():
            overlap = foreground & M
            k = len(overlap)
            expected = nn * mm / N if N else np.nan
            fold = (k / expected) if expected and expected > 0 else np.nan
            table = [[k, nn - k], [mm - k, N - nn - (mm - k)]]
            try:
                _, p = fisher_exact(table, alternative="two-sided")
            except Exception:
                p = np.nan
            row = dict(pathway=term, description=desc, top_frac=frac,
                      layer_size=mm, foreground_size=nn, observed=k,
                      expected=expected, fold_enrichment=fold, p_value=p,
                      direction=("enriched" if fold and fold >= 1 else "depleted"),
                      # N as actually used for THIS row's Fisher test (respects
                      # per-database coverage restriction, e.g. COSMIC/PAI3D) --
                      # exposed so report consumers can run their own honest
                      # position-level hypergeometric tests (e.g. the Overall
                      # dotplot's position-native FADE-overlap component)
                      # against the exact same universe, not a re-derived guess.
                      background_n=N)
            row.update(annotate_overlap(overlap, annot or {}, flag_names or []))
            # Kept only transiently for the leading-edge export; stripped in
            # main() before the row goes into the main wide table (sig is not
            # known yet here -- it depends on the BH pass over the whole batch).
            row["_overlap"] = overlap
            out.append(row)
    return out


def restrict_background_to_coverage(bg_set, coverage_genes):
    """Restrict a Gene:Position background set to positions whose gene is in
    `coverage_genes`. Used for GMTs derived from external, incompletely-covered
    databases (COSMIC, PAI3D): a gene that database could never have annotated
    should not count as a "true negative" in that database's own enrichment
    test, or the test is diluted by genes it structurally could never see."""
    return {p for p in bg_set if p.rsplit(":", 1)[0] in coverage_genes}


def restrict_foreground(foreground_by_frac, restricted_bg):
    """Recompute (size, set) foreground pairs against a restricted background.
    Required alongside restrict_background_to_coverage: the Fisher 2x2 table
    needs foreground_size <= background_size, so once a database's background
    (N) shrinks to a coverage subset, its foreground must shrink to the same
    subset too — not just the background/mm/expected terms — or the table
    becomes inconsistent (foreground_size can exceed the shrunken N)."""
    return {frac: (len(fg & restricted_bg), fg & restricted_bg)
            for frac, (nn, fg) in foreground_by_frac.items()}


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
    bg_set = set(background)
    N = len(background)
    print(f"[posenrich] background positions: {N} (background.output∩universe={n_bg_only}) "
          f"| scored: {(obs['CAAS_score']>0).sum()}", flush=True)

    # Per-database background restriction for GMTs derived from external,
    # incompletely-covered databases (COSMIC, PAI3D) — see
    # restrict_background_to_coverage()'s docstring. Every other GMT (internal
    # deterministic computations with complete coverage by construction) keeps
    # the full honest background untouched.
    coverage_restricted_bg = {}
    for db_name, coverage_path in (("cosmic_orthogroups", args.cosmic_coverage),
                                    ("pai3d_orthogroups", args.pai3d_coverage)):
        if coverage_path and os.path.exists(coverage_path):
            coverage_genes = load_universe_genes(coverage_path)
            restricted_bg = restrict_background_to_coverage(bg_set, coverage_genes)
            coverage_restricted_bg[db_name] = restricted_bg
            print(f"[posenrich] {db_name}: background restricted to "
                  f"{len(coverage_genes)} coverage genes "
                  f"({len(restricted_bg)}/{N} positions)", flush=True)

    gmts = load_gmts(args.gmt_dir)
    char_layers = read_charset(args.characterization)
    annot, flag_names = load_annot(args.annot_file)
    if flag_names:
        print(f"[posenrich] cross-module flags: {', '.join(flag_names)} "
              f"({len(annot)} genes annotated)", flush=True)
    # Unify GMT-style term collections and the broad characterization layers
    # into one sources dict: database -> (terms, descs, apply_size_filter).
    # Both are structurally the same object (name -> position members), so one
    # Fisher-loop handles both; only the min/max-size filter differs (GMT terms
    # are size-bounded as before, the 4 broad layers never were).
    sources = {db: (terms, descs, True) for db, (terms, descs) in gmts.items()}
    if char_layers:
        char_terms = {name: members for name, (desc, members) in char_layers.items()}
        char_descs = {name: desc for name, (desc, members) in char_layers.items()}
        sources["characterization"] = (char_terms, char_descs, False)
    print(f"[posenrich] sources: {len(sources)} ({', '.join(sources)})", flush=True)

    # 'global' restored alongside top/bottom: FCS scoring (12.FCS_general_report.Rmd)
    # and AMI (13.AMI_analysis.Rmd) both carry a global/top/bottom triad, and
    # posenrich should offer the same non-directional view (every position,
    # no change_side filter) rather than being the odd one out. direction_scores()
    # already implements "global" -- see its docstring -- it just wasn't in
    # this loop.
    directions = ["global", "top", "bottom"]
    rows = []
    leading_edge_rows = []
    for direction in directions:
        obs_scores = direction_scores(obs, direction)
        n_scored = sum(1 for s in obs_scores.values() if s > 0)
        if n_scored == 0:
            continue
        # Foreground comes exclusively from SCORING's published
        # position_lists/slice_<direction><pct>.tsv (the same canonical
        # top/bottom/global x 25/10/5/1% cutoff every other report/module
        # reads) -- no local re-ranking; see load_published_foreground()'s
        # docstring for why there is no fallback.
        foreground_by_frac = load_published_foreground(args.position_lists_dir, direction, char_fracs)
        print(f"[posenrich] {direction}: {n_scored} scored positions ranked "
              f"(foreground from published position_lists)", flush=True)

        for db, (terms, descs, apply_size_filter) in sources.items():
            if db in coverage_restricted_bg:
                db_bg_set = coverage_restricted_bg[db]
                db_N = len(db_bg_set)
                db_foreground_by_frac = restrict_foreground(foreground_by_frac, db_bg_set)
            else:
                db_bg_set, db_N, db_foreground_by_frac = bg_set, N, foreground_by_frac
            res = run_fisher_for_terms(
                terms, descs, db_foreground_by_frac, db_bg_set, db_N,
                args.min_size if apply_size_filter else 0,
                args.max_size if apply_size_filter else 0,
                annot=annot, flag_names=flag_names)
            if not res:
                continue
            # BH per (ranking, database): a handful of huge broad layers and up
            # to thousands of small GMT terms carry very different multiple-
            # testing burdens, so one global family across everything would
            # conflate them.
            padj = bh_adjust([r["p_value"] for r in res])
            for i, r in enumerate(res):
                r["p_adj"] = padj[i]
                # Scored-position population for THIS direction (constant across
                # every term/database in this direction's loop iteration) --
                # top_frac cutoffs are fractions of n_scored, NOT of background_n
                # (the honest tested-position pool, which also includes untested/
                # zero-score positions) -- exposed so report consumers computing
                # their own top-fraction-recurrence tests use the matching universe.
                r["n_scored"] = n_scored
                # No fold-enrichment gate here by design — p_adj alone decides
                # significance; fold_enrichment stays in the output for the
                # caller to filter on downstream if they want a magnitude cutoff.
                r["sig"] = bool(r["p_adj"] < args.padj_thr)
                overlap = r.pop("_overlap")
                if r["sig"]:
                    for pos_id in sorted(overlap):
                        leading_edge_rows.append(dict(
                            ranking=direction, database=db, pathway=r["pathway"],
                            top_frac=r["top_frac"], gene=pos_id.rsplit(":", 1)[0],
                            gene_position=pos_id))
                rows.append(dict(ranking=direction, database=db, **r))

    results = pd.DataFrame(rows)
    if not results.empty:
        results = results.sort_values(["p_adj", "p_value"], na_position="last")
    out_path = os.path.join(args.output_dir, "posenrich_characterization.tsv")
    results.to_csv(out_path, sep="\t", index=False)
    print(f"[posenrich] wrote {out_path} ({len(results)} rows)", flush=True)

    # Observed gene:position members of the overlap, one row per member, for
    # significant (ranking, database, pathway, top_frac) terms only — mirrors
    # fcs_leading_edge.tsv in the gene-level FCS report.
    leading_edge = pd.DataFrame(leading_edge_rows)
    le_path = os.path.join(args.output_dir, "posenrich_leading_edge.tsv")
    leading_edge.to_csv(le_path, sep="\t", index=False)
    print(f"[posenrich] wrote {le_path} ({len(leading_edge)} rows)", flush=True)


if __name__ == "__main__":
    main()
