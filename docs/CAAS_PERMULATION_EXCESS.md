# CAAS permulation-excess enrichment — implementation plan

> Status: **planned** (not started). Captured 2026-06-25.
> Goal: give CAAS a genome-wide *excess* null so FCS pathway enrichment can be
> permulation-corrected (a `p.perm` column), exactly like RER already has — by
> running a small N (~1000) of permuted phenotypes through the FULL position pool
> and feeding the resulting null gene-rankings into the existing FCS permulation
> path in `subworkflows/ENRICHMENT/local/fcs_enrich.R`.

## Why it's feasible (the load-once invariant)

ASR ancestral posteriors are **phenotype-invariant**: `load_precomputed_asr()`
(`subworkflows/CT_DISAMBIGUATION/local/src/utils/gene_wrapper.py:341`) is a pure
function of (alignment, gene tree). The phenotype enters only *downstream*, in the
pure scoring functions of `src/convergence/path_scores.py`, as which species form
each pair and which tip states sit on each side. So a permuted phenotype changes
**which nodes are walked**, never the posteriors themselves.

⇒ Per gene: load `node_posteriors` ONCE, then replay N permuted labelings over the
cached object in memory. 1000× JSON re-parse collapses to 1× I/O + 1000× CPU
traversal. Because the null reuses the *same* `path_scores` functions as the real
run, the null is calibrated by construction (no code-path drift).

## Architecture

Per gene (one gene in memory at a time):
1. `load_precomputed_asr()` → `node_posteriors` (once).
2. For each of N (~1000) permuted FG/BG labelings (reuse the same permulation
   resamples already generated for CAAS significance — `permulations.R`
   `resample_NNN.tab`, but applied to the FULL pool, not positives-only):
   a. **Discovery (full pool)** under the permuted labeling → participating
      positions/pairs. Cheap, embarrassingly parallel.
   b. **Postproc: cluster-purger only** → collapse clusters, emit minimal
      gene×score table. No reports, no JSON dumps.
   c. **Score** via existing `path_scores` functions on cached posteriors →
      per-gene null `CAAS_score` (same 0.90-quantile aggregation as the real run,
      `scoring_compute.R:636`).
3. Emit a **genes × N** null-score matrix (storage-light: gene scores only, never
   per-position — this is what makes CAAS lighter than the accumulation version,
   which needed `export_individual` position dumps).

Hand-off: genes×N matrix → `fcs_enrich.R` permulation path
(`fcs_get_enrich_perms_parallel` / `fcs_permpvalenrich_vectorized` / the `p.perm`
column, `fcs_enrich.R:134-314`). Enrichment side is already implemented — free.

## Components to build

1. **Null-mode `gene_wrapper`** — new entry point / flag in
   `subworkflows/CT_DISAMBIGUATION/local/src/utils/gene_wrapper.py`: load
   `node_posteriors` once, accept a list of N labelings, return N per-gene scores.
   The load/score seam already exists at `node_posteriors`; this is the main work.
2. **Full-pool discovery in null mode** — discovery currently runs the real
   labeling; needs to run per permuted labeling over the full position pool (not
   the positives-only fast path used for significance permulations today). MUST
   mirror the real run's directives exactly (grouping scheme, thresholds) or the
   null is miscalibrated.
3. **Headless cluster-purger mode** in CT_POSTPROC — emit gene×score table only,
   skip report rendering and any JSON-dump side effects.
4. **Permulation driver (Nextflow)** — orchestrate N labelings → genes×N matrix.
   Mirror the RER perms plumbing in `workflows/rerconverge.nf` (perms_file →
   `RER_FCS_REPORT`). Add a `caas_permulation_enrichment` toggle (default false),
   paralleling `params.rer_permulation_enrichment`.
5. **Wire matrix → FCS** — pass the genes×N null matrix into the CAAS
   `SCORING_FCS_REPORT` permulation path (same shape RER's perms take).

## Calibration / stats notes

- Null scores MUST be computed identically to observed (same code path). Reuse
  `path_scores` functions verbatim; do not reimplement.
- N≈1000 → p.perm floor = 1/(N+1) ≈ 0.001; p.perm takes only N+1 discrete values
  and ties at the floor → coarse robustness filter, not a fine scale. Rank within
  the floored set by the AUC `stat` / analytic `p.adj`.
- p.perm is NOT BH-corrected in `fcs_enrich.R` (only analytic `p.val`→`p.adj` is).
  Defensible filter = `p.adj < 0.05` ∩ `p.perm` at/near floor. If BH-ing p.perm
  later, the floor bites (need a co-floored block) → only then consider larger N.
- Directionality is free: a permuted labeling reshuffles top/bottom, so directional
  nulls (top={top,both}, bottom={bottom,both}) drop out of the same loop.

## Open question to resolve first

Confirm the cluster-purger can run **headless** (gene table only, no report/JSON
side effects). If it currently couples to reporting, that's the one real plumbing
task before prototyping.

## Effort estimate

Weekend prototype IF the purger headless-mode is clean. Dominant work = null-mode
`gene_wrapper` (small, seam exists) + headless purger (unknown coupling). Discovery
is cheap/parallel. FCS hand-off is done.

## Decision deferred

Wait and see how RER FCS permulation behaves first (already running). If the
analytic-vs-permulation gap is informative and worth the compute, build this for
CAAS. Accumulation stays BH-only (per-randomization position dumps too heavy).
