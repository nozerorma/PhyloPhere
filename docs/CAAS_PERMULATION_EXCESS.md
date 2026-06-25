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

## Alignment with RERconverge (validated against msae210, the categorical paper)

The categorical RERconverge paper (MBE 2024, msae210; vignette + supplement in
~/Downloads/msae210_supplementary_data/) confirms our FCS design IS their method:
- Engine: `fastwilcoxGMTall(getStat(...), annotlist, outputGeneVals=T)` — same as ours.
- **One-sided** `alternative="greater"` for magnitude/omnibus tests ("only more positive
  statistics indicate strength"); two-sided for signed pairwise. ⇒ TODO: set
  alternative="greater" in `fcs_run_ranking` for the non-negative rankings
  (CAAS/FADE/accum, RER accel/decel zero-floored); keep two-sided for RER signed global.
- Pathway-level permulations (permute phenotype → recompute RER → re-enrich → empirical
  p) — exactly what fcs_enrich.R does by re-enriching the corStat perm matrix.
- Significance = **adjusted p < 0.15 ∩ permulation p < 0.025** (relaxed BH because
  underpowered + strict perm-p as the non-independence filter; their own example has NO
  pathway < 0.05 adj-p). p.perm is NOT BH'd (matches fcs_enrich.R). ⇒ adopt 0.15/0.025
  as the default reporting convention.
The vectorization (below) must reproduce fastwilcoxGMT EXACTLY incl. the one-sided
alternative and ties.method="average" — guarded by the equivalence test.

## FCS permulation performance + progress — DONE (2026-06-25)

Implemented in `fcs_enrich.R`: the permulation null is now a vectorized sparse
`M %*% Rk` per GMT (`fcs_null_enrichstat_vectorized`), replacing 3×N
`fastwilcoxGMTall` calls. Equivalence to `fastwilcoxGMT` proven at max|Δstat|=0
by `fcs_enrich_equivtest.R` (clean/ties/NA/num.g/multi-col). Benchmark (16k genes,
16 GMTs, 35,520 sets): null **~100 s for N=1000** vs ~5 h before; total fcs_run_all
~95 s (the rest is the N-independent observed enrichment). One-sided `alternative`
auto-detected per ranking (signed→two.sided, magnitude→greater) incl. recomputed
one-sided observed analytic p. Progress → stderr + `fcs_progress.log`. The mclapply
ncores path was removed entirely (no longer needed), so the ncores env bug is moot.

Original analysis (kept for context):

## FCS permulation performance + progress (PREREQUISITE — affects RER too)

The current FCS permulation path (`subworkflows/ENRICHMENT/local/fcs_enrich.R`) is
the bottleneck for the existing RER run (≈4 h on 16k genes) and will be reused by
the CAAS excess null — so fix it before scaling either up.

**Root causes**
1. `fcs_get_enrich_perms_parallel` calls `RERconverge::fastwilcoxGMTall` 3 rankings ×
   N permulations times (`fcs_enrich.R:149`), re-ranking 16k genes + re-looping all
   ~30k sets every call. `fcs_permpvalenrich_vectorized` only vectorized the final
   p-step, not null generation.
2. ncores fallback bug (`fcs_enrich.R:258`): `Sys.getenv("SLURM_CPUS_PER_TASK","1")`
   → off-SLURM defaults to "1"; the `detectCores()` fallback never fires because 1L
   is not NA/<1. Local runs silently use ONE core. Fix:
   `n <- Sys.getenv("SLURM_CPUS_PER_TASK"); ncores <- if (nzchar(n)) as.integer(n) else parallel::detectCores()`.

**The fix — vectorize (stay in R; do NOT switch language)**
The Wilcoxon-AUC for a set depends only on its members' rank-sum. Replace the
per-permulation `fastwilcoxGMTall` calls with:
- `Rk = matrixStats::colRanks(corStat, ties.method="average")` — rank each perm col once.
- `M` = sparse set×gene membership matrix (`Matrix`), built once.
- `M %*% Rk` → rank-sums for all sets × all perms in one BLAS call.
- Analytic AUC: `U = ranksum - n_S(n_S+1)/2; stat = U/(n_S*n_other) - 0.5`; empirical p.
- Compute the OBSERVED stat with the same formula for the comparison (keep
  fastwilcoxGMTall's value for display only) so observed/null share a scale.
Expect hours → seconds-minutes. Memory ~130 MB (ranks) + ~240 MB (rank-sums) at N=1000.

**Progress tracking (R-side, Nextflow-log friendly)**
Chunk the N perms; after each chunk emit a timestamped line to stderr AND append to a
tailable `fcs_progress.log`:
`[FCS][<ranking>] perms 400/1000 (40%) | 12.3/s | elapsed 0.5m | ETA 0.8m`.
(`pbmcapply::pbmclapply` is a drop-in alternative but a bar renders poorly in
non-TTY pipeline logs.) Once vectorized, a per-ranking message largely suffices.

## Effort estimate

Weekend prototype IF the purger headless-mode is clean. Dominant work = null-mode
`gene_wrapper` (small, seam exists) + headless purger (unknown coupling). Discovery
is cheap/parallel. FCS hand-off is done.

## Decision deferred

Wait and see how RER FCS permulation behaves first (already running). If the
analytic-vs-permulation gap is informative and worth the compute, build this for
CAAS. Accumulation stays BH-only (per-randomization position dumps too heavy).
