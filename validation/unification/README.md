# Unification safety net (Phase 0)

Checks that the permulation-null code path reproduces the observed results when it is fed the
real labeling (`b_0`). Spec: the "Unified CAAS core" plan.

- `compare_b0.py`: observed vs `b_0` at checkpoints A (discovery rows), B (post-filter survivors),
  C (`asr_path_score`), D (`caas_row`/`CAAS_score`), E (gene score). Exact key-set equality, |delta| <= 1e-12.
  `test_compare_b0.py` is its plumbing test (synthetic run).
- `--caas_b0_diagnostic true` makes the pipeline write the `b_0` slice to `<outdir>/caas_permulation/b0/`
  (never part of the null). Needs `ct_disambig_perms_batch_size = 1`.
- `golden/pepc_c4_complete/`: frozen observed outputs of the Tier 1 PEPC genotypic run
  (`validation/tier1/output/pepc/results/c4_complete`, code at 7892848).
- `baseline/`: harness reports on the code before unification (see below).

Usage:
```bash
python validation/unification/compare_b0.py --run <results_dir> [--out report.json]
```

## Baselines (code before unification; `baseline/*.json`)

`compare_b0.py` on runs that also replay b_0 through the null path (`--caas_b0_diagnostic true`).

| Fixture | A discovery | B survivors | C asr_path_score | D CAAS_score | E gene score |
|---|---|---|---|---|---|
| PEPC, FOP (100 hyp, 1 gene) | pass | pass | pass (max delta 1.1e-16) | pass (1.1e-16) | pass |
| PEPC, plain null (H1 only) against the FOP-pooled observed | pass (H1 only) | fail | fail | fail | fail |
| Cancer toy, FOP (12 hyp, 28 genes, batched) | pass | pass | pass (1.1e-16) | pass (1.2e-16) | pass |

- PEPC plain fails by design: the observed score pools 100 hypotheses, the plain null replays H1 only.
  That configuration (`caas_perms_fop=false` with `multi_hypothesis=true`) no longer exists: the null now
  mirrors the observed design, following `multi_hypothesis`; the row is kept as the historical baseline.
- What the fixtures do not exercise: gene removal (no gene was removed on either side; M6),
  `remove_caas_clusters=false` (M7), and, if the alignments have no missing data, `miss_pair` (M1).
  Clusters are exercised (23 positions, one gene, identical on both sides).
- Local test of the train (ctrain) grain, minlen 3 / maxcaas 0.7 on `discovery.tab`: toy, union and per-hypothesis
  flag the same 23 positions (all hypotheses share one position set); PEPC (100 hyp), union flags 25 and
  per-hypothesis 15, with 10 positions flagged only on the union.

## Post-processing (`core/postproc.py`)

One implementation of cluster trains and gene removal serves the observed chain
(`filter_caas_clusters-param.py`, `filter_caas_genes.py`) and the null (`gene_wrapper.py`, passes A and B).

- Trains: `train_flags` files positions under a key (the null uses `(cycle, caap_group)` per gene; the observed
  scripts use `(caap_group)` per gene over the pooled rows). This is the single switch point for the train grain.
- Gene removal: per labeling and `caap_group`, over the pooled scored rows, with no per-hypothesis grain.
  `dubious` = distinct positions above `Q3 + k*IQR` and at least one train position, calibrated over all units;
  `extreme` = density above the percentile, calibrated over units whose gene has a positive length.
- `remove_caas_clusters=false` reaches the null as `--keep-clusters`: train positions stay in the scored pool and
  still count for the dubious test.
- Tests: `test_postproc.py` (units, pandas oracle for the thresholds, mutation-sensitive boundaries) and
  `test_postproc_wiring.py` (observed script and null pass B agree with the core).
- Checked on the toy (stored code-before outputs): `clust` recomputed on 451 `b_0` positions, 0 mismatches;
  pass B (`reaggregate_perm_scores.py`) reproduces `gene_cycle_scores`, `perm_pos_quantiles`, `perm_pos_sample`
  and `perm_pos_cycle_caas` exactly; `filtered_discovery.tsv` reproduced byte for byte.
- Not exercised by any fixture: a gene actually removed at scale, genes without a length, and the effect of
  moving the observed gene removal from hypothesis to labeling grain (toy hypotheses share one position set,
  PEPC has one gene).

## Scores (`core/scores.py`)

Position score, side collapse and gene score are defined once and called by the null (pass B of
`gene_wrapper.py`); the observed side adopts them in the scoring step.

- Position score: mean of `caas_row` (= `asr_path_score`) over the schemes that scored it, per side. The sum is
  `math.fsum` (correctly rounded, independent of scheme order).
- Directions: `top` / `bottom` use only that side; `all` keeps one entry per position, its best side.
- Gene score: `(#{pool <= max + 1e-12} / |pool|) ** n`. The tolerance (`TIE_TOL`) makes ties deterministic: position
  scores are means of a few values and the pool is heavily tied (toy, 1000 genes: 3016 rows, ~887 distinct values),
  so means equal in exact arithmetic that differ by an ulp would otherwise move the count by several positions.
  On that toy the gaps between distinct scores were either ~1e-17 or >= 1e-7. It is None (written NA) when the gene has no scored position in
  the direction or the pool is empty. `scoring_caas_perms.R` fills those cells with 0 when it builds the dense
  genes x cycles matrix, so `caas_perms.rds` and the FCS null do not change; only `gene_cycle_scores.tsv` shows
  NA where it used to show 0.0.
- Checked: `size_adj_max` reproduces `gene_caas_score{,_top,_bottom}` written by the R pipeline (PEPC golden and
  toy, abs 1e-12). On the stored toy `b_0` detail, pass B outputs equal the previous ones (delta 0) except 10
  cells that went from 0 to NA (8 top, 2 bottom). The position score differs from R's `mean()` in the last bit
  on 46 of 165 rows (1.2e-16); it disappears once the observed consumes these scores.

## Observed scores (`observed_core_scores.py`, `scoring_compute.R`)

The observed position and gene scores are the `b_0` slice of the same functions as the null.
`SCORING_COMPUTE` runs `observed_core_scores.py` (stdlib + `core/scores.py`) on `filtered_discovery.tsv`, then
`scoring_compute.R` reads `core_positions.tsv` / `core_genes.tsv` and integrates them (FADE, RER, accumulation,
`p.emp`, BH, SAM). R no longer computes `mean(caas_row)` or `size_adj_max`.

- `p.emp` and SAM count null values within `TIE_TOL` of the observed score as ties (`>=`); the constant is defined in
  `core/scores.py` and repeated in `scoring_compute.R` (a test keeps them equal).
- The observed chain reads floats exactly (`float_precision="round_trip"` in `prepare_postproc_input.py`,
  `filter_caas_genes.py`); the pandas default alters ~1/3 of doubles by an ulp. On the local toy, positions
  differing from `b_0` in the last bit went from 50 to 16 of 165. The remainder appears only in rows pooled over
  several hypotheses (`n_hypotheses` >= 3); its origin (PSS weights of the null's `b_0`) was not traced.
- Checked: `scoring_compute.R` on the toy (50 genes) gives position and gene tables equal to the previous script
  (Δ <= 1.1e-16, `p.emp` and `p.adj_*` unchanged). Not checked with a null that has ties at scale (toy 1000 genes).

## Null per-cycle position scores (`perm_pos_cycle_caas.tsv.gz`)

Columns: `Gene, Position, side, cycle, caas_score, n_schemes`. `caas_score` is the `core.scores` position score of
the null cycle (empty when no scheme scored it), written once by `_finalize_perm_scores`. `scoring_compute.R`
(`p.emp`, SAM), `posenrich_enrich.py` and `posenrich_prep_caas_null.py` read it and no longer divide. A file
without `caas_score` (an earlier null) is rejected with a message to regenerate it. `readr` misparses ~13 % of
17-digit doubles by an ulp, so R reads `caas_score` as text and converts with `as.numeric`; pandas readers use
`float_precision="round_trip"`.

## Order-independent pooling (`fop_pool.pool_domains`)

The hypothesis and domain sums in `pool_domains` use `math.fsum`, so the pooled score is the same whatever order
the hypotheses arrive in (the null feeds them in labeling arrival order, the observed in its own order). A naive
sum of M >= 3 terms depends on that order; with M = 2 it does not, which matches where observed and null `b_0`
differed in the last bit (rows pooled over >= 3 hypotheses). `test_pool_domains_order.py` shuffles the hypotheses
(M = 3, 12, 100) and requires identical bits. Against the frozen PEPC master (100 hypotheses) the float columns now
differ by <= 1.7e-15 and every other cell is identical, so `test_observed_pepc.py` compares floats with a
1e-12 tolerance. Whether observed and `b_0` now agree bit for bit is checked on the next cluster run.

## Harness details (`compare_b0.py`)

- C reads with `float_precision="round_trip"` and reports `n_bitwise_different` (informational: rows whose
  `asr_path_score` differs in any bit between the observed master and the `b_0` detail; the pass criterion stays
  |delta| <= 1e-12).
- E requires the NA pattern of the gene scores to match: a gene with no scored position in a direction is NA on both
  sides. A null that writes 0 there fails.
