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
- Gene score: `(#{pool <= max} / |pool|) ** n`. It is None (written NA) when the gene has no scored position in
  the direction or the pool is empty. `scoring_caas_perms.R` fills those cells with 0 when it builds the dense
  genes x cycles matrix, so `caas_perms.rds` and the FCS null do not change; only `gene_cycle_scores.tsv` shows
  NA where it used to show 0.0.
- Checked: `size_adj_max` reproduces `gene_caas_score{,_top,_bottom}` written by the R pipeline (PEPC golden and
  toy, abs 1e-12). On the stored toy `b_0` detail, pass B outputs equal the previous ones (delta 0) except 10
  cells that went from 0 to NA (8 top, 2 bottom). The position score differs from R's `mean()` in the last bit
  on 46 of 165 rows (1.2e-16); it disappears once the observed consumes these scores.
