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
| PEPC, plain (`caas_perms_fop=false`) | pass (H1 only) | fail | fail | fail | fail |
| Cancer toy, FOP (12 hyp, 28 genes, batched) | pass | pass | pass (1.1e-16) | pass (1.2e-16) | pass |

- PEPC plain fails by design: the observed score pools 100 hypotheses, the plain null replays H1 only.
- What the fixtures do not exercise: gene removal (no gene was removed on either side; M6),
  `remove_caas_clusters=false` (M7), and, if the alignments have no missing data, `miss_pair` (M1).
  Clusters are exercised (23 positions, one gene, identical on both sides).
- Local test of the train (ctrain) grain, minlen 3 / maxcaas 0.7 on `discovery.tab`: toy, union and per-hypothesis
  flag the same 23 positions (all hypotheses share one position set); PEPC (100 hyp), union flags 25 and
  per-hypothesis 15, with 10 positions flagged only on the union.
