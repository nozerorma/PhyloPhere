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
