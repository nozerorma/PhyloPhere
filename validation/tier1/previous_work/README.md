# Previous work

Material produced with the pipeline that preceded the unified CAAS core. It is kept for comparison and as a record, not as a current result.

| Path | Content | Status |
|---|---|---|
| `output/pepc_pre_unification/` | results of the two PEPC runs (genotypic `c4`, phenotypic `c4_phenotypic`) | the runs in `../output/pepc/` reproduce them: `compare_pepc_runs.py --require-equal` reports equal positions, scores, `p.emp`, BH and SAM in both traits |
| `input/pepc_negctrl/` | the 20 negative-control traits (`my_traits.tsv`), their generator `build_negctrl_traits.R`, the launcher `run_negctrl_local.sh` and the summary script `analyze_negctrl.py` | not repeated with the current pipeline |
| `output/pepc_negctrl/` | results of those controls and their summary `negctrl_summary.tsv` | produced by the previous pipeline; section 6e of `../reports/pepc_genotypic_vs_phenotypic.md` reports them |

Output directories are ignored by git.

To repeat the controls with the current pipeline, run `input/pepc_negctrl/run_negctrl_local.sh` (about 10 minutes per control on a laptop). Its results go to `../output/pepc_negctrl/`, not here, so that old and new results stay apart; `analyze_negctrl.py <results_dir> <overlap_tsv> <out_tsv>` summarizes either set.
