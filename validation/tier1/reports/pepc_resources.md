# Tier 1 PEPC: Resource & Runtime Report

Two PhyloPhere runs on the shared PEPC fixture (`validation/tier1/input/pepc/`), one per trait definition, identical except for the trait and the pruning it implies, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` |

## Run configuration

`output/pepc/work/<run>_nxf_run/params.json`, compared key by key. Besides run-scoped paths (`outdir`, `gmt_dir`, `tax_id`) the two runs differ only in the trait (`traitname` / `secondary_trait` swapped) and the pruning parameters (`prune_data`, `prune_list`, `prune_list_secondary`: off vs the *Eleocharis* list). Shared settings relevant to the results: `seed` 1998, `caas_full_perms` 1000, `max_fop` 100, `min_contrasts` 3, `filter_minlen` 3, `filter_maxcaas` 0.7, `scoring_p_emp_thr` 0.05, `fade_direction` top, `fade_background_scope` all, `fade_internal_nodes` all_descendants (from `conf/fade.config`; not written to `params.json` by the runner), `fcs_enabled` false.

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Disk (`/home`) | 232 GB total, 55 GB available |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` |

RAM, swap and disk are a snapshot at report time.

## PhyloPhere runs: totals

Both runs executed every task (`cachedCount=0`), sequentially, while the negative-control batch (`output/pepc_negctrl/`, one Nextflow run at a time, up to 8 CPUs) was running on the same machine. Wall times and %CPU below therefore include contention and are an upper bound on the pipeline's own cost.

| | genotypic | phenotypic |
|---|---|---|
| Nextflow run name | `c4_complete_8a176b39` | `c4_phenotypic_complete_97c9f7ca` |
| Launch → completion (`.nextflow.log`) | 2026-09-27 23:14:53 → 23:32:35 (17m 42s) | 2026-09-27 23:33:40 → 23:51:26 (17m 46s) |
| Processes succeeded / failed | 37 / 0 | 37 / 0 |
| Nextflow `succeedDuration` (cumulative task time) | 2h 20m 52s | 2h 11m 43s |
| `peakRunning` / `peakCpus` | 5 / 26 | 5 / 26 |
| `peakMemory` (requested, not resident) | 28 GB | 28 GB |
| Largest single-task `peak_rss` | 3.0 GB (`CAAS_PERMS_DISAMBIGUATE`) | 2.4 GB (`CAAS_PERMS_DISAMBIGUATE`) |
| Permulation cycles × hypotheses | 1000 × 100 | 1000 × 100 |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; resident usage is bounded by the per-task `peak_rss` column.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | 1m 5s, 1.9 GB | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION, NAME_CURATION:TREE_CLEANUP}` | run | not run | non-pruning path |

## Dominant cost centres

| process | genotypic realtime | %cpu | phenotypic realtime | %cpu |
|---|---|---|---|---|
| `CT:RESAMPLE` | 7m | 377.3 % | 5m 9s | 412.5 % |
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 5m 11s | 22.1 % | 6m 1s | 19.3 % |
| `FADE:FADE_BATCHED` | 3m 47s | 77.1 % | 3m 13s | 79.7 % |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 1m 14s | 101.0 % | 1m 35s | 100.5 % |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 1m 13s | 74.3 % | 41.9 s | 77.9 % |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 1m 8s | 72.2 % | 1m 19s | 65.5 % |

- **`CAAS_PERMS_DISAMBIGUATE` runs at about 20 % CPU** in both runs: most of its wall time is not spent computing. Whether this is I/O, worker start-up or contention with the concurrent negative-control batch was not profiled.
- **`CT:RESAMPLE` breakdown** (from its `.command.log`):

| phase | genotypic | phenotypic |
|---|---|---|
| pool harvest, 1000 Tier-1 cycles | 1m 24s (1012 draws) | 26 s (1002 draws) |
| design matching: FOP harvest of every candidate | 2m 27s + 27 s top-up (1147 candidates, 89.0 % reach 100 hypotheses) | 2m 11s + 12 s (1066 candidates, 94.4 %) |
| FOP mirror, 1000 cycles × 100 hypotheses, 6 workers | 2m 32s | 2m 16s |

Design matching harvests each candidate cycle's FOP hypotheses once to count them and again in the mirror; together they are the largest share of `RESAMPLE` in both runs. Differences between the runs combine the 4-pair vs 3-pair design, the tip set, the lower genotypic match rate (more candidates to harvest) and machine contention; they were not separated.

## Per-process trace, side by side

Sources: `results/<run>/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `CONTRAST_SELECTION:NAME_CURATION:TREE_CLEANUP` | 313ms | 71.5% | 28 MB | 149ms | 87.4% | 26.9 MB |
| `CONTRAST_SELECTION:REPORTING:NAME_CURATION:TREE_CLEANUP` | 304ms | 70.8% | 28.6 MB | not run | |  |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 12.9s | 71.5% | 432.8 MB | not run | |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 1m 7s | 96.8% | 2 GB | not run | |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 20.9s | 94.2% | 902.6 MB | 13.9s | 102.8% | 903.5 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 95ms | 83.3% | 3.8 MB | 106ms | 74.5% | 18.2 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 55.5s | 90.8% | 1.2 GB | 35.1s | 100.9% | 1.2 GB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 56ms | 60.0% | 0 | 32ms | 65.5% | 0 |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 289ms | 88.6% | 7.9 MB | 446ms | 65.3% | 7.9 MB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 215ms | 94.6% | 7.9 MB | 172ms | 98.6% | 8 MB |
| `CT:DISCOVERY_BATCHED` | 3.2s | 84.6% | 68.3 MB | 2.1s | 99.2% | 68 MB |
| `FADE:FADE_BATCHED` | 3m 47s | 77.1% | 477.9 MB | 3m 13s | 79.7% | 484.6 MB |
| `CT:CONCAT_BACKGROUND` | 84ms | 43.8% | 0 | 45ms | 64.9% | 0 |
| `CT:CONCAT_DISCOVERY` | 78ms | 30.2% | 0 | 26ms | 100.0% | 0 |
| `CT:RESAMPLE` | 7m | 377.3% | 2.1 GB | 5m 9s | 412.5% | 2 GB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 15.7s | 46.5% | 250.1 MB | 5.1s | 113.0% | 249.3 MB |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 1m 13s | 74.3% | 1.2 GB | 41.9s | 77.9% | 1.1 GB |
| `CT_POSTPROC:ASR_ROBUSTNESS:ASR_ROBUSTNESS_REPORT` | 14.3s | 54.0% | 323.5 MB | 15s | 50.9% | 325.7 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 2.6s | 58.2% | 105.2 MB | 2.8s | 52.4% | 132.4 MB |
| `CT_POSTPROC:CT_FILTER` | 1.8s | 67.4% | 121.2 MB | 1.8s | 70.2% | 122.4 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 1.4s | 78.0% | 24.7 MB | 3.9s | 30.3% | 118.3 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 5s | 70.5% | 129.1 MB | 6.2s | 60.0% | 126.5 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 1.4s | 79.6% | 24.8 MB | 2.2s | 52.6% | 96.9 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 1m 8s | 72.2% | 781.1 MB | 1m 19s | 65.5% | 773.6 MB |
| `FADE:FADE_JSON_TO_CSV_TOP` | 1.3s | 58.4% | 35.1 MB | 678ms | 133.3% | 30.4 MB |
| `FADE:FADE_REPORT_TOP` | 8.3s | 55.0% | 241.5 MB | 5.5s | 82.6% | 240.9 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 363ms | 138.7% | 41.2 MB | 486ms | 107.9% | 28.8 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 39ms | 55.2% | 0 | 47ms | 37.9% | 0 |
| `CT:CONCAT_RESAMPLE` | 37ms | 107.1% | 0 | 83ms | 45.8% | 4.4 MB |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 1.3s | 101.6% | 4 MB | 1.2s | 94.8% | 3.9 MB |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 1m 14s | 101.0% | 1.8 GB | 1m 35s | 100.5% | 1.7 GB |
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 5m 11s | 22.1% | 3 GB | 6m 1s | 19.3% | 2.4 GB |
| `CAAS_PERMULATION:CAAS_PERMS_AGGREGATE` | 753ms | 195.2% | 3.9 MB | 2.1s | 75.2% | 81.5 MB |
| `SCORING:SCORING_COMPUTE` | 1.7s | 153.6% | 142.3 MB | 3.9s | 80.5% | 139.9 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 4.5s | 113.9% | 254.1 MB | 9.7s | 67.2% | 376.8 MB |
| `SCORING:SCORING_REPORT` | 12.1s | 104.8% | 462.5 MB | 25.5s | 78.2% | 461.5 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 4.7s | 116.1% | 374.2 MB | 11.6s | 67.5% | 394.1 MB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run | |  | 7s | 105.1% | 431.4 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | |  | 1m 5s | 99.5% | 1.9 GB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run | |  | 51.9s | 99.6% | 2 GB |

`peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
