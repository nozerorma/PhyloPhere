# Tier 1 PEPC: Resource & Runtime Report

Two PhyloPhere runs on the shared PEPC fixture (`validation/tier1/input/pepc/`), one per trait definition, identical except for the trait and the pruning it implies, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` |

## Run configuration

`output/pepc/work/<run>_nxf_run/params.json`, compared key by key. Besides run-scoped paths (`outdir`, `gmt_dir`) the two runs differ only in the trait (`traitname` / `secondary_trait` swapped) and the pruning parameters (`prune_data`, `prune_list`, `prune_list_secondary`: off vs the *Eleocharis* list).

| setting | value |
|---|---|
| `seed` | 1998 |
| `caas_full_perms` (null cycles) | 1000 |
| `max_fop` / `min_contrasts` | 100 / 3 |
| `filter_minlen` / `filter_maxcaas` | 3 / 0.7 |
| `scoring_p_emp_thr` | 0.05 |
| `tax_id` | `input/pepc/taxid.tsv` (fixture table; no NCBI lookup) |
| `fade_direction` / `fade_background_scope` / `fade_internal_nodes` | top / all / all_descendants |
| `fcs_enabled` | false |

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Disk (`/home`) | 232 GB total, 54 GB available |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` |

RAM, swap and disk are a snapshot at report time. The two runs were the only pipeline runs on the machine; load from other processes was not recorded.

## PhyloPhere runs: totals

Both runs executed every task (`cachedCount=0`), sequentially.

| | genotypic | phenotypic |
|---|---|---|
| Nextflow run name | `c4_complete_0431ed9d` | `c4_phenotypic_complete_4cd1fe06` |
| Launch → completion (`.nextflow.log`) | 2026-09-28 19:01:16 → 19:13:15 (11m 59s) | 2026-09-28 19:13:18 → 19:23:13 (9m 55s) |
| Processes succeeded / failed | 37 / 0 | 37 / 0 |
| Nextflow `succeedDuration` (cumulative task time) | 1h 26m 13s | 1h 13m 28s |
| `peakRunning` / `peakCpus` | 5 / 26 | 5 / 26 |
| `peakMemory` (requested, not resident) | 28 GB | 28 GB |
| Largest single-task `peak_rss` | 3.0 GB (`CAAS_PERMS_DISAMBIGUATE`) | 2.4 GB (`CAAS_PERMS_DISAMBIGUATE`) |
| Permulation cycles × hypotheses | 1000 × 100 | 1000 × 100 |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; resident usage is bounded by the per-task `peak_rss` column.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | 34.4 s, 1.9 GB | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION, NAME_CURATION:TREE_CLEANUP}` | run | not run | non-pruning path |

Process count is 37 in both: the pruning and non-pruning paths contribute three processes each.

## Dominant cost centres

| process | genotypic realtime | %cpu | phenotypic realtime | %cpu |
|---|---|---|---|---|
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 4m 28s | 24.9 % | 3m 1s | 30.5 % |
| `CT:RESAMPLE` | 4m 16s | 446.2 % | 3m 15s | 483.9 % |
| `FADE:FADE_BATCHED` | 1m 40s | 100.8 % | 1m 49s | 98.8 % |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 46.2 s | 101.5 % | 43.8 s | 101.8 % |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 38.7 s | 83.0 % | 25.3 s | 87.7 % |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 36.8 s | 100.7 % | 51 s | 97.6 % |

- **`CAAS_PERMS_DISAMBIGUATE` runs at 25 to 31 % CPU**: most of its wall time is not spent computing on the requested cores. Whether the gap is I/O or worker start-up was not profiled.
- **`CT:RESAMPLE` breakdown** (from its `.command.log`):

| phase | genotypic | phenotypic |
|---|---|---|
| pool harvest, 1000 Tier-1 cycles | 45 s (1012 draws) | 21 s (1002 draws) |
| design matching: FOP harvest of every candidate | 1m 33s + 15 s top-up (1147 candidates, 89.0 % reach 100 hypotheses) | 1m 17s + 5 s (1066 candidates, 94.4 %) |
| FOP mirror, 1000 cycles × 100 hypotheses, 6 workers | 1m 38s | 1m 27s |

Design matching harvests each candidate cycle's FOP hypotheses once to count them and again in the mirror; together they are the largest share of `RESAMPLE` in both runs. The genotypic run harvests more candidates (lower match rate) with 4 rather than 3 pairs per hypothesis; the two contributions were not separated.

## Per-process trace, side by side

Sources: `results/<run>/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `CONTRAST_SELECTION:NAME_CURATION:TREE_CLEANUP` | 245ms | 76.5% | 16.6 MB | 142ms | 85.7% | 17.9 MB |
| `CONTRAST_SELECTION:REPORTING:NAME_CURATION:TREE_CLEANUP` | 222ms | 87.6% | 17.7 MB | not run | |  |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 5.5s | 99.6% | 433.6 MB | not run | |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 36.9s | 100.4% | 2 GB | not run | |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 9.7s | 107.8% | 1 GB | 8.6s | 114.2% | 897.9 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 98ms | 74.4% | 15.3 MB | 45ms | 109.1% | 3.9 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 21.7s | 105.3% | 1.2 GB | 21.3s | 106.1% | 1.2 GB |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 245ms | 84.3% | 7.9 MB | 125ms | 96.9% | 7.8 MB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 25ms | 126.3% | 0 | 15ms | 80.0% | 0 |
| `CT:DISCOVERY_BATCHED` | 1.3s | 95.3% | 68.3 MB | 1.1s | 102.8% | 4 MB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 125ms | 113.2% | 8 MB | 110ms | 100.7% | 8 MB |
| `FADE:FADE_BATCHED` | 1m 40s | 100.8% | 499.8 MB | 1m 49s | 98.8% | 513.3 MB |
| `CT:CONCAT_DISCOVERY` | 24ms | 83.7% | 0 | 21ms | 105.9% | 0 |
| `CT:CONCAT_BACKGROUND` | 21ms | 123.1% | 0 | 24ms | 117.1% | 0 |
| `CT:RESAMPLE` | 4m 16s | 446.2% | 2.2 GB | 3m 15s | 483.9% | 2.1 GB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 5.4s | 116.8% | 249.6 MB | 3.5s | 118.5% | 361.6 MB |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 38.7s | 83.0% | 1.1 GB | 25.3s | 87.7% | 1.1 GB |
| `CT_POSTPROC:ASR_ROBUSTNESS:ASR_ROBUSTNESS_REPORT` | 6.6s | 97.5% | 323.2 MB | 8.9s | 87.5% | 324.9 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 1.2s | 127.4% | 17.7 MB | 1.7s | 93.9% | 124.4 MB |
| `CT_POSTPROC:CT_FILTER` | 734ms | 202.5% | 24.8 MB | 1.8s | 89.5% | 123.1 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 669ms | 192.8% | 28.4 MB | 1s | 104.1% | 27.3 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 1.5s | 142.3% | 139.2 MB | 2.4s | 109.9% | 132 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 792ms | 163.6% | 27 MB | 1.1s | 110.7% | 27.5 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 36.8s | 100.7% | 929.9 MB | 51s | 97.6% | 779.3 MB |
| `FADE:FADE_JSON_TO_CSV_TOP` | 691ms | 148.9% | 41.2 MB | 620ms | 131.2% | 28.7 MB |
| `FADE:FADE_REPORT_TOP` | 4.4s | 105.2% | 241.5 MB | 5.1s | 91.2% | 313.7 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 234ms | 307.1% | 45.4 MB | 324ms | 247.7% | 35.2 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 16ms | 72.7% | 0 | 17ms | 94.7% | 0 |
| `CT:CONCAT_RESAMPLE` | 48ms | 94.7% | 0 | 28ms | 106.7% | 0 |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 646ms | 103.8% | 8.1 MB | 566ms | 103.3% | 8.2 MB |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 46.2s | 101.5% | 1.8 GB | 43.8s | 101.8% | 1.7 GB |
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 4m 28s | 24.9% | 3 GB | 3m 1s | 30.5% | 2.4 GB |
| `CAAS_PERMULATION:CAAS_PERMS_AGGREGATE` | 570ms | 273.1% | 3.9 MB | 551ms | 262.8% | 3.9 MB |
| `SCORING:SCORING_COMPUTE` | 1.3s | 175.6% | 143 MB | 1.1s | 192.8% | 24.2 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 3.8s | 116.4% | 376.3 MB | 3.6s | 119.6% | 379.1 MB |
| `SCORING:SCORING_REPORT` | 9.8s | 107.6% | 473.6 MB | 8.7s | 107.9% | 465.2 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 4.9s | 119.3% | 371.5 MB | 3.9s | 120.7% | 360.9 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | |  | 34.4s | 101.4% | 1.9 GB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run | |  | 3.9s | 125.6% | 429.4 MB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run | |  | 32.5s | 102.8% | 1.8 GB |

`peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
