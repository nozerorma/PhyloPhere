# Tier 1 PEPC: Resource & Runtime Report

Two PhyloPhere runs on the shared PEPC fixture (`validation/tier1/input/pepc/`), identical except for trait definition, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` |

The two `params.json` files (`output/pepc/work/<run>_nxf_run/params.json`) differ only in `outdir`, `tax_id`, `gmt_dir` (all run-scoped paths), `traitname`/`secondary_trait` (swapped), and `prune_data`/`prune_list`/`prune_list_secondary` (off vs. the Eleocharis list). The `resources.override.config` files are byte-identical.

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Disk (`/home`) | 232 GB total, 56 GB available |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` |

RAM/swap/disk values are a snapshot taken at report time, not during the runs; memory pressure during the runs was not recorded.

## PhyloPhere runs: totals

Both runs launched with `-resume` against empty caches (`cachedCount=0` in both logs), so every task executed. They ran sequentially.

| | genotypic (`c4_complete`) | phenotypic (`c4_phenotypic_complete`) |
|---|---|---|
| Nextflow run name | `c4_complete_1d1df5c2` | `c4_phenotypic_complete_7a0d3208` |
| Launch → completion (`.nextflow.log`) | 2026-09-27 13:59:21 → 14:09:19 (9m 58s) | 2026-09-27 14:10:24 → 14:18:52 (8m 28s) |
| Processes succeeded / failed | 37 / 0 | 37 / 0 |
| Nextflow `succeedDuration` (cumulative task time) | 1h 20m 13s | 1h 7m 46s |
| `peakRunning` / `peakCpus` | 5 / 26 | 5 / 26 |
| `peakMemory` (requested, not resident) | 28 GB | 28 GB |
| Largest single-task `peak_rss` | 2.9 GB (`CAAS_PERMS_DISAMBIGUATE`) | 2.5 GB (`CAAS_PERMS_DISAMBIGUATE`) |
| Permulation cycles in `perm_pos_cycle_caas.tsv.gz` | 999 of 1000 | 1000 of 1000 |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; actual resident usage is bounded by the per-task `peak_rss` column. One of the 1000 replayed cycles is absent from the genotypic null file (cause not traced; most plausibly a cycle that re-detected no position). `scoring_compute.R` sets `N` to the number of distinct cycles present, so `p.emp = (k+1)/1000` there and `(k+1)/1001` in the phenotypic run.

Module set per `gui/templates/tier1_pepc_c4.json`: CAAS, CT_DISAMBIGUATION, FADE, SCORING, ENRICHMENT. Within ENRICHMENT, `fcs_enabled=false` and `posenrich=false`, so only `PUBLISH_UNIVERSES` and `SCORING_COMPARE_REPORT` execute; no FCS task appears in either trace.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | 35.2 s, 1.9 GB | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION, NAME_CURATION:TREE_CLEANUP}` | run | not run | non-pruning path |

Process count is 37 in both. The SCORING task tags read `scoring_compute|c4` / `scoring_report|c4` in both traces, including the phenotypic run; the outputs themselves are correctly trait-scoped (`11.Scoring_report_c4_phenotypic.html`).

## Dominant cost centres

| process | genotypic realtime | phenotypic realtime | notes |
|---|---|---|---|
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 3m 55s | 2m 49s | 24 to 32 % CPU; low utilisation, bottleneck not profiled |
| `CT:RESAMPLE` | 2m 41s | 2m 1s | only clearly multi-threaded step (411 % / 447 % CPU) |
| `FADE:FADE_BATCHED` | 2m 8s | 1m 55s | single gene, top side only |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 54.3 s | 44.6 s | |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 52.7 s | 46.8 s | Rmd render |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 40.8 s | 39.3 s | |

The phenotypic run is faster on every permulation-bound step. It has 7 fewer tips and 3 rather than 4 contrast pairs per hypothesis (see `pepc_genotypic_vs_phenotypic.md`), both of which shrink the per-cycle replay; which of the two dominates the saving was not isolated. `peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## Per-process trace, side by side

Sources: `results/c4_complete/pipeline_info/execution_trace.txt`, `results/c4_phenotypic_complete/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `CONTRAST_SELECTION:NAME_CURATION:TREE_CLEANUP` | 160ms | 89.6% | 26.4 MB | 121ms | 101.3% | 21.1 MB |
| `CONTRAST_SELECTION:REPORTING:NAME_CURATION:TREE_CLEANUP` | 151ms | 101.6% | 26.7 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 5s | 110.9% | 432.8 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 42.2s | 97.8% | 2.1 GB | not run |  |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 11.5s | 105.3% | 902.6 MB | 9.3s | 111.6% | 902.7 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 107ms | 78.8% | 15.8 MB | 68ms | 102.1% | 18.6 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 28.2s | 104.1% | 1.2 GB | 22.5s | 105.4% | 1.2 GB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 26ms | 83.7% | 0 | 18ms | 70.6% | 0 |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 259ms | 91.1% | 7.9 MB | 170ms | 106.0% | 7.9 MB |
| `CT:DISCOVERY_BATCHED` | 1.8s | 95.5% | 67.5 MB | 1.2s | 102.7% | 67.7 MB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 144ms | 113.3% | 40.1 MB | 178ms | 103.6% | 8 MB |
| `FADE:FADE_BATCHED` | 2m 8s | 93.9% | 471 MB | 1m 55s | 94.3% | 457.6 MB |
| `CT:CONCAT_DISCOVERY` | 217ms | 30.2% | 3.9 MB | 21ms | 100.0% | 0 |
| `CT:CONCAT_BACKGROUND` | 183ms | 41.4% | 4 MB | 18ms | 75.0% | 0 |
| `CT:RESAMPLE` | 2m 41s | 411.2% | 2.1 GB | 2m 1s | 447.4% | 2 GB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 4.9s | 112.0% | 248.7 MB | 4.8s | 101.6% | 312.1 MB |
| `CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN` | 40.8s | 88.6% | 1.2 GB | 39.3s | 84.4% | 1.1 GB |
| `CT_POSTPROC:ASR_ROBUSTNESS:ASR_ROBUSTNESS_REPORT` | 6.9s | 89.5% | 362.9 MB | 11.4s | 61.7% | 433 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 1.2s | 98.4% | 27.6 MB | 2.3s | 66.7% | 92.6 MB |
| `CT_POSTPROC:CT_FILTER` | 790ms | 192.5% | 23.8 MB | 1.3s | 114.5% | 120.9 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 581ms | 188.8% | 29.4 MB | 1.3s | 77.9% | 24.9 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 1.9s | 122.1% | 129.8 MB | 4.3s | 86.2% | 125.4 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 890ms | 156.2% | 25.7 MB | 3.8s | 36.8% | 118.4 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 52.7s | 88.9% | 776.9 MB | 46.8s | 90.8% | 772.2 MB |
| `FADE:FADE_JSON_TO_CSV_TOP` | 944ms | 75.2% | 19.6 MB | 997ms | 59.7% | 33.4 MB |
| `FADE:FADE_REPORT_TOP` | 5.5s | 83.5% | 241 MB | 5s | 85.9% | 240.5 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 242ms | 318.7% | 57.6 MB | 256ms | 196.9% | 51 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 16ms | 102.9% | 0 | 29ms | 59.0% | 0 |
| `CT:CONCAT_RESAMPLE` | 31ms | 94.1% | 4.4 MB | 25ms | 85.7% | 0 |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 667ms | 103.6% | 3.9 MB | 625ms | 104.1% | 8 MB |
| `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` | 54.3s | 101.7% | 1.7 GB | 44.6s | 101.7% | 1.7 GB |
| `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` | 3m 55s | 24.4% | 2.9 GB | 2m 49s | 31.9% | 2.5 GB |
| `CAAS_PERMULATION:CAAS_PERMS_AGGREGATE` | 595ms | 262.9% | 3.9 MB | 609ms | 254.6% | 3.9 MB |
| `SCORING:SCORING_COMPUTE` | 1.5s | 162.0% | 140.4 MB | 1.3s | 185.8% | 138.1 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 3.8s | 116.6% | 375 MB | 3.5s | 120.5% | 412.3 MB |
| `SCORING:SCORING_REPORT` | 9.4s | 105.0% | 461.1 MB | 9.9s | 109.0% | 590.9 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 3.9s | 119.8% | 347 MB | 3.9s | 123.2% | 347.1 MB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run |  |  | 4.3s | 118.7% | 430.4 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run |  |  | 35.2s | 100.8% | 1.9 GB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run |  |  | 33.3s | 102.7% | 1.9 GB |

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip genotypic fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
