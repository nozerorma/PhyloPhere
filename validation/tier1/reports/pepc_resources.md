# Tier 1 PEPC: Resource & Runtime Report

Two PhyloPhere runs on the shared PEPC fixture (`validation/tier1/input/pepc/`), one per trait definition, identical except for the trait and the pruning it implies, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` |

## Run configuration

The settings are those of `gui/templates/tier1_pepc_c4.json`, from which the two runs were generated; they differ only in the trait (`traitname` / `secondary_trait` swapped) and the pruning parameters (`prune_data`, `prune_list`, `prune_list_secondary`: off vs the *Eleocharis* list). The runs used `cleanup_policy = remove_work`, so the work directories and the `params.json` files were removed after each run; the trace, the Nextflow log and the results remain.

| setting | value |
|---|---|
| `seed` | 1998 |
| `caas_full_perms` (null cycles) | 1000 |
| `max_fop` / `min_contrasts` | 100 / 3 |
| `filter_minlen` / `filter_maxcaas` | 3 / 0.7 |
| `ct_disambig_asr_model` | lg |
| `fade_direction` / `fade_background_scope` / `fade_internal_nodes` | top / all / all_descendants |
| `fade_bf_threshold` | 100 |
| `fcs_enabled` | false |
| `ct_core_batch_size` | 20 |
| `caas_evidence_top_n` | 10 |
| `cleanup_policy` | remove_work |

Launch and end (`pipeline_info/nextflow.log`): 2026-10-14 16:19:22 to 16:29:04 (genotypic) and 16:29:08 to 16:37:02 (phenotypic). Code: branch `scoring_v2`, HEAD `f47de98` when the runs started.

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Disk (`/home`) | 232 GB total |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` |

The runs were the only pipeline runs started from this project, next to the desktop session and other applications; memory availability and swap activity during the runs were not recorded.

## PhyloPhere runs: totals

Both runs executed every task (`cachedCount=0`), sequentially.

| | genotypic | phenotypic |
|---|---|---|
| Nextflow run name | `c4_complete_92e7ffaa` | `c4_phenotypic_complete_81223265` |
| Launch → completion (`.nextflow.log`) | 2026-10-14 16:19:22 → 16:29:04 (9m 42s) | 2026-10-14 16:29:08 → 16:37:02 (7m 54s) |
| Processes succeeded / failed | 33 / 0 | 34 / 0 |
| Cached tasks | 0 | 0 |
| Nextflow `succeedDuration` (cumulative task time) | 1h 1m 45s | 52m 43s |
| `peakRunning` / `peakCpus` | 5 / 18 | 6 / 20 |
| `peakMemory` (requested, not resident) | 44 GB | 48 GB |
| Largest single-task `peak_rss` | 8.5 GB (`CAAS_CORE:CAAS_CORE_BATCHED`) | 4 GB (`CAAS_CORE:CAAS_CORE_BATCHED`) |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; resident usage is bounded by the per-task `peak_rss` column. The largest resident task, `CAAS_CORE_BATCHED` of the genotypic run (8.5 GB), is above half of the RAM of the machine.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | run | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION}` | run | not run | non-pruning path |

The process count is 33 (genotypic) and 34 (phenotypic): the pruning path contributes three processes and the non-pruning path two.

## Dominant cost centres

| process | genotypic realtime | %cpu | phenotypic realtime | %cpu |
|---|---|---|---|---|
| `CT:RESAMPLE` | 4m | 540.5 % | 3m 20s | 576.7 % |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 3m 37s | 33.2 % | 1m 55s | 55.3 % |
| `FADE:FADE_BATCHED` | 1m 8s | 101.1 % | 1m 14s | 97.6 % |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 33s | 93.9 % | not run |  |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 29.2s | 101.3 % | 29.6s | 101.5 % |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 26s | 101.4 % | 26.8s | 101.1 % |

- **`CAAS_CORE_BATCHED` holds the null of the single gene**: the replay of 100 000 labelings and of the 100 labelings of the real design (cycle b_0), and the scoring of the real labeling. Its `%cpu` (33 % and 55 %) is below one core although the replay runs on 4 workers. The internal timings of the replay and of `CT:RESAMPLE` come from the task logs, which were removed with the work directories.
- **The ASR of PEPC is computed once**, inside the genotypic task (empty cache: one `codeml` run, the other workers wait for it); the phenotypic run reads it from the same cache.

## Per-process trace, side by side

Sources: `results/<run>/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `NAME_CURATION:TREE_CLEANUP` | 588ms | 37.4% | 16.9 MB | 182ms | 89.1% | 23.5 MB |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 4.8s | 100.8% | 504.9 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 33s | 93.9% | 1.9 GB | not run |  |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 8.9s | 99.7% | 874.8 MB | 9.6s | 102.3% | 875.7 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 63ms | 84.7% | 18.5 MB | 44ms | 97.3% | 3.9 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 29.2s | 101.3% | 1.6 GB | 29.6s | 101.5% | 1.6 GB |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 235ms | 86.9% | 8 MB | 166ms | 109.1% | 7.9 MB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 19ms | 102.9% | 0 | 27ms | 106.7% | 0 |
| `CT:RESAMPLE` | 4m | 540.5% | 2.6 GB | 3m 20s | 576.7% | 2.6 GB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 150ms | 96.0% | 7.9 MB | 236ms | 106.1% | 7.9 MB |
| `FADE:FADE_BATCHED` | 1m 8s | 101.1% | 102.2 MB | 1m 14s | 97.6% | 102.8 MB |
| `FADE:FADE_TABLES_TOP` | 48ms | 101.4% | 3.8 MB | 115ms | 80.4% | 3.9 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 268ms | 344.4% | 27.8 MB | 467ms | 126.7% | 54.5 MB |
| `FADE:FADE_REPORT_TOP` | 2.3s | 93.9% | 241.7 MB | 4.3s | 84.7% | 241.7 MB |
| `CT:CONCAT_RESAMPLE` | 24ms | 126.3% | 0 | 28ms | 78.3% | 0 |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 833ms | 116.6% | 8.2 MB | 894ms | 119.2% | 8.1 MB |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 3m 37s | 33.2% | 8.5 GB | 1m 55s | 55.3% | 4 GB |
| `CAAS_CORE_OBSERVED` | 300ms | 430.1% | 3.9 MB | 381ms | 315.4% | 4 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 636ms | 201.3% | 31.2 MB | 848ms | 174.8% | 24.3 MB |
| `CT_POSTPROC:ASR_DIAGNOSTICS:ASR_DIAGNOSTICS_REPORT` | 9.8s | 96.3% | 373 MB | 9.4s | 95.2% | 371.5 MB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 3.9s | 84.5% | 281.2 MB | 4.8s | 83.6% | 340.5 MB |
| `CT_POSTPROC:CT_FILTER` | 499ms | 249.1% | 25.5 MB | 665ms | 209.9% | 27.6 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 641ms | 202.8% | 20.9 MB | 811ms | 165.5% | 26.1 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 662ms | 184.3% | 28.2 MB | 793ms | 135.4% | 30.4 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 589ms | 217.9% | 29.5 MB | 713ms | 174.8% | 29.9 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 26s | 101.4% | 766.9 MB | 26.8s | 101.1% | 765.8 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 12ms | 92.3% | 0 | 12ms | 92.3% | 0 |
| `CAAS_CORE_MERGE` | 4.5s | 100.6% | 205.7 MB | 4s | 100.4% | 183.7 MB |
| `SCORING:SCORING_COMPUTE` | 2s | 106.1% | 139 MB | 1.7s | 108.3% | 136.3 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 6s | 100.6% | 385.7 MB | 5.3s | 101.2% | 255.6 MB |
| `CAAS_EVIDENCE` | 2.7s | 39.6% | 626.9 MB | 2.3s | 45.1% | 276.8 MB |
| `SCORING:SCORING_REPORT` | 16.6s | 100.0% | 483.7 MB | 15.7s | 99.9% | 365.3 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 6.2s | 111.0% | 362.8 MB | 6.1s | 110.8% | 364 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run |  |  | 34.7s | 100.4% | 1.9 GB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run |  |  | 3.9s | 123.5% | 425.3 MB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run |  |  | 33s | 98.8% | 1.8 GB |

`peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
