# Tier 1 PEPC: Resource & Runtime Report

Four PhyloPhere runs (two traits, lenient and strict discovery) on the shared PEPC fixture (`validation/tier1/input/pepc/`), one per trait definition, identical except for the trait and the pruning it implies, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` (lenient), `output/pepc_nocons/results/c4_complete/` (strict) |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` (lenient), `output/pepc_nocons/results/c4_phenotypic_complete/` (strict) |

## Run configuration

Two modes, each with the two traits, generated from `gui/templates/tier1_pepc_c4.json`: **lenient** (`output/pepc/`) and **strict** (`output/pepc_nocons/`). The modes share every setting except the divergence requirement of the discovery: lenient uses `min_divergent_fraction` 0.5 (up to 2 of 4 pairs or 1 of 3 may keep the ancestral state), strict requires divergence in every pair (`n_conserved_pairs` is 0 in all positions; the value of `min_divergent_fraction` of the strict runs is not recorded). Within a mode the traits differ only in `traitname` / `secondary_trait` and in the pruning parameters (`prune_data`, `prune_list`, `prune_list_secondary`: off vs the *Eleocharis* list). The runs used `cleanup_policy = remove_work`, so the work directories and the `params.json` files were removed after each run; the trace, the Nextflow log and the results remain.

| setting | value |
|---|---|
| input | 970 columns, 78 sequences; `c4`: 23 C₄ / 54 C₃ (77 tips); `c4_phenotypic`: 20 C₄ / 50 C₃ (70 tips, 7 *Eleocharis* pruned); ordinal binary trait |
| `seed` | 1998 |
| contrast selection | `min_contrasts` 3, `max_contrasts` 0 (until Dunn < 1), `max_fop` 100; BM selected by AIC in both traits |
| contrast hypotheses | 100 × 4 pairs (129 distinct pairs; 23 top and 25 bottom species) and 100 × 3 pairs (130 distinct pairs; 20 top and 29 bottom species) |
| discovery | patterns 1, 2, 3; `max_*_gaps_fraction` and `max_*_miss_fraction` 0.0, `miss_pair` true; scheme US only |
| ASR | LG, `ct_disambig_posterior_threshold` 0.1 |
| null | `caas_full_perms` 1000, `perm_strategy` auto (BM/OU by AICc), `perm_match_pss` true with tolerance 0.25 |
| filtering | cluster filter `filter_minlen` 3, `filter_maxcaas` 0.7, `remove_caas_clusters` true; `gene_filter_mode` none; `caas_perms_postproc` true (the null follows the filters) |
| FADE | direction top, background `all`, `all_descendants`, BF ≥ 100, Variational-Bayes, LG |
| other modules | scoring and enrichment on; accumulation, RER, VEP, FCS and POSENRICH off |
| `ct_core_batch_size` / `caas_evidence_top_n` | 20 / 10 |

Funnel of positions (discovered, discarded by the cluster filter, with a score above 0):

| run | discovered | discarded | passing the filter | with score > 0 |
|---|---|---|---|---|
| genotypic, lenient | 62 | 10 | 52 | 35 |
| genotypic, strict | 10 | 0 | 10 | 10 |
| phenotypic, lenient | 62 | 16 | 46 | 34 |
| phenotypic, strict | 18 | 0 | 18 | 17 |

The positions that pass the filter and have no score have `asr_path_score` 0 and side `none`: no domain changed independently (16 and 12 positions in lenient mode, 0 and 1 in strict mode).

Launch and end (`pipeline_info/nextflow.log`): lenient 2026-10-15 01:41:49 to 01:48:13 (genotypic) and 01:48:16 to 01:53:48 (phenotypic); strict 10:43:31 to 10:49:40 and 10:49:44 to 10:55:18. Code: branch `scoring_v2`, HEAD `9ae9bd2` when the runs started.

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` (8 CPUs and 16 GB allowed to the pipeline) |

The runs were the only pipeline runs started from this project, next to the desktop session and other applications; memory availability and swap activity during the runs were not recorded.

## PhyloPhere runs: totals

All runs executed every task (`cachedCount=0`), sequentially (genotypic, then phenotypic).

| | genotypic, lenient | genotypic, strict | phenotypic, lenient | phenotypic, strict |
|---|---|---|---|---|
| Nextflow run name | `c4_complete_eaa07c32` | `c4_complete_335d3a11` | `c4_phenotypic_complete_01e1555b` | `c4_phenotypic_complete_964b14e2` |
| Launch → completion (`.nextflow.log`) | 2026-10-15 01:41:49 → 01:48:13 (6m 24s) | 2026-10-15 10:43:31 → 10:49:40 (6m 09s) | 2026-10-15 01:48:16 → 01:53:48 (5m 32s) | 2026-10-15 10:49:44 → 10:55:18 (5m 34s) |
| Processes succeeded / failed | 33 / 0 | 33 / 0 | 34 / 0 | 34 / 0 |
| Cached tasks | 0 | 0 | 0 | 0 |
| Nextflow `succeedDuration` (cumulative task time) | 45m 34s | 43m 35s | 37m 13s | 38m 2s |
| `peakRunning` / `peakCpus` | 5 / 18 | 5 / 18 | 5 / 18 | 5 / 18 |
| `peakMemory` (requested, not resident) | 44 GB | 44 GB | 44 GB | 44 GB |
| Largest single-task `peak_rss` | 4.4 GB (`CAAS_CORE:CAAS_CORE_BATCHED`) | 2.7 GB (`CT:RESAMPLE`) | 2.6 GB (`CT:RESAMPLE`) | 2.5 GB (`CT:RESAMPLE`) |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; resident usage is bounded by the per-task `peak_rss` column.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | run | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION}` | run | not run | non-pruning path |

The process count is 33 (genotypic) and 34 (phenotypic) in both modes: the pruning path contributes three processes and the non-pruning path two.

## Dominant cost centres

Realtime and peak resident memory (`peak_rss`) of the processes that differ or dominate:

| process | geno, lenient | geno, strict | pheno, lenient | pheno, strict |
|---|---|---|---|---|
| `CT:RESAMPLE` | 3m 37s, 2.7 GB | 3m 48s, 2.7 GB | 2m 39s, 2.6 GB | 2m 43s, 2.5 GB |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 59s, 4.4 GB | 21.2s, 2.1 GB | 35.8s, 2.5 GB | 19.9s, 2 GB |
| `CAAS_CORE_MERGE` | 3s, 141.8 MB | 1.4s, 75.3 MB | 2.3s, 150.6 MB | 2.2s, 86 MB |
| `SCORING:SCORING_COMPUTE` | 1.3s, 130.5 MB | 1.3s, 121.6 MB | 1.1s, 128 MB | 1.1s, 20 MB |
| `SCORING:SCORING_REPORT` | 13.5s, 361.6 MB | 13.7s, 354.8 MB | 13s, 438.2 MB | 13.3s, 368.8 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 17.5s, 734.7 MB | 17.1s, 728.7 MB | 18.1s, 677.6 MB | 15.8s, 727.8 MB |
| `FADE:FADE_BATCHED` | 57.5s, 102.5 MB | 52s, 101.3 MB | 51.3s, 100.2 MB | 59.2s, 102 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 28.2s, 1.6 GB | 29.5s, 1.6 GB | 28.3s, 1.6 GB | 29.8s, 1.6 GB |

- **`CT:RESAMPLE`** (the design matching and the null labelings, about 6 cores) takes 2.7 to 3.8 minutes and 2.5 to 2.7 GB in every run, and does not depend on the mode.
- **`CAAS_CORE_BATCHED`** holds the replay of the null of the single gene and depends on the mode: 59 s and 4.4 GB (genotypic) and 36 s and 2.5 GB (phenotypic) in lenient mode, 21 s and 2.1 GB and 20 s and 2.0 GB in strict mode. The null of strict mode has 513 and 3 200 detections over 1000 cycles, against 24 449 and 18 409 in lenient mode, and the null table is 7 KB and 45 KB against 405 KB and 302 KB (compressed).
- **The ASR of PEPC is computed once**, inside the genotypic task of each mode (empty cache: one `codeml` run, the other workers wait for it); the phenotypic run of the mode reads it from the same cache. The internal timings of `CT:RESAMPLE` and of the replay come from task logs that were removed with the work directories.
- **Result directories** take 196 MB and 127 MB (genotypic and phenotypic, lenient) and 112 MB and 122 MB (strict).

## Per-process trace, side by side

Sources: `results/<run>/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

### Lenient

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `NAME_CURATION:TREE_CLEANUP` | 143ms | 72.9% | 18.5 MB | 137ms | 95.5% | 25.2 MB |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 4.8s | 104.0% | 496.8 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 31.4s | 96.8% | 1.9 GB | not run |  |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 8.1s | 101.9% | 874.8 MB | 8.2s | 102.4% | 875.5 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 77ms | 83.2% | 16.6 MB | 56ms | 91.1% | 18.7 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 28.2s | 101.1% | 1.6 GB | 28.3s | 101.3% | 1.6 GB |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 221ms | 82.6% | 7.9 MB | 184ms | 89.1% | 7.8 MB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 17ms | 116.1% | 0 | 15ms | 133.3% | 0 |
| `CT:RESAMPLE` | 3m 37s | 565.8% | 2.7 GB | 2m 39s | 573.3% | 2.6 GB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 125ms | 105.7% | 7.9 MB | 145ms | 80.0% | 7.9 MB |
| `FADE:FADE_BATCHED` | 57.5s | 100.8% | 102.5 MB | 51.3s | 100.6% | 100.2 MB |
| `FADE:FADE_TABLES_TOP` | 114ms | 100.0% | 3.9 MB | 62ms | 96.6% | 3.8 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 452ms | 136.8% | 28 MB | 258ms | 214.5% | 41.4 MB |
| `FADE:FADE_REPORT_TOP` | 4.3s | 94.1% | 241.4 MB | 3s | 95.9% | 241.3 MB |
| `CT:CONCAT_RESAMPLE` | 22ms | 87.8% | 0 | 22ms | 133.3% | 0 |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 759ms | 117.5% | 8.1 MB | 650ms | 119.8% | 8.1 MB |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 59s | 40.6% | 4.4 GB | 35.8s | 51.1% | 2.5 GB |
| `CAAS_CORE_OBSERVED` | 222ms | 580.4% | 3.8 MB | 201ms | 624.0% | 3.9 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 553ms | 228.9% | 24.5 MB | 665ms | 187.1% | 27 MB |
| `CT_POSTPROC:ASR_DIAGNOSTICS:ASR_DIAGNOSTICS_REPORT` | 4.9s | 100.0% | 359 MB | 4.7s | 100.9% | 359.1 MB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 2.4s | 95.6% | 355.1 MB | 1.9s | 97.6% | 201.3 MB |
| `CT_POSTPROC:CT_FILTER` | 518ms | 226.3% | 24 MB | 496ms | 279.1% | 23.8 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 590ms | 229.9% | 17.3 MB | 509ms | 198.4% | 28.1 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 526ms | 246.8% | 15.5 MB | 479ms | 165.5% | 26 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 604ms | 233.5% | 25.4 MB | 419ms | 269.8% | 29.7 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 17.5s | 101.3% | 734.7 MB | 18.1s | 101.3% | 677.6 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 13ms | 138.5% | 0 | 11ms | 156.5% | 0 |
| `CAAS_CORE_MERGE` | 3s | 100.7% | 141.8 MB | 2.3s | 101.3% | 150.6 MB |
| `SCORING:SCORING_COMPUTE` | 1.3s | 112.4% | 130.5 MB | 1.1s | 112.2% | 128 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 2.9s | 103.9% | 244.7 MB | 2.8s | 104.8% | 254.8 MB |
| `CAAS_EVIDENCE` | 1.4s | 41.9% | 588.5 MB | 1.4s | 43.5% | 581.9 MB |
| `SCORING:SCORING_REPORT` | 13.5s | 101.6% | 361.6 MB | 13s | 101.1% | 438.2 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 4s | 116.2% | 355.8 MB | 3.6s | 116.9% | 357.4 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run |  |  | 31.1s | 102.0% | 2 GB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run |  |  | 3.7s | 123.0% | 424.8 MB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run |  |  | 30.3s | 97.6% | 1.8 GB |

### Strict

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `NAME_CURATION:TREE_CLEANUP` | 192ms | 90.4% | 18.4 MB | 119ms | 95.4% | 18.9 MB |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 3.9s | 122.5% | 425.2 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 33.8s | 97.4% | 1.9 GB | not run |  |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 8.9s | 100.3% | 874.5 MB | 8.5s | 102.6% | 875.5 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 98ms | 85.7% | 15.2 MB | 85ms | 88.1% | 17.5 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 29.5s | 100.3% | 1.6 GB | 29.8s | 101.6% | 1.6 GB |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 259ms | 83.3% | 7.8 MB | 111ms | 108.3% | 7.9 MB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 18ms | 112.5% | 0 | 14ms | 133.3% | 0 |
| `CT:RESAMPLE` | 3m 48s | 584.2% | 2.7 GB | 2m 43s | 562.3% | 2.5 GB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 173ms | 83.3% | 7.9 MB | 180ms | 106.2% | 7.9 MB |
| `FADE:FADE_BATCHED` | 52s | 101.0% | 101.3 MB | 59.2s | 100.2% | 102 MB |
| `FADE:FADE_TABLES_TOP` | 62ms | 111.6% | 18.6 MB | 112ms | 75.9% | 18.5 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 174ms | 383.6% | 57.5 MB | 411ms | 159.2% | 25.8 MB |
| `FADE:FADE_REPORT_TOP` | 2.5s | 94.3% | 355.5 MB | 4.4s | 84.3% | 241.3 MB |
| `CT:CONCAT_RESAMPLE` | 37ms | 113.2% | 0 | 31ms | 120.0% | 0 |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 799ms | 119.5% | 8.1 MB | 788ms | 114.4% | 3.9 MB |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 21.2s | 89.1% | 2.1 GB | 19.9s | 86.6% | 2 GB |
| `CAAS_CORE_OBSERVED` | 202ms | 623.8% | 3.9 MB | 217ms | 588.8% | 4 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 592ms | 196.1% | 31.2 MB | 681ms | 180.8% | 21.6 MB |
| `CT_POSTPROC:ASR_DIAGNOSTICS:ASR_DIAGNOSTICS_REPORT` | 4.9s | 98.2% | 325.2 MB | 4.9s | 99.4% | 325.2 MB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 2.4s | 95.1% | 362 MB | 2s | 99.6% | 197.7 MB |
| `CT_POSTPROC:CT_FILTER` | 468ms | 260.3% | 25 MB | 525ms | 246.8% | 25.1 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 555ms | 178.6% | 27.7 MB | 611ms | 199.7% | 21.8 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 468ms | 173.2% | 28.4 MB | 428ms | 178.4% | 27.3 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 629ms | 209.0% | 18.3 MB | 608ms | 230.5% | 21.3 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 17.1s | 100.8% | 728.7 MB | 15.8s | 101.6% | 727.8 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 9ms | 133.3% | 0 | 10ms | 128.6% | 0 |
| `CAAS_CORE_MERGE` | 1.4s | 100.5% | 75.3 MB | 2.2s | 100.7% | 86 MB |
| `SCORING:SCORING_COMPUTE` | 1.3s | 109.4% | 121.6 MB | 1.1s | 112.4% | 20 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 3.6s | 101.4% | 253.7 MB | 3.4s | 102.9% | 378.3 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 3s | 124.8% | 216.9 MB | 3.4s | 113.8% | 372.7 MB |
| `CAAS_EVIDENCE` | 1.2s | 48.0% | 261.1 MB | 1.6s | 48.9% | 255.9 MB |
| `SCORING:SCORING_REPORT` | 13.7s | 101.1% | 354.8 MB | 13.3s | 101.0% | 368.8 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run |  |  | 32.7s | 101.4% | 1.9 GB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run |  |  | 3.8s | 125.3% | 425.1 MB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run |  |  | 31.2s | 98.3% | 2 GB |

`peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
