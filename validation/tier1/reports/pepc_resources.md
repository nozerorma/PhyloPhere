# Tier 1 PEPC: Resource & Runtime Report

Two PhyloPhere runs on the shared PEPC fixture (`validation/tier1/input/pepc/`), one per trait definition, identical except for the trait and the pruning it implies, plus the trait-independent OC/FUBAR run used as the site-level baseline in `pepc_results.md`.

| run | trait | tips (in-group + outgroup) | C4 / C3 | results dir |
|---|---|---|---|---|
| genotypic | `c4` | 77 + 1 | 23 / 54 | `output/pepc/results/c4_complete/` |
| phenotypic | `c4_phenotypic` | 70 + 1 (7 *Eleocharis* pruned via `prune/eleocharis_intermediate.txt`) | 20 / 50 | `output/pepc/results/c4_phenotypic_complete/` |

## Run configuration

`output/pepc/work/<run>_nxf_run/params.json`, compared key by key. Besides the run-scoped path (`outdir`) the two runs differ only in the trait (`traitname` / `secondary_trait` swapped) and the pruning parameters (`prune_data`, `prune_list`, `prune_list_secondary`: off vs the *Eleocharis* list).

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
| `ct_core_batch_size` | 20 |
| `caas_evidence_top_n` | 30 |

## System

| | |
|---|---|
| CPU | 13th Gen Intel Core i5-1335U, 10 cores (2P+8E) / 12 threads, 1 socket, 400 to 4600 MHz |
| RAM | 14 GiB |
| Swap | 21 GiB |
| Disk (`/home`) | 232 GB total, 51 GB available |
| Nextflow | 25.10.4, `-profile local`, env `phylophere` |

RAM, swap and disk are a snapshot at report time. The two runs were the only pipeline runs on the machine, next to the desktop session: `MemAvailable` was 5.7 GB at launch and between 5.5 and 8.8 GB during the runs, with about 16 GB of swap in use and no swap activity in `vmstat`. Other load was not recorded.

## PhyloPhere runs: totals

Both runs executed every task (`cachedCount=0`), sequentially.

| | genotypic | phenotypic |
|---|---|---|
| Nextflow run name | `c4_complete_2fff0887` | `c4_phenotypic_complete_72786aa9` |
| Launch → completion (`.nextflow.log`) | 2026-10-06 00:26:02 → 00:38:08 (12m 06s) | 2026-10-06 00:38:11 → 00:47:32 (9m 21s) |
| Processes succeeded / failed | 34 / 0 | 34 / 0 |
| Nextflow `succeedDuration` (cumulative task time) | 1h 19m 47s | 1h 3m 22s |
| `peakRunning` / `peakCpus` | 5 / 18 | 5 / 18 |
| `peakMemory` (requested, not resident) | 44 GB | 44 GB |
| Largest single-task `peak_rss` | 3.2 GB (`CAAS_CORE:CAAS_CORE_BATCHED`) | 3.1 GB (`CAAS_CORE:CAAS_CORE_BATCHED`) |
| Permulation cycles × hypotheses | 1000 × 100 | 1000 × 100 |

`peakMemory` is the sum of *requested* memory across concurrently running tasks and exceeds physical RAM; resident usage is bounded by the per-task `peak_rss` column.

## Structural difference between the two DAGs

| process | genotypic | phenotypic | reason |
|---|---|---|---|
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run | 29.8s, 1.9 GB | `prune_data=true` |
| `CONTRAST_SELECTION:DATASET_EXPLORATION`, `:PHENOTYPE_EXPLORATION` | not run | run | pruning path |
| `CONTRAST_SELECTION:REPORTING:{DATASET_EXPLORATION, PHENOTYPE_EXPLORATION, NAME_CURATION:TREE_CLEANUP}` | run | not run | non-pruning path |

Process count is 34 in both: the pruning and non-pruning paths contribute three processes each.

## Dominant cost centres

| process | genotypic realtime | %cpu | phenotypic realtime | %cpu |
|---|---|---|---|---|
| `CAAS_CORE:CAAS_CORE_BATCHED` | 6m 33s | 31.5 % | 4m 26s | 40.8 % |
| `CT:RESAMPLE` | 3m 50s | 609.5 % | 2m 45s | 655.5 % |
| `FADE:FADE_BATCHED` | 2m 9s | 98.6 % | 2m 1s | 99.8 % |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 30.9s | 96.0 % | not run |  |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 22.9s | 101.2 % | 21.8s | 101.1 % |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 19.5s | 100.8 % | 19.3s | 101.3 % |

- **`CAAS_CORE_BATCHED` holds the null of the single gene**: the perm-replay kernel and set-up (58 s genotypic, 50 s phenotypic, until the replay of the labelings starts), pass A over 100 000 labelings on 4 workers (330 s and 212 s; 161 768 and 127 521 pooled rows over 1000 cycles) and the scoring of the real labeling (a few seconds). Utilisation of this task was not profiled: its `%cpu` (31.5 % and 40.8 %) is below one core although pass A runs on 4 workers.
- **The ASR of PEPC is computed once**, inside the genotypic task (empty cache: one `codeml` run, the other workers wait for it); the phenotypic run reads it from the same cache (no `codeml` run).
- **`CT:RESAMPLE` breakdown** (from its `.command.log`):

| phase | genotypic | phenotypic |
|---|---|---|
| pool harvest | 32 s (1012 draws) | 17 s (1002 draws) |
| design matching: FOP harvest of every candidate | 1m 26s + 20 s top-up (1147 candidates, 89.0 % reach 100 hypotheses) | 1m 06s + 5 s top-up (1066 candidates, 94.4 % reach 100 hypotheses) |
| FOP mirror | 1m 30s (1000 cycles × 100 hypotheses, 8 workers) | 1m 14s (1000 cycles × 100 hypotheses, 8 workers) |

Design matching harvests each candidate cycle's FOP hypotheses once to count them and again in the mirror; together they are the largest share of `RESAMPLE` in both runs. The genotypic run harvests more candidates (lower match rate) with 4 rather than 3 pairs per hypothesis; the two contributions were not separated.

## Per-process trace, side by side

Sources: `results/<run>/pipeline_info/execution_trace.txt`. Rows ordered by genotypic submit time; phenotypic-only processes last.

| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |
|---|---|---|---|---|---|---|
| `CONTRAST_SELECTION:NAME_CURATION:TREE_CLEANUP` | 131ms | 94.0% | 24.1 MB | 137ms | 88.1% | 17.7 MB |
| `CONTRAST_SELECTION:REPORTING:NAME_CURATION:TREE_CLEANUP` | 137ms | 91.2% | 25.6 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION` | 3.7s | 117.9% | 431.2 MB | not run |  |  |
| `CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION` | 30.9s | 96.0% | 2 GB | not run |  |  |
| `CONTRAST_SELECTION:CI_COMPOSITION_REPORT` | 8.2s | 101.1% | 901.8 MB | 7.9s | 101.9% | 895.3 MB |
| `SELECTION_PREP:EXTRACT_EXTREME_SPECIES` | 69ms | 96.0% | 16.8 MB | 34ms | 90.6% | 3.8 MB |
| `CONTRAST_SELECTION:CONTRAST_ALGORITHM` | 19.5s | 100.8% | 1.2 GB | 19.3s | 101.3% | 1.2 GB |
| `SELECTION_PREP:PREP_ALIGNMENTS_BATCHED` | 203ms | 98.4% | 10.2 MB | 152ms | 104.6% | 8 MB |
| `CONTRAST_SELECTION:CHECK_MIN_CONTRASTS` | 45ms | 110.8% | 0 | 15ms | 128.6% | 0 |
| `CT:RESAMPLE` | 3m 50s | 609.5% | 2.5 GB | 2m 45s | 655.5% | 2.6 GB |
| `FADE:ANNOTATE_TREE_FG_BATCHED` | 189ms | 81.3% | 7.9 MB | 223ms | 69.1% | 7.9 MB |
| `FADE:FADE_BATCHED` | 2m 9s | 98.6% | 474.4 MB | 2m 1s | 99.8% | 489.8 MB |
| `FADE:FADE_JSON_TO_CSV_TOP` | 837ms | 121.8% | 25.8 MB | 710ms | 147.2% | 31.8 MB |
| `FADE:FADE_REPORT_TOP` | 5.6s | 92.3% | 242 MB | 5.1s | 87.0% | 348 MB |
| `FADE:FADE_GENE_LISTS_TOP` | 499ms | 134.0% | 28 MB | 227ms | 318.6% | 55.1 MB |
| `CT:CONCAT_RESAMPLE` | 22ms | 133.3% | 0 | 44ms | 109.1% | 4.4 MB |
| `CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS` | 758ms | 116.2% | 8.2 MB | 765ms | 104.9% | 8.1 MB |
| `CAAS_CORE:CAAS_CORE_BATCHED` | 6m 33s | 31.5% | 3.2 GB | 4m 26s | 40.8% | 3.1 GB |
| `CAAS_CORE_OBSERVED` | 254ms | 488.5% | 3.9 MB | 213ms | 604.5% | 3.9 MB |
| `CT_POSTPROC:ASR_ROBUSTNESS:ASR_ROBUSTNESS_REPORT` | 4.8s | 98.8% | 321.9 MB | 4.7s | 97.3% | 323.4 MB |
| `CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT` | 622ms | 214.0% | 27.6 MB | 581ms | 231.9% | 27 MB |
| `CT_META_CAAS:CAAS_META_CAAS_REPORT` | 3.1s | 99.7% | 236.5 MB | 2.8s | 94.0% | 235.9 MB |
| `CT_POSTPROC:CT_FILTER` | 548ms | 250.8% | 27.8 MB | 677ms | 224.7% | 24.2 MB |
| `CT_POSTPROC:CT_FILTER_SUMMARY` | 750ms | 197.6% | 23.7 MB | 577ms | 208.3% | 25.6 MB |
| `CT_POSTPROC:CAAS_FILTER_GENES` | 528ms | 228.7% | 29.4 MB | 731ms | 181.5% | 29.2 MB |
| `CT_POSTPROC:CAAS_BACKGROUND_CLEANUP` | 756ms | 182.5% | 25.5 MB | 560ms | 229.4% | 25 MB |
| `CT_POSTPROC:CT_POSTPROC_REPORT` | 22.9s | 101.2% | 774.5 MB | 21.8s | 101.1% | 774.4 MB |
| `ENRICHMENT:PUBLISH_UNIVERSES` | 14ms | 109.1% | 0 | 15ms | 72.7% | 0 |
| `CAAS_CORE_MERGE` | 3.8s | 100.9% | 209.6 MB | 3.5s | 100.9% | 192.6 MB |
| `SCORING:SCORING_COMPUTE` | 1.2s | 110.3% | 142.4 MB | 1.2s | 110.7% | 138.5 MB |
| `CAAS_SIGNIFICANCE_REPORT` | 4.8s | 101.2% | 375.2 MB | 4.1s | 103.4% | 373.7 MB |
| `CAAS_EVIDENCE` | 4.1s | 25.3% | 651.4 MB | 3.1s | 35.8% | 661.3 MB |
| `SCORING:SCORING_REPORT` | 9.8s | 101.3% | 469.7 MB | 10.4s | 100.5% | 461.8 MB |
| `ENRICHMENT:SCORING_COMPARE_REPORT` | 4.9s | 111.8% | 370.7 MB | 5s | 115.6% | 362.7 MB |
| `CONTRAST_SELECTION:DATASET_PRUNE` | not run |  |  | 29.8s | 100.5% | 1.9 GB |
| `CONTRAST_SELECTION:DATASET_EXPLORATION` | not run |  |  | 3.5s | 124.7% | 429.3 MB |
| `CONTRAST_SELECTION:PHENOTYPE_EXPLORATION` | not run |  |  | 29.1s | 97.3% | 1.8 GB |

`peak_vmem` values near 1 TB on Rmd-render steps in the raw traces are shared-library mmap accounting, not resident memory.

## OC / FUBAR baseline run

Trait-independent (site-level dN/dS over the whole tree), so it was run once, on the 78-tip fixture. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`), MEME disabled (`--meme_min_sites 999999`). Trace: `input/pepc/oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED.

| | |
|---|---|
| Wall clock (`time`) | 57.0 s real, 1m 34.3s user, 10.6 s sys |
| `HYPHY_BATCH` (FUBAR, 970 codons / 78 sequences) | 34.4 s, 62 MB peak RSS |
| `RENDER_TRANSLATION_REPORT` + `RENDER_PSEL_REPORT` | 14.4 s + 11.8 s, ≤ 496 MB peak RSS |
| All other steps | < 1.5 s each |
