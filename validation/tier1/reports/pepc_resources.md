# Tier 1 PEPC — Resource & Runtime Report

Fixture: `validation/tier1/input/pepc/` (78 tips, 23 C4 / 55 C3). Both runs executed locally, sequentially (not concurrently), on the system below.

## System specs

| | |
|---|---|
| CPU | 13th Gen Intel(R) Core(TM) i5-1335U — 10 cores (2P+8E) / 12 threads, 1 socket |
| CPU max/min MHz | 4600 / 400 |
| RAM | 14 GiB total |
| Swap | 21 GiB total |
| Disk (repo volume) | 110 GB total, 24 GB available (79% used) |
| Nextflow | 25.10.4 (`phylophere` env); OC pipeline run under `bmge-tools` env |

Memory pressure during these runs: free RAM fluctuated between ~330 MiB and ~5 GiB across the two runs, with swap usage between 9–10 GiB throughout. Neither run was observed to fail or be OOM-killed.

## Run 1 — ortholog_characterizator (OC), FUBAR only

Cold start: `oc_run/` removed before launch, no `-resume`. Launcher: `validation/tier1/input/pepc/scripts/run_ortholog_characterizator.sh pepc` (env `bmge-tools`). MEME disabled (`--meme_min_sites 999999`).

**Total wall-clock** (`time nextflow run ...`):

| | |
|---|---|
| real | 0m57.046s |
| user | 1m34.260s |
| sys | 0m10.576s |

**Per-process trace** (`oc_run/nextflow_trace_pepc_20260922_110520.tsv`, 10 processes, all COMPLETED):

| name | status | duration | realtime | %cpu | peak_rss | peak_vmem | rchar | wchar |
|---|---|---|---|---|---|---|---|---|
| TRANSLATION:TRIM_TRANSLATE:BMGE_BATCH (1) | COMPLETED | 919ms | 863ms | 174.4% | 7.1 MB | 26.1 MB | 1.1 MB | 4.1 MB |
| TRANSLATION:TRIM_TRANSLATE:BMGE_MAP_BATCH (1) | COMPLETED | 406ms | 365ms | 59.2% | 7 MB | 26.1 MB | 5.3 MB | 25.8 KB |
| POSITIVE_SELECTION:HYPHY:DEDUP_IDENTICAL_BATCH (1) | COMPLETED | 698ms | 654ms | 199.2% | 17.3 MB | 1 GB | 8.9 MB | 228.9 KB |
| POSITIVE_SELECTION:HYPHY:PRUNE_TREE_BATCH (1) | COMPLETED | 154ms | 101ms | 113.4% | 9.2 MB | 39.1 MB | 1.6 MB | 4.9 KB |
| TRANSLATION_REPORT:QUALITY_REPORT_DATA | COMPLETED | 1.2s | 1.2s | 105.1% | 27.3 MB | 2.7 GB | 25.7 MB | 2.5 KB |
| TRANSLATION_REPORT:GENERATE_METADATA | COMPLETED | 1.3s | 1.2s | 132.2% | 20.7 MB | 2.7 GB | 25.6 MB | 3 KB |
| TRANSLATION_REPORT:RENDER_TRANSLATION_REPORT | COMPLETED | 14.4s | 14.4s | 99.0% | 490.9 MB | 1.1 TB | 130.4 MB | 24 MB |
| POSITIVE_SELECTION:HYPHY:HYPHY_BATCH (1) — **FUBAR, 970 codons / 78 taxa** | COMPLETED | 34.4s | 34.3s | 109.2% | 62 MB | 484.1 MB | 1.3 MB | 25 MB |
| POSITIVE_SELECTION:HYPHY:AGGREGATE_FIRST_PASS | COMPLETED | 948ms | 891ms | 192.5% | 3.6 MB | 13 MB | 33.3 MB | 51.3 KB |
| PSEL_REPORT:RENDER_PSEL_REPORT (1) | COMPLETED | 11.8s | 11.7s | 100.7% | 495.7 MB | 1.1 TB | 139.7 MB | 33 MB |

FUBAR itself (`HYPHY_BATCH`) is 34.4s of the 57s total; the two Rmd-render steps (`RENDER_TRANSLATION_REPORT`, `RENDER_PSEL_REPORT`) account for most of the rest. `peak_vmem` values of ~1 TB on the Rmd-render steps are R/pandoc shared-library mmap accounting, not real resident memory (see `peak_rss` for actual usage, ≤500 MB throughout).

## Run 2 — PhyloPhere pipeline (CAAS → CT_DISAMBIGUATION → FADE → SCORING → ENRICHMENT)

Cold start: `output/pepc/results/` and `output/pepc/work/` removed before launch, `RESUME=0`. Launcher: `validation/tier1/input/pepc/scripts/run_tier1_pepc_local_complete.sh` (env `phylophere`), which dispatches `run_tier1_pepc_single_complete.sh c4 ... ordinal`. Modules enabled per `gui/templates/tier1_pepc_c4.json`: CAAS, disambiguation, FADE, SCORING, ENRICHMENT (FCS only, POSENRICH off); Accumulation/RER/VEP off.

**Total wall-clock** (Nextflow's own summary):

| | |
|---|---|
| Duration | 11m 40s |
| CPU hours | 1.5 |
| Processes succeeded | 40 / 40 |

**Per-process trace** (`results/c4_complete/nextflow_trace_c4_complete_8bde5a69.tsv`, 40 processes, all COMPLETED), in execution order:

| name | status | duration | realtime | %cpu | peak_rss | peak_vmem | rchar | wchar |
|---|---|---|---|---|---|---|---|---|
| CONTRAST_SELECTION:NAME_CURATION:TREE_CLEANUP (curate tree tip labels) | COMPLETED | 286ms | 226ms | 99.3% | 23.3 MB | 46.6 MB | 4.5 MB | 9.1 KB |
| CONTRAST_SELECTION:REPORTING:NAME_CURATION:TREE_CLEANUP (curate tree tip labels) | COMPLETED | 404ms | 252ms | 85.2% | 22.7 MB | 46.4 MB | 4.5 MB | 9.1 KB |
| CONTRAST_SELECTION:REPORTING:DATASET_EXPLORATION (dataset_exploration) | COMPLETED | 5s | 4.9s | 114.3% | 434.2 MB | 2.3 GB | 26.6 MB | 3.2 MB |
| CONTRAST_SELECTION:REPORTING:PHENOTYPE_EXPLORATION (phenotype_exploration) | COMPLETED | 36.5s | 36.5s | 98.6% | 2 GB | 1 TB | 50.8 MB | 32.9 MB |
| CONTRAST_SELECTION:CI_COMPOSITION_REPORT (CI_COMPOSITION_REPORT) | COMPLETED | 10.3s | 10.2s | 108.0% | 903.6 MB | 2.8 GB | 48 MB | 23.8 MB |
| SELECTION_PREP:EXTRACT_EXTREME_SPECIES (extract_extreme_species) | COMPLETED | 140ms | 78ms | 91.4% | 16.9 MB | 39.3 MB | 1.3 MB | 2.1 KB |
| CONTRAST_SELECTION:CONTRAST_ALGORITHM (CONTRAST_ALGORITHM) | COMPLETED | 23.9s | 23.8s | 103.5% | 1.2 GB | 3.1 GB | 63.4 MB | 49.2 MB |
| CONTRAST_SELECTION:CHECK_MIN_CONTRASTS (CHECK_MIN_CONTRASTS) | COMPLETED | 94ms | 20ms | 145.5% | 0 | 0 | 245.9 KB | 81.2 KB |
| SELECTION_PREP:PREP_ALIGNMENTS_BATCHED (prep_batch_00001) | COMPLETED | 329ms | 269ms | 90.9% | 7.8 MB | 25.4 MB | 4.7 MB | 154.7 KB |
| FADE:ANNOTATE_TREE_FG_BATCHED (annotate_batch_top_00001 (1 genes)) | COMPLETED | 323ms | 245ms | 102.3% | 7.9 MB | 25.4 MB | 4.8 MB | 81.9 KB |
| CT:DISCOVERY_BATCHED (discovery_batch_00001 (1 genes)) | COMPLETED | 30.1s | 3s | 98.4% | 68.6 MB | 343.2 MB | 8.8 MB | 2.4 MB |
| CT:CONCAT_DISCOVERY (Concatenating discovery outputs) | COMPLETED | 199ms | 47ms | 91.1% | 0 | 0 | 7.5 MB | 2.4 MB |
| CT:CONCAT_BACKGROUND (Concatenating background outputs) | COMPLETED | 176ms | 72ms | 79.2% | 0 | 0 | 229.5 KB | 6.1 KB |
| CT_META_CAAS:CAAS_META_CAAS_REPORT | COMPLETED | 6.5s | 6.4s | 106.9% | 300.2 MB | 1 TB | 24.6 MB | 3 MB |
| CT_DISAMBIGUATION:CT_DISAMBIGUATION_RUN (ct_disambiguation) | COMPLETED | 27.2s | 27s | 52.3% | 1.1 GB | 7.3 GB | 23.8 GB | 13.5 MB |
| CT_POSTPROC:CAAS_PREPARE_POSTPROC_INPUT (prepare_postproc_input) | COMPLETED | 1.7s | 1.6s | 122.9% | 124.8 MB | 2.9 GB | 26.4 MB | 322.1 KB |
| CT_POSTPROC:CT_FILTER (filter:3x70) | COMPLETED | 1.1s | 981ms | 163.7% | 23.9 MB | 1.7 GB | 24 MB | 822.8 KB |
| CT_POSTPROC:CT_FILTER_SUMMARY (filter_summary) | COMPLETED | 889ms | 822ms | 163.3% | 28 MB | 1.7 GB | 23.7 MB | 543 B |
| CT_POSTPROC:CAAS_FILTER_GENES (gene_filter:none) | COMPLETED | 1.1s | 1s | 147.1% | 25.7 MB | 1.7 GB | 24 MB | 322.1 KB |
| CT_POSTPROC:CAAS_BACKGROUND_CLEANUP (bg_cleanup) | COMPLETED | 1.1s | 988ms | 159.9% | 25.5 MB | 1.7 GB | 23.6 MB | 456 B |
| CT_POSTPROC:ASR_ROBUSTNESS:ASR_ROBUSTNESS_REPORT (asr_robustness) | COMPLETED | 11.7s | 11.6s | 98.2% | 593.3 MB | 1 TB | 31.7 MB | 3.6 MB |
| CT_POSTPROC:CT_POSTPROC_REPORT (caas_postproc_report) | COMPLETED | 49.7s | 49.6s | 101.4% | 772.7 MB | 2.9 GB | 28.5 MB | 12.6 MB |
| FADE:FADE_BATCHED (fade_batch_top_00001 (1 genes, top)) | COMPLETED | 2m 29s | 2m 29s | 100.5% | 485.3 MB | 734.6 MB | 1.2 MB | 1.6 GB |
| FADE:FADE_JSON_TO_CSV_TOP (fade_json_to_csv\|top) | COMPLETED | 839ms | 757ms | 131.6% | 36 MB | 1.7 GB | 10 MB | 1.2 KB |
| FADE:FADE_REPORT_TOP (fade_report\|top) | COMPLETED | 6.3s | 6.1s | 91.6% | 363.2 MB | 1 TB | 23.6 MB | 1.9 MB |
| FADE:FADE_GENE_LISTS_TOP (fade_gene_lists\|top) | COMPLETED | 431ms | 351ms | 212.7% | 33.6 MB | 1.7 GB | 6 MB | 4.6 KB |
| ENRICHMENT:PUBLISH_UNIVERSES | COMPLETED | 161ms | 26ms | 115.4% | 0 | 0 | 134.4 KB | 284 B |
| CT:RESAMPLE (nw_tree.nwk) | COMPLETED | 3m | 3m | 420.8% | 2.3 GB | 14.3 GB | 77.3 MB | 95.8 MB |
| CT:CONCAT_RESAMPLE (Concatenating resample outputs) | COMPLETED | 90ms | 28ms | 111.6% | 0 | 0 | 914.4 KB | 250.8 KB |
| CT:CAAS_PERMS_PREP:SUBSET_RESAMPLE_PERMS (caas_perms_subset\|N=1000) | COMPLETED | 708ms | 660ms | 103.0% | 8.1 MB | 27.3 MB | 112.4 MB | 61.1 MB |
| CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED (perm_replay_batch_00001 (1 genes)) | COMPLETED | 1m | 55.4s | 101.7% | 1.7 GB | 2.1 GB | 25.8 MB | 80.6 MB |
| CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE (caas_perms_disambiguate) | COMPLETED | 5m 21s | 5m 21s | 22.3% | 2.9 GB | 12.1 GB | 96.3 MB | 2.2 GB |
| CAAS_PERMULATION:CAAS_PERMS_AGGREGATE (caas_perms_aggregate) | COMPLETED | 641ms | 551ms | 281.2% | 3.9 MB | 12.7 MB | 21.5 MB | 8.2 MB |
| SCORING:SCORING_COMPUTE (scoring_compute\|c4) | COMPLETED | 1.3s | 1.2s | 187.6% | 142.6 MB | 2.7 GB | 15.1 MB | 1.9 MB |
| CAAS_SIGNIFICANCE_REPORT (1) | COMPLETED | 3.7s | 3.6s | 120.7% | 249.6 MB | 2.7 GB | 24.7 MB | 3 MB |
| ENRICHMENT:FCS_COMPUTE_SCORING:FCS_COMPUTE_BATCHED (fcs_batch_001 (4 GMTs)) | COMPLETED | 6.8s | 6.7s | 107.7% | 439.3 MB | 2.4 GB | 50.8 MB | 139.3 KB |
| ENRICHMENT:FCS_COMPUTE_SCORING:FCS_CONCAT (Concatenating FCS batch outputs) | COMPLETED | 115ms | 14ms | 112.5% | 0 | 0 | 136.6 KB | 497 B |
| ENRICHMENT:SCORING_FCS_REPORT (scoring_fcs\|c4) | COMPLETED | 4.9s | 4.8s | 118.9% | 422.9 MB | 2.7 GB | 58 MB | 1.9 MB |
| SCORING:SCORING_REPORT (scoring_report\|c4) | COMPLETED | 15.6s | 15.5s | 106.2% | 330.9 MB | 2.9 GB | 34.6 MB | 15.6 MB |
| ENRICHMENT:SCORING_COMPARE_REPORT (scoring_compare\|c4) | COMPLETED | 4.5s | 4.4s | 124.6% | 395 MB | 1 TB | 27.7 MB | 2.8 MB |

**Dominant cost centers** (by duration): `CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE` (5m 21s, the CAAS-null replay disambiguation pass — 1,000,000-draw permulation budget per `N_RANDOMIZATIONS`), `CT:RESAMPLE` (3m, peak RSS 2.3 GB, peak `%cpu` 420.8% — the only clearly multi-threaded step), `FADE:FADE_BATCHED` (2m 29s), `CT:CAAS_PERMS_PREP:PERM_REPLAY_BATCHED` (1m). Together these four account for ~11m40s of the ~11m40s total, i.e. essentially the entire run. Peak single-process RSS across the run was 2.9 GB (`CAAS_PERMS_DISAMBIGUATE`); `peak_vmem` values of ~1 TB (four Rmd-render steps) are again shared-library mmap accounting, not real usage. `CT_DISAMBIGUATION_RUN`'s `rchar` of 23.8 GB against a 27.2s runtime (≈875 MB/s sustained) is far above what a single-pass read of this fixture's inputs would produce; not independently verified against the process's I/O pattern, so reported as-is rather than attributed to a specific cause.
