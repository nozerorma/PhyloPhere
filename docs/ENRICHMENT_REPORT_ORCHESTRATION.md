# Enrichment report orchestration — Option A + COMPARE cross-FCS overview

> Status: **planned** (not started). Captured 2026-06-26. Approved direction:
> Option A (per-module FCS reports downstream of SCORING, auto-annotated) plus a
> generalized COMPARE_REPORT giving a cross-FCS overview of enriched terms.

## Why

Cross-module leading-edge flags (CAAS gates, top/bottom percentiles, FADE, accum)
are produced by SCORING (`scoring_compute.R` fuses all modules into `fcs_stats.tsv`).
But `RER_FCS_REPORT` / `FADE_FCS_REPORT` / `ACCUMULATION_FCS_REPORT` run inside their
own tool branches — BEFORE SCORING — so they can't see that annotation. Today the
only fix is the manual `--rer_gene_scores` param (commit c811311) pointing at a prior
`fcs_stats.tsv` (two-pass). The dependency is structural: the annotation IS scoring's
output, so the consumer must move later, not the data earlier.

Confirmed ordering (main.nf): `RER_MAIN` (runs RER FCS internally, line ~338) →
`SCORING` (line ~364, consumes `RER_MAIN.out.summary_tsv`).

## Option A — per-module FCS reports downstream of SCORING (auto-annotated)

When SCORING is in the run, defer the per-module FCS reports to AFTER it and auto-feed
`SCORING.out.fcs_stats` as the `gene_scores` annotation input — reusing the optional
cross-module join already built in `RER_GENE_LISTS`. No `--rer_gene_scores` needed.

- Reorder: `SCORING` → `{RER,FADE,ACCUMULATION}_FCS_REPORT` (annotated) in main.nf.
- Extend the optional `gene_scores` join (from RER_GENE_LISTS) to `FADE_GENE_LISTS`
  and `ACCUMULATION_GENE_LISTS` (same `want` column list, same NO_FILE sentinel).
- Auto-wire `SCORING.out.fcs_stats` into all three when scoring runs; keep the
  in-branch RER/FADE/accum-only path for STANDALONE runs (tool alone, no scoring).
- Touch points: `main.nf` (orchestration order + pass fcs_stats), workflow take/emit
  signatures (rerconverge.nf / fade.nf / ct_accumulation.nf), the FADE/accum gene-list
  producers (add the join), conf comments.

## COMPARE_REPORT — cross-FCS overview of enriched terms (runs LAST)

Generalize `SCORING_COMPARE_REPORT` (`COMPARE_scoring.Rmd`) from "CAAS top vs bottom"
to a META-overview across ALL module FCS outputs, downstream of every FCS report.

- Input: every module's `fcs_all_results.tsv` (CAAS, RER, FADE, accum), tagged by module.
- Output: a term × module enrichment view — which pathways are enriched (stat > 0,
  passing the dual threshold) in which module; HIGHLIGHT multi-module convergence
  (terms significant across ≥2 modules = strongest, independent-evidence pathways).
  Suggested: heatmap/dotplot of -log10(p.adj or p.perm) per term per module + a table
  of convergent terms ranked by how many modules support them.
- Touch points: `scoring_enrichment.nf` (collect all module `fcs_all_results` into the
  COMPARE process), `COMPARE_scoring.Rmd` (the cross-module overview logic), `main.nf`.

## Final sequencing

SCORING → per-module FCS reports (auto cross-annotated) → COMPARE_REPORT (meta-overview).
All downstream of the scoring join point. Standalone single-tool runs keep their
module-only reports unchanged.

## Reuses (already built)

- Optional `gene_scores` cross-module join — `RER_GENE_LISTS` (c811311); generalize.
- Vectorized FCS engine + one-sided alternative + dual threshold (e73ecd1).
- Verified column names: flag_gate_sig/_fdr, score_top/bottom, flag_fade[_top/_bottom],
  flag_accum live in SCORING's `fcs_stats.tsv` (scoring_compute.R:990).

## Not to break

- Standalone runs (e.g. `--rer_tool` without `--scoring`) must still produce the RER
  FCS report with RER-only annotation (gates/fade/accum simply absent).
- Don't double-run the per-module FCS report (in-branch AND post-scoring) — gate the
  in-branch one on "scoring not in this run".
