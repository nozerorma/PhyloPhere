# PhyloPhere validation suite

Tiered validation of the PhyloPhere pipeline against known-truth data. The tiers
fail fast and cheap: each answers a different question, and a failure at tier *n*
makes tier *n+1* uninterpretable, so they run in order and gate on each.

| Tier | Question | Truth | Cost | Status |
|------|----------|-------|------|--------|
| ~~0~~ | ~~Under a known generative process, is the null calibrated and what is the power?~~ | simulated | — | **DEMOTED** — see below |
| 1 | Does it recover experimentally/manually supported convergent **sites**? | curated site lists | minutes | in progress (PEPC fixture done) |
| 2 | Does the integrated pipeline recover the biology other tools agree on? | published tool outputs (marine mammals, CAAP/RER/CSUBST) | one cluster run | not started |
| 3 | Does the CAAS core reproduce a published analysis given the same input? | CAAStools `/test` lifespan set; RERconverge marine/subterranean | low | not started |

**Tier 0 is demoted** (`.demoted/`, decision D-DIR-01). The planted-simulation
design is irreducibly circular: making planted signal detectable forces the
observed trait to point at the planted lineages, which makes contrast-recovery
tautological; leaving it non-circular collapses detection (the pipeline's
contrast selector prefers shallow congeneric pairs that miss the planted deep
lineages). Precedent: CAAP (Chen 2025) uses only a minimal neutral-sim null and
does all recovery on real data. The suite is now Tiers 1–3 on real datasets,
pipeline unchanged. Code kept, unmaintained, under `validation/.demoted/`.

Full rationale and dataset provenance: `docs/DESIGN.md`.

## Why a shared harness

The pipeline makes three qualitatively different claims and each needs its own
truth set and metric:

- **site-level** — a CAAS/CAAP at a position, a FADE site, a position-FCS hit
- **gene-level** — a scored/ranked gene, an RER acceleration
- **set-level** — an enriched pathway / PPI module

`harness/` provides the pieces every tier reuses:

- `phenotype.py`  — `PhenotypeSpec`: the per-tip phenotype for one dataset, plus loaders
- `trait_matrix.py` — expands a `PhenotypeSpec` into the concrete input files each
  PhyloPhere trait mode consumes (see **Trait matrix** below)
- `metrics.py`    — scoring metrics (precision/recall/MCC, PR-AUC, ROC-AUC, rank of
  each known positive) **and** null-calibration metrics (KS uniformity of the
  permulation/bootstrap p-values, QQ data, observed type-I error)
- `truthset.py`   — truth-set schema + coordinate mapping (which reference
  sequence, 1-based, gap handling)

## Trait matrix

PhyloPhere's trait handling has three independent axes. Every tier must cover the
full cross-product, because bugs live in specific cells (the project history has
several: local multi-phenotype loop, contrast_selection, prune_data/prune_list,
percent_rank tie-breaking).

| axis | modes | PhyloPhere surface |
|------|-------|--------------------|
| **type** | categorical / continuous | `caas_config` groups + `perm_strategy=FGBG` vs `traitvalues` + `perm_strategy=BM\|lambda`; RER continuous vs binary |
| **count** (`n_trait`) | single / multi | one `caas_config` file vs a directory of them (`load_cfg(mode="multi")`); one vs many phenotype columns |
| **grouping** | direct / percentilized | native fg/bg groups vs top-q/bottom-q of a continuous vector fed through the contrast selection + `phen_score` percent-rank path |

`trait_matrix.py` emits, from one `PhenotypeSpec`:

```
<out>/
  categorical/
    direct/            single/  <trait>.cfg              # species \t 1|0 \t pair_id
                       multi/   <traitA>  <traitB> ...    # dir -> load_cfg(mode="multi")
    percentilized/     single/  <trait>.q<q>.cfg
                       multi/   ...
  continuous/
    direct/            single/  <trait>.traitvalues.tsv   # species \t value  (BM/lambda)
                       multi/   <trait>.traitvalues.tsv   # species \t vA \t vB ...
    percentilized/     single/  <trait>.traitvalues.tsv   # + percentile column (phen_score input)
  tree/  pruned to the tips present in the spec
  manifest.json                                           # what was emitted + provenance hash
```

Pairing note: paired mode is mandatory in this caastools fork (`pindex.load_cfg`
rejects <3-column configs). For `direct` categorical the pairs come from the
spec. For `percentilized` the harness rank-pairs by default (highest fg with
lowest bg, …) and records that in the manifest — phylogenetically-informed
pairing (what `contrast_selection` does in production) is preferred and can be
supplied explicitly in the spec.

## Layout

```
validation/
  README.md              — this file
  docs/DESIGN.md         — full tiered rationale, dataset provenance, limits
  harness/               — shared, tier-agnostic code
  truthsets/             — curated known-positive tables (version controlled, small)
    tier1/
  fixtures/              — alignments + trees (NOT committed; built/fetched, see each tier's README)
  tier1/                 — site-level truth sets: PEPC C4, RH1 spectral tuning, echolocation
  runs/                  — pipeline outputs (NOT committed)
  .demoted/              — Tier 0 (simulation): demoted, unmaintained (D-DIR-01)
```

## Running

Each tier has its own README with the exact commands. The harness is plain
Python 3.12 (stdlib + numpy); no pipeline runtime needed to compute metrics from
an existing `runs/` output.
