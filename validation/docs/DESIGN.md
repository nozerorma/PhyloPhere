# Validation design

> **Tier 0 is DEMOTED** (decision D-DIR-01, see DECISIONS.md). Its
> planted-simulation design is irreducibly circular. Code lives unmaintained
> under `validation/.demoted/`. All Tier 0 rows/sections below are kept for the
> record and marked ~~struck~~ or "(demoted)". The live suite is Tiers 1–3 on
> real data, pipeline unchanged.

## Principle

Tiered, fail-fast. Each tier answers a distinct question; a failure at tier *n*
makes tier *n+1* uninterpretable. Live run order: **3 → 1 → 2**.

- **3 first** — nearly free, tells us immediately whether the CAAS core is sound.
- **1 next** — locks site-level accuracy; CI-scale, runs on every CAAS/FADE change.
- **2 last** — the expensive genome-wide integration run.
- ~~**0 before 2** — the null must be calibrated before Tier 2's numbers mean anything.~~ (demoted; a minimal neutral-sim null check may return later, no planting)

## What each tier validates in PhyloPhere

| Tier | Modules exercised | Gate |
|------|-------------------|------|
| ~~0~~ (demoted) | CAAS permulation/bootstrap null, RER permulation, FADE calibration, ASR path-score edge cases, SCORING degradation paths | ~~KS uniformity of null p-values; type-I error within MC error; monotone power curves; zero crashes on pathological tree~~ — demoted (D-DIR-01) |
| 1 | CAAS/CAAP position calls, position-FCS, FADE sites, CT_DISAMBIGUATION, position→gene aggregation | validated positives recovered at default threshold and in top decile of the position score; echolocation genome-wide inflation suppressed by the null |
| 2 | full DAG incl. SCORING joins, per-module FCS universes, ENRICHMENT (FCS/STRING/DOMINO) | rank concordance vs RERconverge/CSUBST (Spearman ρ > 0.3); consensus biology in top enrichment; per-module gene counts correct; subsampled null still ~uniform |
| 3 | CAAS discovery + signification (anchor 1); RER module (anchor 2) | CAAS calls on the CAAStools lifespan `/test` set overlap the reference at Jaccard ≥ 0.9 after matching settings; RER accelerated-gene lists overlap published lists; eye/lens genes rank high in subterranean run |

## Datasets

| Tier | Dataset | Phenotype axes covered | Truth | Source |
|------|---------|------------------------|-------|--------|
| ~~0~~ (demoted) | simulated (WAG + site-specific Dirichlet profiles, Gillespie; planted convergence) | — | generative | ~~primate/mammal/star/ladder trees~~. Demoted, D-DIR-01. Code: `validation/.demoted/tier0/` |
| 1 | PEPC / C4 in sedges (Cyperaceae) | categorical · single | 7 tiered sites (780/665 mutagenesis) + Besnard-16 | Besnard 2009 via PCOC/ConDor `test_data`. Fixture built (`fixtures/tier1/pepc/`). Comparators: PCOC, ConDor, FADE (Morel 2024) |
| 1 | RH1 spectral tuning in fish | continuous · single · direct **and** percentilized | Yokoyama 2008 mutagenesis λmax-shift sites (bovine numbering) | broad teleost RH1 alignment w/ bovine RHO ref, many independent dim-light transitions (option #3 — pick a specific alignment). NOT "Hauser 2017" (single origin) |
| 1 | echolocation (bats + toothed whales) | categorical · single; multi (n_trait=2) | 7 Parker 2013 hearing genes + neutral-inflation expectation | **CAAP setup** (Chen 2025): OrthoMaM v10c, BP1/BP2 branch pairs, Bovidae negative control |
| 2 | marine mammals (Foote 2015) | categorical · single; continuous (dive depth) · percentilized | RERconverge (Chikina 2016) / CSUBST / **CAAP (Chen 2025)** published gene rankings + GO terms | **CAAP setup**: Foote 2015, 22 sp, walrus/dolphin/orca/manatee, terrestrial-sister negative control |
| 3 | mammalian lifespan (Farré 2021) | continuous · single · percentilized (native top/bottom) | CAAStools `/test` CAAS calls | CAAStools repo |
| 3 | subterranean mammals (eye degeneration) | categorical · single | published RER accelerated genes; crystallin/lens expectation | Partha 2017 |

## Trait matrix coverage

Every tier must cover the cross-product of {categorical, continuous} ×
{single, multi} × {direct, percentilized}. ~~Tier 0 covers it exhaustively by
construction;~~ Tiers 1–3 cover it across datasets (see the table above — each
dataset is chosen partly for which cells it fills). `harness/trait_matrix.py`
is the single code path that emits all of it, so the coverage lives in the
fixtures, not in duplicated glue.

## Limits (state these in any write-up)

- Echolocation "truth" is contested — it is a good test of the null, **not** a
  clean precision benchmark.
- Site-level mutagenesis truth sets are small (n ≈ 5–20); recall CIs are wide.
- Simulation truth depends on the planted-signal model resembling real convergent
  phenotypes, which is unknown. Mitigation: simulate under a **different model
  family** than the scorer uses, and sweep misspecification.
- None of these tiers validates a *novel* PhyloPhere hit. They validate that the
  machinery reproduces what is already established and controls error. That is a
  defensible validation story for a methods paper / thesis chapter; novel-hit
  validation is a separate (wet or orthogonal-data) problem.

## Review that frames the field

Allard JB, Kumar S. 2026. "The genetic foundations of convergent traits."
Nat Rev Genet 27(7):563–578. doi:10.1038/s41576-026-00933-7. Cites CAAStools,
PCOC, CSUBST, RERconverge, and the Kumar-lab ESL line. (Full-text package table
not yet reconciled against this suite — paywalled, no preprint.)
