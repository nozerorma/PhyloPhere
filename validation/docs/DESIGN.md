# Validation design

## Principle

Tiered, fail-fast. Each tier answers a distinct question; a failure at tier *n*
makes tier *n+1* uninterpretable. Run in order: **3 → 1 → 0 → 2**.

- **3 first** — nearly free, tells us immediately whether the CAAS core is sound.
- **1 next** — locks site-level accuracy; CI-scale, runs on every CAAS/FADE change.
- **0 before 2** — the null must be calibrated before Tier 2's numbers mean anything.
- **2 last** — the expensive genome-wide integration run.

## What each tier validates in PhyloPhere

| Tier | Modules exercised | Gate |
|------|-------------------|------|
| 0 | CAAS permulation/bootstrap null, RER permulation, FADE calibration, ASR path-score edge cases (star topology, empty private segment), SCORING degradation paths | KS uniformity of null p-values (p > 0.05); type-I error within MC error of nominal; monotone power curves; zero crashes on pathological tree |
| 1 | CAAS/CAAP position calls, position-FCS, FADE sites, CT_DISAMBIGUATION, position→gene aggregation | validated positives recovered at default threshold and in top decile of the position score; echolocation genome-wide inflation suppressed by the null |
| 2 | full DAG incl. SCORING joins, per-module FCS universes, ENRICHMENT (FCS/STRING/DOMINO) | rank concordance vs RERconverge/CSUBST (Spearman ρ > 0.3); consensus biology in top enrichment; per-module gene counts correct; subsampled null still ~uniform |
| 3 | CAAS discovery + signification (anchor 1); RER module (anchor 2) | CAAS calls on the CAAStools lifespan `/test` set overlap the reference at Jaccard ≥ 0.9 after matching settings; RER accelerated-gene lists overlap published lists; eye/lens genes rank high in subterranean run |

## Datasets

| Tier | Dataset | Phenotype axes covered | Truth | Source |
|------|---------|------------------------|-------|--------|
| 0 | simulated (empirical tree + planted shifts) | all cells of the trait matrix, by construction | generative | pyvolve/INDELible over a pruned Zoonomia tree + a star-topology pathological tree |
| 1 | PEPC / C4 in sedges & grasses | categorical · single | ~21 manual sites (pos 780 validated) | Besnard 2009; Christin 2007; PCOC repo alignment |
| 1 | RH1 spectral tuning in fish | continuous · single · direct **and** percentilized; categorical recode | ~12–30 mutagenesis/manual sites | Hauser 2017; Yokoyama 2008 |
| 1 | echolocation (bats + toothed whales) | categorical · single; multi (n_trait=2) when paired with marine | 6 hearing genes + neutral-inflation expectation | Parker 2013 + the two reanalyses |
| 2 | marine mammals (Zoonomia subset / Chikina 2016) | categorical · single; continuous (dive depth) · percentilized | RERconverge / CSUBST / CAAP published gene rankings | Chikina 2016; Foote 2015; Zoonomia |
| 3 | mammalian lifespan (Farré 2021) | continuous · single · percentilized (native top/bottom) | CAAStools `/test` CAAS calls | CAAStools repo |
| 3 | subterranean mammals (eye degeneration) | categorical · single | published RER accelerated genes; crystallin/lens expectation | Partha 2017 |

## Trait matrix coverage

Every tier must cover the cross-product of {categorical, continuous} ×
{single, multi} × {direct, percentilized}. Tier 0 covers it exhaustively by
construction; Tiers 1–3 cover it across datasets (see the table above — each
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
