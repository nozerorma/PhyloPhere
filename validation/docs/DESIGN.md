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
| 1 | CAAS/CAAP position calls, position-FCS, FADE sites, CT_DISAMBIGUATION, position→gene aggregation | externally-validated convergent genes/sites recovered at default threshold and in the top decile of the position score (recovery only — no null calibration at Tier 1) |
| 2 | full DAG incl. SCORING joins, per-module FCS universes, ENRICHMENT (FCS/STRING/DOMINO) | rank concordance vs RERconverge/CSUBST (Spearman ρ > 0.3); consensus biology in top enrichment; per-module gene counts correct; subsampled null still ~uniform |
| 3 | CAAS discovery + signification (anchor 1); RER module (anchor 2) | CAAS calls on the CAAStools lifespan `/test` set overlap the reference at Jaccard ≥ 0.9 after matching settings; RER accelerated-gene lists overlap published lists; eye/lens genes rank high in subterranean run |

## Datasets

| Tier | Dataset | Phenotype axes covered | Truth | Source |
|------|---------|------------------------|-------|--------|
| ~~0~~ (demoted) | simulated (WAG + site-specific Dirichlet profiles, Gillespie; planted convergence) | — | generative | ~~primate/mammal/star/ladder trees~~. Demoted, D-DIR-01. Code: `validation/.demoted/tier0/` |
| 1 | PEPC / C4 in sedges (Cyperaceae) | categorical · single | 7 tiered sites (780/665 mutagenesis) + Besnard-16 | Besnard 2009 via PCOC/ConDor `test_data`. Fixture built (`fixtures/tier1/pepc/`). Comparators: PCOC, ConDor, FADE (Morel 2024) |
| 1 | **haemoglobin / high altitude** (Sino-Himalayan tits) | **continuous · single · direct + percentilized**; categorical · single · direct | αA **A34T** + αA **P119A** (site-directed mutagenesis, Storz lab) | **Zhu et al. 2018 PNAS** (D-DIR-03). GenBank MG772099–MG772439, αA/αD/βA globin, 16 tips, elevation. Continuous `elev_mid` is our encoding (Zhu analyse H/L). Fixture: `fixtures/tier1/hb_altitude/`. No PCOC/ConDor comparator here. |
| 1 | ~~RH1 in fish~~ | — | — | **DROPPED** (D-DIR-04). Hill et al. 2019 F261Y is real (~23 origins) but Otophysi is all-freshwater → a paired-contrast method reaches only ~4 origins. "Can't test what we wanted." |
| ~~1~~ | ~~echolocation~~ | — | Prestin N7T/I384T (parked in `truthsets/tier1/echolocation.*.tsv`) | **DROPPED from Tier 1** (D-DIR-05): 3 origins = at the `min_contrasts` floor; Thomas & Hahn data has only 2; no accessible 3-origin set. Revisit at Tier 2 via explicit branch pairs. |
| 2 | marine mammals (Foote 2015) | categorical · single; continuous (dive depth) · percentilized | RERconverge (Chikina 2016) / CSUBST / **CAAP (Chen 2025)** published gene rankings + GO terms | **CAAP setup**: Foote 2015, 22 sp, walrus/dolphin/orca/manatee, terrestrial-sister negative control |
| 3 | mammalian lifespan (Farré 2021) | continuous · single · percentilized (native top/bottom) | CAAStools `/test` CAAS calls | CAAStools repo |
| 3 | subterranean mammals (eye degeneration) | categorical · single | published RER accelerated genes; crystallin/lens expectation | Partha 2017 |

## Trait matrix coverage

`harness/trait_matrix.py` is the single code path that emits the cross-product of
{categorical, continuous} × {single, multi} × {direct, percentilized} from one
spec. **Tier 1 (PEPC + Hb)** covers: categorical·single·direct (both);
categorical·single·percentilized + continuous·single·{direct,percentilized} (Hb
`elev_mid`); n_trait=2 (PEPC genotypic vs phenotypic C4). The
continuous·**multi** and categorical·**multi** cells are **not** covered at
Tier 1 — they move to Tier 2/3 (marine-mammal dive depth; lifespan). This is a
narrower Tier 1 than first planned: the echo/RH1 fixtures that would have added
`n_trait=2` and more categorical origins were dropped (D-DIR-04/05) because the
pipeline's contrast-selection floor (≥3 well-separated origins with accessible
sisters) excludes them.

## Limits (state these in any write-up)

- Echolocation *genome-wide* convergence is contested (Parker 2013 vs Thomas &
  Hahn 2015 / Zou & Zhang 2015), but the *specific* prestin sites (N7T, I384T)
  are functionally validated (Liu 2014) — those are the Tier 1 targets. Genes
  like TMC1/CDH23/OTOF enter only as weak reanalysis-tier positives.
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
