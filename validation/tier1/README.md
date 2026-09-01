# Tier 1 — site-level truth sets

**Question:** on compact, well-studied datasets with a curated list of supported
convergent sites, does PhyloPhere recover them and where do they rank?

**Gate:** validated positives recovered at the default significance threshold and
sitting in the top decile of the position score.

Tier 1 is deliberately simple: take datasets where a published study demonstrated
convergence and *externally validated* specific genes / sites (mutagenesis,
kinetics, expression assays), run the pipeline unchanged, and check we recover
them and where they rank. No null-calibration or genome-wide-inflation analysis
here — that debate (Thomas & Hahn 2015, Zou & Zhang 2015) lives in the CAAP
comparison track, not Tier 1.

## Datasets and the trait cells they fill

| dataset | spec | type | count | grouping |
|---------|------|------|-------|----------|
| PEPC / C4 | `fixtures/tier1/pepc/pepc.spec.json` | categorical | single | direct (+ genotypic/phenotypic → n_trait=2) |
| **Hb / high altitude** | `fixtures/tier1/hb_altitude/hb_altitude.spec.json` | **continuous** + categorical | single | direct + percentilized |

Two fixtures, each anchored by mutagenesis-validated convergent sites — **PEPC**
780/665, **Hb** αA34T/αA119A — and each with real contrast headroom (23 C4
origins / ~5 independent high-altitude tit lineages).

**Hb / high altitude** (Zhu et al. 2018, Sino-Himalayan tits; D-DIR-03) fills the
continuous cells — elevational-range midpoint as a metric trait, then top/bottom-q
percentilized — and also `categorical · single · direct` via a 2500 m H/L cut.
The `n_trait=2` code path is exercised by **PEPC's genotypic vs phenotypic C4
annotation** (two real trait columns, Task 4 runs both), not a separate fixture.

Dropped from Tier 1: the fish-RH1 fixture (D-DIR-04 — Otophysi all-freshwater, a
paired-contrast method reaches only ~4 of ~23 F261Y origins); echolocation
(D-DIR-05 — 3 origins, at the `min_contrasts` floor, no accessible 3-origin
dataset); the `echo_marine` n_trait=2 fixture (redundant with Tier 2). Echo +
low-origin convergence revisited when Tier 2 is scoped, via explicit branch
pairs. `../truthsets/tier1/echolocation.*.tsv` kept parked for that.

## Fixtures (not committed — build or fetch)

Each `fixtures/tier1/<name>/` needs: `align/` (codon-aware CDS alignment, one file
per gene, tips = species), `tree.nwk` (species tree with branch lengths, same
tips), `<name>.spec.json` (the `PhenotypeSpec`; commit this — it is small).

| dataset | alignment source | tree |
|---------|------------------|------|
| PEPC | PCOC repo `data/` (sedge PEPC alignment, Besnard 2009); or realign the Besnard 2009 supplement with MACSE | gene tree from the alignment (IQ-TREE) |
| RH1 | Hauser et al. 2017 supplementary alignment; or realign fish RH1 CDS from NCBI with MACSE | species tree (fish ToL / actinopt) pruned to tips |
| echolocation | Parker et al. 2013 alignments (Dryad); the two reanalyses use the same set | mammal species tree pruned to the ~22 taxa |

Put a `bovine_RHO` / `maize_PEPC` reference row in the respective alignment so
`truthset.map_site_to_alignment` can place the truth-set coordinates.

## Run

```bash
# 1. emit the trait matrix from the spec
python -m validation.harness.cli emit \
    validation/fixtures/tier1/rh1/rh1.spec.json \
    --out validation/runs/tier1/rh1/traits

# 2. run PhyloPhere per trait cell (categorical -> CAAS+FADE, continuous -> RER+CAAS/BM)
#    model runner scripts on run_scripts/run_integrated_toy.sh; one invocation per
#    cell under runs/tier1/rh1/traits/**. TODO: tier1/run_rh1.sh

# 3. score
python -m validation.harness.cli score \
    --truth validation/truthsets/tier1/rh1_spectral.sites.tsv \
    --run   validation/runs/tier1/rh1/<cell>/ \
    --ref-row bovine_RHO
```

`harness/cli.py` is a TODO stub — the `emit` path is wired, `score` needs the
PhyloPhere-output adapter (which output table, which column is the position, which
the score). Fill it once the first real run exists so the adapter matches reality.
