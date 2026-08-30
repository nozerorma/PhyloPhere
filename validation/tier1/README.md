# Tier 1 — site-level truth sets

**Question:** on compact, well-studied datasets with a curated list of supported
convergent sites, does PhyloPhere recover them and where do they rank?

**Gate:** validated positives recovered at the default significance threshold and
sitting in the top decile of the position score. For echolocation additionally:
the permulation/bootstrap null suppresses the genome-wide inflation
(Thomas & Hahn 2015 / Zou & Zhang 2015 expectation).

## Datasets and the trait cells they fill

| dataset | spec | type | count | grouping |
|---------|------|------|-------|----------|
| PEPC / C4 | `fixtures/tier1/pepc/pepc.spec.json` | categorical | single | direct |
| RH1 spectral tuning | `fixtures/tier1/rh1/rh1.spec.json` | continuous | single | direct + percentilized |
| echolocation | `fixtures/tier1/echolocation/echo.spec.json` | categorical | single | direct |
| echolocation + marine | `fixtures/tier1/echolocation/echo_marine.spec.json` | categorical | multi (n_trait=2) | direct |

RH1 fills the continuous cells (photic depth as a metric trait, then top/bottom-q
percentilized into a contrast); the echo+marine spec fills `n_trait=2` with
partially overlapping foregrounds, which is where the local multi-phenotype loop
has broken before.

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
