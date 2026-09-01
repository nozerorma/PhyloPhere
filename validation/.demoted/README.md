# Demoted — Tier 0 (simulation)

**Status: demoted, unmaintained.** Decision D-DIR-01 (2026-08-31), see
`../docs/DECISIONS.md`.

## Why

Tier 0 planted convergent substitutions on simulated genes and asked whether the
pipeline recovered them and kept its null calibrated. The design could not be
made both non-circular and informative at the same time:

- To make planted signal **detectable**, the observed trait has to point clearly
  at the planted lineages — which makes "did contrast selection recover the true
  foreground" tautological.
- Leaving it **non-circular** (latent Brownian trait + observation noise)
  collapses detection: the pipeline's contrast selector prefers shallow
  congeneric high/low pairs that miss the planted deep lineages entirely. In the
  last smoke run (`primate_half`, λ=0.5), binary and rate archetypes had
  AUC ≈ 0.5 and zero site recall.

Precedent: CAAP (Chen et al. 2025, Mol Ecol Resour, doi:10.1111/1755-0998.70052)
validates with only a **minimal** neutral-simulation null (confirm the
convergence statistic stays at its neutral expectation) and does all
recovery/power evidence on real datasets. The suite now follows that: Tiers 1–3
on real data, pipeline unchanged.

## Contents

```
tier0/              the simulation package (model, trees, simulate, pheno,
                    replicate, build_project, run_replicates, plots, groups,
                    grid.json, reproduce.sh) + gate_result/ (last reference run)
tier0_adapter.py    the runs/ -> gate-metrics adapter (was harness/tier0_adapter.py)
fixtures_tier0/     tree-fixture README
test_tier0.py.disabled   the unit tests (renamed so pytest skips them)
```

Package imports are broken by the move (`.demoted` is not an importable name).
Restore under `validation/tier0/` and `validation/harness/` to run it.

## If Tier 0 is revived

The realistic scope is the **minimal neutral-sim null only** — N neutral
alignments through CAAS discovery + permulation, confirm the genome-wide p-value
distribution stays flat. No phenotype, no planting, no power claim. The
`D-T0-*` decisions in `DECISIONS.md` record every planting variant that was
tried; read them before attempting a power test again.
