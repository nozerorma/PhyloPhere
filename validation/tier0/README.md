# Tier 0 — simulation: does the pipeline recover planted convergence?

**Question.** Plant convergent substitutions on lineages a phylogenetic process
picked out, hand the pipeline a noisy trait for those lineages, and ask: does it
rank those genes and sites to the top, assign them to the right direction, and
keep genes with nothing planted at the bottom — and how does all of that hold up
as the trait becomes more phylogenetically clustered?

It is a **signal-recovery check on a prioritisation engine** (`../docs/DECISIONS.md`
D-T0-P; the pipeline's permulation p-values are foreground-specificity scores,
never Uniform(0,1) — no KS gate).

## The design — evolutionary-model planting (D-T0-Q)

1. **Latent trait = Brownian motion on the tree rescaled by Pagel's lambda.**
   Same family the permulation draws its nulls from (`permulations.R` `simpermvec`).
   - `lambda = 0` — internal branches collapse to a star → latent is iid → the
     trait extremes are **scattered** across the tree.
   - `lambda = 1` — full BM → the extremes **clump into clades**.
   - `lambda = 0.5` — in between (the real cancer-prevalence trait fitted ≈ 0.74).
2. **True foreground** = the independent top- and bottom-tail origins of that
   latent draw, reduced to non-nested clusters that are mutually Dunn ≥ 1. The
   number **used is fixed** (`--n-pairs 4`, matching the pipeline's
   `contrast_max_iter = 3` → seed + 3 = 4 pairs) so replicates are comparable
   across lambda; the latent is resampled until it supports at least that many.
   How many it *could* support (`n_possible_pairs`, ~20 at lambda 0, ~5 at
   lambda 1) is **reported** as the structure-dependence diagnostic.
3. **Observed trait** handed to the pipeline = latent + sampling noise:
   - `binary` — threshold the latent → 0/1 (`--trait_type ordinal`), a few bits flipped.
   - `rate` — `c ~ Binomial(n, rate(latent_percentile))`, `n ~ U(25, 70)` → CLASS 1
     (Jeffreys CI). This is the shape of the real cancer run.
4. **Molecular convergence is planted on the TRUE foreground lineages** — not on
   whatever contrast selection outputs. Each planted gene is planted in one
   direction, `top` or `bottom` (the pipeline scores a convergent change in the
   high-trait and the low-trait species symmetrically).
5. Contrast selection runs on the *observed* trait and produces its own operative
   pairs, which differ from the true foreground. That gap, **as a function of
   lambda**, is the non-circular measurement.

Gate axes: `primate` tree × {binary, rate} × {null, power} × **lambda ∈ {0, 0.5, 1}**.
~120 genes, 350 sites, `concentration ≈ 2`, production permulation settings.
Modules: `contrast_selection` + `CAAS` + `CT_DISAMBIGUATION` + `SCORING`.

## Planted mechanisms (D-T0-02)

Per planted gene, ~12 + ~12 sites:

- **`identical_aa`** — near-point-mass on one residue → every planted lineage
  converges to the *same* residue → a strict-CAAS (`US`) positive, and a
  grouped-CAAP (`GS1–GS4`) hit for free.
- **`grouped_caap`** — ~uniform over one GS1 class, a distinct residue seeded per
  origin → origins share the *class*, differ in *residue*. GS should call it;
  whether US also calls it depends on the background (see gate result).

Generative model differs from the pipeline's inference model on purpose (WAG
exchangeabilities, per-site Dirichlet profiles, continuous Γ, Gillespie).

## Scoring (`harness score`)

`python -m validation.harness.cli score --run validation/runs/tier0`

Reported **per (archetype, lambda) cell**:

1. **Null vs power — occurrence** (the "type-I" analogue): AUC of the
   **detected-CAAS count** of planted genes vs every gene in the matched null
   replicates. `separated` = AUC ≥ 0.95 and planted p50 > null p95. Score-level
   AUCs (position `CAAS_score`, position `−log10 p`, `gene_caas_score`) are
   **reported, not gated** — `gene_caas_score` is a within-run percentile rank
   and does not separate across replicates by design.
2. **Prioritisation** — `precision@k` (top-n_planted genes by CAAS count that are
   planted), planted rank percentile, `slice_global25` membership.
3. **Site recovery** — planted-site detection recall, `directional` recall
   (detected *and* on the right top/bottom side), split by mechanism × scheme.
4. **Contrast recovery** — operative pairs (`traitfile.tab`) vs the true
   lambda-drawn pairs: anchor recall, partner recall, exact-pair rate, lineage
   Jaccard. **Reported as a lambda-curve, not gated** — expected to degrade as
   the trait clumps.

`VERDICT: PASS` when, at **every** cell, the pipeline manufactures no false
signal (occurrence separates) and, **at the lowest lambda**, `precision@k ≥ 0.85`
and `site_recall ≥ 0.8`.

## Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then one command:
validation/tier0/reproduce.sh          # -> validation/runs/tier0_gate/
#   reproduce.sh [OUT_DIR] [N_REPLICATES] [N_GENES];  LAMBDAS= JOBS= env vars

# step by step:
python -m validation.tier0.run_replicates --out validation/runs/tier0_gate \
    --archetypes binary,rate --sets null,power --trees primate \
    --lambdas 0,0.5,1 --n-replicates 10 --n-genes 120 --run --jobs 2
python -m validation.harness.cli score --run validation/runs/tier0_gate --json .../score.json
python -m validation.tier0.plots       --run validation/runs/tier0_gate
```

`3 lambda × 2 archetypes × 2 sets × 10 reps = 120 replicates`, ~10–20 min each,
memory-bound (each is a full Nextflow run).

## Tests

```bash
python -m pytest validation/tests -q
Rscript subworkflows/TRAIT_ANALYSIS/local/src/ordinal_trait_test.R
```

## Status

| piece | state |
|-------|-------|
| `model.py` / `trees.py` / `simulate.py` / `groups.py` / `pheno.py` (lambda) | done |
| `replicate.py` (top/bottom direction) / `run_replicates.py` (lambda axis) | done |
| `tier0_adapter.py` + `cli.py score` (per-cell, directional, score-level) | done |
| `plots.py` (separation, recovery-vs-lambda, score-separation, site) | done |
| binary-trait pipeline support | landed on `main` (`1fbdfbf`) |
| **lambda gate run** | **pending** — external drive was disconnected; `reproduce.sh` when back |

`gate_result/` holds the previous (no-lambda) reference — superseded, kept until
the lambda gate replaces it.
