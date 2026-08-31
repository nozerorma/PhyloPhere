# Tier 0 — simulation: does the pipeline recover planted convergence?

**Question.** If I plant convergent substitutions in genes/sites I choose, on
foreground lineages I choose, does PhyloPhere rank those genes and sites to the
top — and do genes with nothing planted stay at the bottom?

That is the whole test. It is a **signal-recovery check on a prioritisation
engine** (the pipeline's permulation p-values are foreground-specificity scores,
never Uniform(0,1) — see `../docs/DECISIONS.md` D-T0-M and the project memory).
There is no KS uniformity gate.

## The experiment (the gate)

One configuration, two runs per archetype:

| | genes | planted |
|---|---|---|
| **power** | ~150 | convergence in ~20% of them |
| **null** | ~150 | none |

Everything else identical: full 237-tip primate tree, `n_pairs = 4` foreground/
background contrast pairs, ~400 sites, `concentration ≈ 2`, production permulation
settings. Modules: `contrast_selection` + `CAAS` + `CT_DISAMBIGUATION` + `SCORING`
(RER / FADE / accumulation / enrichment / VEP off).

### Two archetypes

Both build the foreground the same way — `n_pairs` phylogenetically independent
anchor tips, each paired with its nearest background neighbour — and run through
contrast selection in the loop. They differ only in the trait column:

- **`binary`** — a 0/1 presence/absence code (`--trait_type ordinal`, CLASS 2).
  Anchors = 1, everyone else = 0. Discrete top/bottom candidate path.
- **`rate`** — a `c / n` proportion with `n_pop` / `n_cases` count columns
  (CLASS 1). Anchors get a high rate + tight Jeffreys CI, partners a low rate +
  tight CI, ~30 other species a mid rate + wide CI, so only the anchor/partner
  pairs are CI-separated candidates. This is the shape of the real cancer run
  (`malignant_prevalence`: 3 pairs from a `malignant_count / necropsy_count`
  trait).

The planted molecular signal is placed on the anchor lineages' terminal edges
(each anchor is its own single-tip origin). For a `null` replicate the trait and
pairs are still built — the pipeline runs identically — but nothing is planted.

## Planted mechanisms (D-T0-02)

Per planted gene, ~12 + ~12 sites:

- **`identical_aa`** — near-point-mass equilibrium on one target residue (the
  least-favoured under the site's background profile). Every foreground lineage
  converges to the *same* residue → a strict-CAAS (`US`) positive, and a
  grouped-CAAP (`GS1–GS4`) hit for free.
- **`grouped_caap`** — equilibrium ~uniform over one GS1 physicochemical class
  (the class most disfavoured in the background); a fresh residue from the class
  is seeded per origin, so origins share the *class* but usually differ in
  *residue* → a **GS-only** positive that `US` should not call. This is what
  tests whether the grouped schemes add signal.

Generative model differs from the pipeline's inference model on purpose: WAG
exchangeabilities, per-site Dirichlet profiles (CAT-like), continuous Γ, exact
Gillespie evolution.

## Scoring (`harness score`)

`python -m validation.harness.cli score --run validation/runs/tier0`

Per archetype:

1. **Null vs power separation** (the "type-I" analogue) — pooled: AUC of the
   **detected-CAAS count** (`n_positions`) of planted genes in power replicates
   vs *every* gene in the matched null replicates. This is absolute and
   replicate-comparable; `gene_caas_score` is a within-run percentile rank and
   is reported only as a secondary read (`caas_score AUC`). `separated` = AUC ≥
   0.95 **and** planted p50 > null p95.
2. **Prioritisation** — `precision@k` (of the top-n_planted genes, fraction
   planted), median planted-gene rank percentile, `slice_global25` /
   `slice_global5` membership.
3. **Site recovery** — planted-site detection recall, split by mechanism and by
   scheme: `identical_aa → US`, `grouped_caap → GS` (`US also` = fraction of
   grouped sites that also tripped US).
4. **Contrast recovery** — operative foreground from `traitfile.tab` vs the
   planted pairs (Jaccard + per-pair recovery).

`VERDICT: PASS` when every archetype separates on detected-CAAS count,
`precision@k ≥ 0.8`, and `site_recall ≥ 0.8`.

## Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then one command:
validation/tier0/reproduce.sh                       # -> validation/runs/tier0_gate/
#   stage + run (binary,rate x null,power x primate) + score + figures
# args: reproduce.sh [OUT_DIR] [N_REPLICATES] [N_GENES]   (JOBS=N env var)
```

Or step by step:

```bash
python -m validation.tier0.run_replicates --out validation/runs/tier0 \
    --archetypes binary,rate --sets null,power --trees primate \
    --n-replicates 10 --n-genes 120 --run --jobs 2
python -m validation.harness.cli score --run validation/runs/tier0 --json validation/runs/tier0/score.json
python -m validation.tier0.plots       --run validation/runs/tier0
```

`~40 replicates × ~10–20 min` each; memory-bound (each replicate is a full
Nextflow run). The reference result is committed under `gate_result/` — see its
`SUMMARY.md` to re-score/re-plot without re-running.

## Tests

```bash
python -m pytest validation/tests -q                # harness + tier0 unit tests
Rscript subworkflows/TRAIT_ANALYSIS/local/src/ordinal_trait_test.R   # binary-trait contrast selection
```

## Gate result — `VERDICT: PASS`

Full numbers, figures and reproduce instructions: **`gate_result/SUMMARY.md`**.
Reference run: `primate` × {binary, rate} × {null, power}, 10 reps/cell, 120
genes, 350 sites, 2026-08-31.

| metric | binary | rate |
|---|---|---|
| null-vs-power separation, AUC(detected-CAAS count) | **1.00** | **1.00** |
| — planted p50 / null p95 (null max) | 24 / 1 (2) | 24 / 1 (2) |
| precision@k, by detected-CAAS count | 1.00 | 1.00 |
| precision@k, by `gene_caas_score` | 0.82 | 0.95 |
| site recall / `identical_aa → US` / `grouped_caap → GS` | 0.99 / 1.00 / 0.99 | 0.99 / 1.00 / 0.99 |
| contrast recovery (Jaccard, pairs) | 1.00, all | 1.00, all |

`separation.png`: planted genes carry 22–26 detected CAAS, every non-planted
gene (~170) carries ≤ 2 — zero overlap. **Not worth a bigger run** (the gap is
saturated; cycles don't touch it — it's discovery-driven, not permulation-driven).

Two characterisation findings (not gate blockers):

- **`gene_caas_score` does not separate null from power across replicates**
  (`caas_score AUC ≈ 0.45`) — a within-run percentile rank, size-decorrelated by
  design (`F(max)^n`). A within-study prioritisation device, not a cross-study
  effect size. Detected-CAAS count is the separating quantity.
- **`grouped_caap` sites also trip `US` ~86%.** With `patterns=1,2,3` and a
  conserved disjoint background, a foreground with class-shared / residue-divergent
  substitutions satisfies US pattern 3. A clean GS-only positive needs a
  class-heterogeneous background — deferred. GS recall (0.99) still exceeds US (0.86).

## Status

| piece | state |
|-------|-------|
| `model.py` / `trees.py` / `simulate.py` / `groups.py` / `pheno.py` | done |
| `replicate.py` / `build_project.py` / `run_replicates.py` | done |
| `harness/tier0_adapter.py` + `cli.py score` | done |
| binary-trait pipeline support (contrast selection) | landed on `main` (`1fbdfbf`) |
| **gate run** | **PASS** — Tier 0 closed |

## Characterisation (optional, D-T0-O)

Only after the gate passes, for the paper's curves: the other trees
(`primate_x5`, `mammal`, `star`, `ladder` in `grid.json`), lower
`--planted-fraction`, `--min-divergent-fraction 1.0`. Not run by default.
