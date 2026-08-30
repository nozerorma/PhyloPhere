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

1. **Prioritisation** — planted genes in `slice_global1` / `slice_global5`;
   planted-gene ranks; `precision@k`.
2. **Null vs power separation** (the "type-I" analogue) — pooled: AUC of
   planted-gene `gene_caas_score` in power replicates against *every* gene's
   score in the matched null replicates. `separated` = AUC ≥ 0.9 **and** planted
   p50 > null p95. If a no-signal replicate reaches the same scores, the ranking
   cannot tell signal from noise.
3. **Site recovery** — planted-site detection recall, split by mechanism and by
   scheme: `identical_aa → US`, `grouped_caap → GS` (with `US` leakage, want ~0).
4. **Contrast recovery** — operative foreground from `traitfile.tab` vs the
   planted pairs (Jaccard + per-pair recovery).

`VERDICT: PASS` when every archetype separates and planted genes land in
`slice_global5` ≥ 80% of the time.

## Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then:
python -m validation.tier0.run_replicates --out validation/runs/tier0 \
    --archetypes binary,rate --sets null,power --trees primate \
    --n-replicates 15 --n-genes 150 --run --jobs 2
python -m validation.harness.cli score --run validation/runs/tier0 --json validation/runs/tier0/score.json
```

`~60 replicates × ~5–15 min` each; use `--jobs` to parallelise (memory-bound —
each replicate is a full Nextflow run).

## Status

| piece | state |
|-------|-------|
| `model.py` / `trees.py` / `simulate.py` / `groups.py` | done |
| `pheno.py` — `make_paired_foreground` (binary + rate) | done |
| `replicate.py` / `build_project.py` / `run_replicates.py` | done |
| `harness/tier0_adapter.py` + `cli.py score` — prioritisation + separation | done |
| binary-trait pipeline support (contrast selection) | landed on `main` (`1fbdfbf`) |
| **gate run** (`primate` × {binary, rate} × {null, power}) | **pending** |

## Characterisation (optional, D-T0-O)

Only after the gate passes, for the paper's curves: the other trees
(`primate_x5`, `mammal`, `star`, `ladder` in `grid.json`), lower
`--planted-fraction`, `--min-divergent-fraction 1.0`. Not run by default.
