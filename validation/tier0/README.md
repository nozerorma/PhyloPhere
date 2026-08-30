# Tier 0 — simulation: null calibration + planted-positive recovery

**Question:** under a known generative process, does PhyloPhere control type-I
error, and does it rank the planted convergent signal to the top?

**Gate (see `docs/DESIGN.md`):**
- null CAAS `pvalue_boot` (`harness.cli calibrate`, default `--pcol meta_caas_boot`)
  passes KS uniformity (p > 0.05); type-I error near nominal α. `perm_pos_boot`
  is the diagnostic secondary.
- power: planted genes/positions ranked to the top of SCORING; `identical_aa`
  recovered by US, `grouped_caap` recovered by GS (not US)
- contrast selection recovers the planted origins (D-T0-A)

## Design

### Module set (DECISIONS.md D-T0-D, revised)

`contrast_selection` + `CAAS` (discovery, resample, bootstrap, permulation) +
`CT_DISAMBIGUATION` (+ CT_POSTPROC + ASR robustness) + `SCORING`.

**OFF:** RER + FADE (published, author-validated — not what Tier 0 certifies),
CT_ACCUMULATION, ENRICHMENT / STRING / DOMINO, VEP, `reporting` (the exploration
Rmds assume real taxonomy / PhyloPic; contrast_selection runs its own
`dataset_exploration` and emits `trait_stats.csv` regardless).

`ct_disambig_asr_mode = "compute"` — the pipeline gate + compute-on-miss
(commits `8af864b`, `fe5279b`) make that complete now: the perms replay waits
for the live disambiguation, then computes ASR for any null-only gene the
observed run never cached.

### Trees (`grid.json`)

- **`primate`** (the gate): `primates_233_subst.tree` used **whole** (237 tips =
  233 primates + 4 outgroups, total length 1.30 subst/site) — exactly what
  production runs on. A 50-tip prune left contrast selection unable to form ≥
  `min_contrasts` Dunn-independent pairs (D-T0-01).
- **characterisation-only**: `primate_x5` (branch lengths ×5, synthetic power
  curve), `mammal` (`speciesTree_speciesname_pruned.nh` pruned to 60, deep-homoplasy
  null stress), `star`/`ladder` (12/16 tips, built in `trees.py`, the
  degenerate-MRCA `asr_path_score` fixtures).

`sample_foreground` picks origins farthest-point (pairwise patristic separation)
so the planted clades are Dunn-independent, and auto-excludes the 4 outgroups
(terminal branch > 5× the 95th pct) from ever being foreground.

### Simulator — differs from the inference model

The pipeline's ASR defaults to **LG, one profile per site** (`ct_disambig_asr_model
= "lg"`, codeml model 3). Tier 0 simulates under:
- **WAG exchangeabilities** (`data/wag.dat`, vendored) — a different matrix
- **site-specific Dirichlet profiles** ~ Dirichlet(`concentration` · π_WAG) —
  CAT-like heterogeneity the single-profile inference does not model; this is
  what makes non-adaptive ("background") convergence arise
- continuous Γ across sites (`gamma_alpha`)
- exact **Gillespie** evolution along each branch — no matrix exponential;
  branch-specific rate matrices (planted positives) are free

Two planted-positive mechanisms (D-T0-02), 12 sites each per planted gene:
- `identical_aa` — near-point-mass equilibrium on one target residue (the
  least-favoured under the site background) → every fg lineage converges to the
  *same* residue. A **US and GS1-GS4** positive.
- `grouped_caap` — equilibrium ~uniform over one GS class (default GS1, the
  class most disfavoured in the site background); a fresh residue from the class
  is seeded per origin, so origins share the *class* but usually differ in
  *residue*. A **GS-only** positive — US should not call it. This is what
  exercises the grouped schemes. GS tables vendored in `groups.py`.
- `profile_shift` is **stripped** (dormant `SimConfig` hook).

### Phenotype archetypes (`pheno.py`, DECISIONS.md D-T0-C)

- **`echo`** — exact 0/1 presence/absence. RER (when on) auto-detects binary;
  contrast selection's discrete path categorises trivially; the CAAS permulation
  works because of the median-guard fix (`8af864b`). `top`/`bottom_quantile` set
  from the fg fraction so the cut lands on the 1s / 0s.
- **`bodysize`** — a Brownian-motion continuous trait on the tree; foreground =
  its top-quantile tips (emergent, not pre-sampled).

### Replicate structure (`replicate.py`)

One replicate = one (tree, phenotype) draw. Within it, `n_genes` genes are
simulated on the same tree + phenotype; `frac_planted_genes` carry planted
convergent substitutions, the rest are pure null. This gives SCORING a gene
population, a meaningful `percent_rank(pvalue_boot)`, and many null p-values.

Output per replicate (all synthetic support files the pipeline needs but real
data would carry): `align/gNNNN.fasta` (gene ids have **no `_`** — see the
`disambiguation_main.py` note), `my_traits.tsv` (+ `family` column), `tree.nwk`,
`ali_sp_names.txt`, `taxid.tsv` (fake taxids — production renames tips to numeric
ids before PAML), `gene_ensembl.tsv` (fake coords), `truth.json`.

### Scope (D-T0-O)

**The gate** — what certifies the pipeline — is a single configuration:
`primate` tree × {`echo`, `bodysize`} × {`null`, `power`} at
`min_divergent_fraction 0.5`, `planted_fraction 0.25`, 20 replicates each. That
is `run_replicates.py`'s default.

**Characterisation** (optional, only if the gate passes and we want curves for
the paper): the other trees (`primate_x5`, `mammal`, `star`, `ladder`) and the
extra axis values (`--divergent-fractions 0.5,1.0`, `--planted-fractions 0.1,0.25`).
Not run by default.

### Compute

Simulation is ~0.02 s/gene. The cost is the pipeline: CAAS discovery + codeml ASR
(main disambiguation + the perms replay's null-only tail) + the permulation pool.
On the full 237-tip primate tree, 15 genes / 200 sites ≈ 4.5 min per replicate
locally. `perm_pool_size` / `caas_full_perms` scale the null.

## Status

| piece | state |
|-------|-------|
| `model.py` / `trees.py` / `simulate.py` / `pheno.py` / `replicate.py` | **done**, 14 tests |
| `build_project.py` — Tier 0 `ProjectConfig` → `gui.generation.render` | **done** |
| `run_replicates.py` — stage + `--run` | **done** |
| **end-to-end pipeline run** | **GREEN** (2026-08-30). `echo_power`, full primate tree, 15 genes, `asr_mode=compute`. 28 processes, 4.5 min, exit 0. `gene_scores.tsv` = exactly the planted genes ranked by `gene_caas_score`; top `position_scores` all planted sites. Permulation null covers every gene with per-cycle CAAS. |
| `harness/cli.py score` / `calibrate` — adapter (`harness/tier0_adapter.py`) | **done**, tests |
| contrast recovery (D-T0-A) + grouped-CAAP scheme recovery (D-T0-N) in `score` | **done** |
| **gate run** (`primate` × echo,bodysize × null,power, 20 reps) | **running** |

## Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then stage (gate config
# is the default) + run:
python -m validation.tier0.run_replicates --out validation/runs/tier0 \
    --archetypes echo,bodysize --sets null,power --trees primate \
    --n-replicates 20 --n-genes 40 --run --jobs 4
# then score:
python -m validation.harness.cli calibrate --run validation/runs/tier0   # KS gate
python -m validation.harness.cli score     --run validation/runs/tier0 --json validation/runs/tier0/score.json
```
