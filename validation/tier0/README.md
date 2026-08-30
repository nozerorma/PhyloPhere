# Tier 0 — simulation: null calibration + planted-positive recovery

**Question:** under a known generative process, does PhyloPhere control type-I
error, and does it rank the planted convergent signal to the top?

**Gate (see `docs/DESIGN.md`):**
- null CAAS `p.perm` / `pvalue_boot` and the FCS `p.perm` pass KS uniformity
  (p > 0.05); observed type-I error within Monte Carlo error of nominal α
- power: planted genes/positions in the top decile of the SCORING ranking
- contrast-selection recovers the planted origins (D-T0-A metric)
- zero crashes on the star / ladder pathological trees

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

- **`primate`**: `primates_233_subst.tree` used **whole** (237 tips = 233
  primates + 4 outgroups, total length 1.30 subst/site) — exactly what production
  runs on. A 50-tip prune left contrast selection unable to form ≥ `min_contrasts`
  Dunn-independent pairs (D-T0-01).
- **`primate_x5`**: same tree, branch lengths ×5 — a non-degenerate power curve.
  **Synthetic**; label it so.
- **`mammal`**: `speciesTree_speciesname_pruned.nh` pruned to 60 — deep enough
  for real homoplasy, the null **stress** case.
- **pathological**: 12-tip `star` + 16-tip `ladder`, built in `trees.py` — the
  star is the empty-private-segment / degenerate-MRCA fixture from
  `asr_path_score_axis_redundancy`.

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

Two planted-positive mechanisms, recorded separately in `truth.json`:
- `profile_shift` — foreground edges use a *different* Dirichlet profile (PCOC's
  preference-shift model)
- `identical_aa` — foreground edges have a near-point-mass equilibrium on one
  target residue → every fg lineage converges to the *same* residue (Zou & Zhang;
  the strict-CAAS positive). Seeded at each origin edge.

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

### Sub-set axes (`run_replicates.py`)

`archetype` {echo, bodysize} × `set` {null, power} × `tree` × `min_divergent_fraction`
{0.5 production, 1.0 strict — D-T0-G} × `planted_fraction` {0.1, 0.25 — D-T0-H,
power only}.

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
| `harness/cli.py score` / `calibrate` — pipeline-output adapter | **in progress** |
| contrast-selection recovery metric (D-T0-A) | **TODO** |
| full `echo`/`bodysize` × `null`/`power` matrix + scale | **TODO** |

## Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then stage:
python -m validation.tier0.run_replicates \
    --out validation/runs/tier0 \
    --archetypes echo,bodysize --sets null,power --trees primate \
    --n-replicates 20 --n-genes 40
# review runs/tier0/**/run.sh, then execute:
python -m validation.tier0.run_replicates ... --run --jobs 4
# then score:
python -m validation.harness.cli calibrate --run validation/runs/tier0
python -m validation.harness.cli score    --run validation/runs/tier0
```
