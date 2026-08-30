# Tier 0 — simulation: null calibration + planted-positive recovery

**Question:** under a known generative process, does PhyloPhere control type-I
error, and what is its power as a function of transition count, branch length,
site rate, and tree size?

**Gate (see `docs/DESIGN.md`):**
- null CAAS `p.perm` / RER `p.perm` / FADE p-values pass KS uniformity (p > 0.05)
- observed type-I error within Monte Carlo error of nominal α
- power curves monotone in signal strength and transition count; convergence >
  divergent > mixed separation preserved
- zero crashes on the star-topology / empty-private-segment fixture

## Design

### Trees (`grid.json`)
- **`primate`**: `primates_233_subst.tree` (233 tips, total length 1.30 subst/site)
  pruned to 50 by farthest-point sampling. The regime we actually analyse; the
  native tree is shallow so the null here is the *easy* case.
- **`primate_x5`**: same prune, branch lengths ×5 — a non-degenerate power curve.
  **Synthetic**; label it so in any write-up.
- **`mammal`**: `speciesTree_speciesname_pruned.nh` (719 tips, total length 29.3)
  pruned to 50. Deep enough for real homoplasy — the null **stress** case.
- **pathological**: 12-tip `star` + 16-tip `ladder`, built in `trees.py`. The
  star is the empty-private-segment / degenerate-MRCA fixture from
  `asr_path_score_axis_redundancy`.

### Simulator — differs from the inference model
The pipeline's ASR / FADE default to **LG, one profile per site** (`ct_disambig_asr_model
= "lg"`; `run_hyphy_fade_batch.sh` `model="LG"`). Tier 0 simulates under:
- **WAG exchangeabilities** (`data/wag.dat`, vendored) — a different matrix, and
- **site-specific Dirichlet profiles** ~ Dirichlet(`concentration` · π_WAG) — CAT /
  C-series-like heterogeneity the single-profile inference does not model. This is
  what makes non-adaptive ("background") convergence actually arise, the way it
  does in real data.
- continuous Γ across sites (`gamma_alpha`).
- exact **Gillespie** evolution along each branch — no matrix exponential, and
  branch-specific rate matrices (for the planted positives) are free.

Two planted-positive mechanisms, recorded separately in `truth.json` so recall is
reported per mechanism:
- `profile_shift` — foreground edges use a *different* Dirichlet profile (PCOC's
  convergent-preference-shift model).
- `identical_aa` — foreground edges have a near-point-mass equilibrium on one
  target residue, so every foreground lineage converges to the *same* amino acid
  (Zou & Zhang style; the strict-CAAS positive). Seeded at each origin edge.

Implementation: `model.py` (matrix + profiles + Γ), `trees.py` (load / depth-
preserving prune / star+ladder / foreground sampling), `simulate.py` (evolve +
plant + write).

### Replicate sets
| set | n | planted signal | foreground | purpose |
|-----|---|----------------|------------|---------|
| `null` | ~1000 | none | random, 2–7 independent transitions | type-I error |
| `power` | ~1000 | 5% of columns, swept | random | power curves |

Sweep grid: transitions ∈ {2,3,5,7} × background branch length ∈ {short, med, long}
× site rate ∈ {slow, med, fast}.

### Trait matrix
`make_replicates.py` calls `harness.trait_matrix` per replicate. The simulated
phenotype is categorical, so it currently emits **categorical/direct/single**
(`sim_fg.cfg`) — the cell CAAS/FADE consume. The remaining cells (continuous
direct/percentilized from a simulated BM phenotype with known λ; multi n_trait=2
with one signal + one null trait to check the multi-phenotype loop does not leak)
need a continuous `PhenotypeSpec` added to the emitted `DatasetSpec` — **TODO**.

### Compute
Simulation is ~0.02 s/replicate; 2000 replicates ≈ 1 min. The cost is entirely in
running PhyloPhere's CAAS/FADE/RER on them. Cap alignment size (~400 codons, 50
tips). Do **not** run genome-scale nulls here — Tier 2 sizing only (project
history flags it as a recurring hang).

## Status

| piece | state |
|-------|-------|
| `model.py` / `trees.py` / `simulate.py` | **done**, tested (`validation/tests/test_tier0.py`) |
| `make_replicates.py` + `grid.json` | **done** — emits `null` / `power` sets + `manifest.jsonl` |
| continuous / multi trait cells | TODO (needs simulated BM phenotype) |
| `tier0/run.sh` — batch replicates through the standalone CAAS/FADE/RER entrypoints | **TODO** (next; needs the module CLI invocations — write against one real run) |
| `harness/cli.py calibrate` — collect null p-values → `NullReport` | **TODO** (shares the output adapter with `score`) |

### Run

```bash
# fetch trees once (validation/fixtures/tier0/README.md), then:
python -m validation.tier0.make_replicates \
    --config validation/tier0/grid.json \
    --out validation/runs/tier0 --set all
# smoke test: add --limit 6
```
