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

### Trees
- **realistic**: Zoonomia mammal tree pruned to ~40 tips (matches target use)
- **pathological**: 8-tip star topology + a 4-tip ladder — exercises the ASR
  path-score edge cases in the project history (`asr_path_score_axis_redundancy`)

### Simulator — must differ from the inference model
Inference (IQ-TREE / codeml in the pipeline) picks an empirical matrix + Γ.
Simulate under a **different family** so we validate the method, not the
implementation:
- background: `pyvolve` or `INDELible`, WAG or JTT + Γ (pipeline tends to select LG)
- planted positives: site-specific equilibrium-profile switch on the branches
  where the phenotype changes (PCOC's recipe) — plus a second positive class of
  identical-AA convergence (Zou & Zhang style), scored separately, because the
  CAAS module and the FADE/PCOC-adjacent logic target different definitions

### Replicate sets
| set | n | planted signal | foreground | purpose |
|-----|---|----------------|------------|---------|
| `null` | ~1000 | none | random, 2–7 independent transitions | type-I error |
| `power` | ~1000 | 5% of columns, swept | random | power curves |

Sweep grid: transitions ∈ {2,3,5,7} × background branch length ∈ {short, med, long}
× site rate ∈ {slow, med, fast}.

### Trait matrix
Tier 0 is the clean place to cover the full matrix because ground truth is
controlled. For each replicate, `harness.trait_matrix` emits:
- categorical direct (the planted foreground/background, phylo pairs known)
- categorical percentilized (from a continuous phenotype simulated with known λ)
- continuous direct + percentilized (BM phenotype on the tree, `sim.bmPhylo`)
- multi (n_trait=2): one planted-signal trait + one null trait sharing the tree,
  to check the multi-phenotype loop does not leak signal between traits

### Compute
Cheap **if alignment size is capped** (~300–500 codons, ~40 tips). ~2000 short
alignments through CAAS + FADE + RER. Budget ~1 cluster-day. Do **not** run
genome-scale nulls here — that belongs in Tier 2 sizing only (project history
flags it as a recurring hang).

## Status: design only

Next steps:
1. `tier0/simulate.py` — tree loader + pyvolve/INDELible wrapper + planted-shift injector, emitting one alignment dir + truth JSON per replicate
2. `tier0/run.sh` — batch the replicates through the standalone module entrypoints
3. `harness/cli.py calibrate` — collect null p-values across replicates → `NullReport`
