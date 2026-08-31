# Tier 0 gate — reference result

Run: `primate` tree × {binary, rate} × {null, power}, **10 replicates/cell**,
120 genes, 350 sites, `n_pairs = 4`, `concentration = 2`, production permulation
settings (`perm_pool_size = 1000`, `caas_full_perms = 200`). 40/40 replicates
exit 0. Date: 2026-08-31.

## `VERDICT: PASS`

| metric | binary | rate |
|---|---|---|
| **null vs power, AUC(detected-CAAS count)** | **1.000** | **1.000** |
| — planted p50 / null p95 (null max) | 24 / 1 (2) | 24 / 1 (2) |
| precision@k, by detected-CAAS count | **1.00** | **1.00** |
| precision@k, by `gene_caas_score` | 0.82 | 0.95 |
| site recall (any scheme) | 0.994 | 0.993 |
| `identical_aa → US` | 0.998 | 0.998 |
| `grouped_caap → GS` | 0.989 | 0.987 |
| `grouped_caap` also trips US | 0.86 | 0.86 |
| contrast recovery (Jaccard, pairs) | 1.00, all | 1.00, all |

`separation.png` is the load-bearing figure: planted genes carry 22–26 detected
CAAS, every non-planted gene (null **or** power replicate, ~170 genes) carries
≤ 2. Zero overlap.

## Is a bigger run worth it?

**No.** The separation is a clean gap with nothing between 3 and 21, over 20
replicates and ~170 non-planted genes. More genes add more of each cluster
without shrinking the gap. Cycles don't touch it at all — the separation is
driven by CAAS discovery + disambiguation, not the permulation (which only
feeds `gene_caas_score`, and that doesn't separate across replicates by design).
The only thing a bigger run tightens is the CI on `precision@k(gene_caas_score)`
≈ 0.82 vs 0.95 — a characterisation of the ranking, not a gate question.

## Two characterisation findings (not gate blockers)

1. **`gene_caas_score` does not separate null from power across replicates**
   (`caas_score AUC ≈ 0.45`). It is a within-run percentile rank whose
   size-adjustment (`F(max)^n`) deliberately decorrelates it from position
   count → a within-study prioritisation device, not a cross-study effect size.
2. **`grouped_caap` sites also trip US ~86%.** With `patterns = 1,2,3` and a
   conserved, disjoint background, a foreground that shares a physicochemical
   class but varies in residue satisfies US pattern 3 (fg variable, bg
   conserved). A clean GS-only positive would need a class-heterogeneous
   background — a mechanism refinement. GS recall (0.99) still exceeds US (0.86).

## Files

| file | what |
|---|---|
| `score.json` | full `harness score` output for this run |
| `figures/*.png` | the four gate figures |
| `manifest.jsonl` | per-replicate config (seed, tree, archetype, set) |
| `replicates.tar.gz` | `truth.json` + the five files the adapter reads per replicate (`scoring/{gene_scores,position_scores}.tsv`, `scoring/gene_lists/slice_global{5,25}.tsv`, `data_exploration/2.CT/1.Traitfiles/traitfile.tab`) — ~500 KB instead of the ~5 GB of full pipeline output |

## Reproduce

```bash
# full re-run (needs the primate tree fixture + phylophere env, ~2-4 h):
validation/tier0/reproduce.sh              # -> validation/runs/tier0_gate/

# or re-score / re-plot this reference without re-running the pipeline:
mkdir /tmp/gr && tar xzf validation/tier0/gate_result/replicates.tar.gz -C /tmp/gr
python -m validation.harness.cli score --run /tmp/gr
python -m validation.tier0.plots        --run /tmp/gr --out /tmp/gr/figures
```
