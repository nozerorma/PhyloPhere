# Tier 0 fixtures — trees  (DEMOTED, D-DIR-01 — see ../README.md)

Not committed (external data). Fetch from the cluster:

```bash
mkdir -p validation/fixtures/tier0/trees
P=mramon@correfoc-01.s.upf.edu:/homes/users/mramon/scratch/4.Generate_alignments_from_codons/external_phylogenies
scp "$P/primates_233_subst.tree" \
    "$P/primates_233_fossil-cal.tree" \
    "$P/speciesTree_speciesname_pruned.nh" \
    "$P/speciesTree_assemblyname_pruned.nh" \
    "$P/tree_prune_report.tsv" \
    validation/fixtures/tier0/trees/
```

| file | tips | total length (subst/site) | role in Tier 0 |
|------|------|---------------------------|----------------|
| `primates_233_subst.tree` | 233 | 1.30 | ML branch lengths; the regime we actually analyse. Pruned to 50 tips (`grid.json` `primate`), and scaled x5 (`primate_x5`) for a non-degenerate power sweep — the scaled variant is **synthetic**, label it so in any write-up. |
| `speciesTree_speciesname_pruned.nh` | 719 | 29.3 | mammal (TOGA/Zoonomia-style) set; deep enough for real homoplasy. Pruned to 50 tips (`mammal`) — the null **stress** case. Tip labels are species names, matching the mammal alignment set. |
| `primates_233_fossil-cal.tree` | 233 | time-calibrated | not used for simulation (branch lengths are Myr, not substitutions); kept for reference / Tier 3. |
| `speciesTree_assemblyname_pruned.nh` | 719 | 29.3 | same topology, assembly-name tips (matches the TOGA alignment set); kept if a Tier 0 mammal run needs to feed the assembly-name pipeline path. |

Pruning is done at run time by `tier0.trees.prune_depth_preserving` (farthest-point
sampling seeded on the maximum-distance tip pair), so a 50-tip prune keeps
~80–100% of the original tree length rather than collapsing into one dense clade.

Pathological trees (`star`, `ladder`) are built in code, nothing to fetch.
