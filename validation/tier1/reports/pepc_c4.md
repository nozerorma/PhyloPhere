# Report — PEPC / C4 photosynthesis (Cyperaceae)

## a) Data source

79 Cyperaceae PEPC amino-acid sequences, projected onto maize PEPC1 coordinates
(UniProt P04711, 970 aa). Originally Besnard et al. 2009 (*MBE* 26:1909),
redistributed via the PCOC (Rey et al. 2018, *MBE* 35:2296) and ConDor (Morel
et al. 2024, *GBE*, doi:10.1093/gbe/evae040) test-data repositories. Trait: `c4`,
the **genotypic** C4 annotation (A780S determinant), 23 C4 / 55 C3 tips, ≥5
independent C4 origins. Truth set: 11 curated sites (`../../truthsets/tier1/pepc_c4.sites.tsv`)
across four confidence tiers (mutagenesis, selection, parallel, weak). Full
provenance: [`../input/pepc/README.md`](../input/pepc/README.md).

## b) Scope

Single-gene, categorical-trait recovery check: does the pipeline, run unchanged,
recover mutagenesis- and selection-supported convergent sites and rank them
sensibly? No null-calibration or genome-wide inflation analysis — a single gene
gives no genome-wide background to calibrate against, and RER's cross-gene
comparison has essentially no power here for the same reason (see
[`../input/pepc/README.md`](../input/pepc/README.md)).

## c) Procedure

**Legacy run** (`../../runs/tier1/pepc/c4_mdf0.5/`, pre-GUI harness scripts, since
retired): CAAS/CAAP discovery+resample → CT_DISAMBIGUATION (ASR + postproc filter)
→ FADE → SCORING. RER/accumulation/enrichment/VEP off. `min_divergent_fraction=0.5`,
`max_fop=100`, seed 1998, local.

**Current template** (`../../../gui/templates/tier1_pepc.json`): same CAAS/CT_DISAMBIGUATION/FADE
chain, **RER added** (single-gene, low power by construction — wired for
completeness against the "run CAAS+RER+FADE" scope, not because a meaningful
result is expected), accumulation/VEP/enrichment still off, `local_lowspec`
resource preset. Not yet re-run under this template — §d's results below are the
legacy CAAS+FADE-only run; re-running with RER on is the immediate next step,
and CAAS/FADE results should not change materially since their config is
unchanged.

## d) Results and conclusions

**6/11 sites recovered** (legacy run; matches ConDor's own comparator run
site-for-site):

| tier | recall | detail |
|------|--------|--------|
| mutagenesis | 2/2 | A780S (rank-pct 0.25), H665N (rank-pct **0.04**) |
| selection | 2/3 | P540T (0.13), E572Q (0.09); M749T missed (4/23 C4, minority convergence) |
| parallel | 1/2 | S761A (0.96); F733V missed (2/23 C4, below the CAAS exclusivity bar) |
| weak | 1/4 | one residue-non-match |

29 scored positions total, **none gate-significant** at the hypergeometric+bootstrap
threshold — expected at this scale (79 tips, one gene), the gate is tuned for
genome-wide multi-gene runs.

**H665N is the clearest cautionary result**, not a clean win: ASR-clean (0.75)
and a real mutagenesis-confirmed convergent substitution, but a His→Asn change
of that strength recurs in ~47% of alignment bootstraps
(`recovery_boot`≈0.472 → `phen_score`≈0.064 → `CAAS_score`≈0.048, rank-pct 0.04).
The score is doing its job — it is *not* a false negative in the sense of a
missed signal, it is a real signal downweighted because it is not phylogenetically
exclusive to this alignment's bootstrap resamples. This is a property of the
scoring mechanism (bootstrap-based specificity discount), not a heuristic guess
about the gene — worth stating plainly since it is the kind of result that looks
like a bug until the mechanism is traced.

The two misses (M749T, F733V) are both minority-lineage convergence (4/23 and
2/23 C4 tips respectively) below CAAS's exclusivity requirements — a known,
documented limitation of the exclusivity-gated CAAS definition itself, not a
pipeline defect specific to this run.
