# Tier 1 PEPC — Results Report

Truth set: `validation/truthsets/tier1/pepc_c4.sites.tsv`, 10 positions in maize PEPC1 (CAA33317) numbering, 1:1 with alignment column (no reference row). Two methods evaluated independently: OC's FUBAR (site-level dN/dS, HyPhy) and PhyloPhere's own CAAS/CT_DISAMBIGUATION → SCORING pipeline (contrast-based convergent-substitution scoring).

## Numbering note

PhyloPhere's `scoring/position_scores.tsv` uses `Position` values that are consistently **one less** than the truth set's maize-PEPC1 numbering (verified by matching `derived_residues` against each truth site's ref/alt amino acids at `truth_position - 1`, e.g. truth position 780 A→S = table position 779, `derived_residues = S/A`). All PhyloPhere positions below are reported in maize-PEPC1 numbering (`table Position + 1`) for direct comparability with the truth set and with FUBAR.

## Method 1 — OC / FUBAR (HyPhy, posterior P[β>α] per site, 970 codons scanned)

Significance threshold used elsewhere in this validation for FUBAR: posterior > 0.9.

| position | ref>alt | tier | P[pos.sel] | rank /970 | sig P>0.9 |
|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | 0.0138 | 650 | no |
| 665 | H→N | mutagenesis | 0.0000 | 944 | no |
| 540 | P→T | selection | 0.0003 | 775 | no |
| 572 | E→Q | selection | 0.0329 | 608 | no |
| 733 | F→V | parallel | 0.0261 | 621 | no |
| 761 | S→A | parallel | 0.0398 | 604 | no |
| 749 | L→T | weak | 0.2727 | 33 | no |
| 505 | F→L | weak | 0.0000 | 825 | no |
| 573 | A→N | weak | 0.0038 | 695 | no |
| 731 | I→V | weak | 0.0004 | 773 | no |

**Overlap: 0/10.** No truth-set position clears P>0.9. Best-ranked truth position is 749 (weak tier) at 33/970; the median rank across all 10 truth positions is ≈622/970 — no systematic enrichment near the top of FUBAR's own distribution, including for the two mutagenesis-tier (strongest-evidence) sites (650/970 and 944/970).

**Beyond the truth set:** FUBAR flags exactly 2 sites genome-wide at P>0.9 — codon 630 (P=0.9412) and codon 474 (P=0.9236). Neither is in the truth set.

## Method 2 — PhyloPhere (CAAS → CT_DISAMBIGUATION → SCORING, `p.emp_adj`)

CAAS/CT_DISAMBIGUATION does not score all 970 alignment columns — contrast-based filtering (gap/missingness thresholds, `min_divergent_fraction`, cluster filters) reduces the candidate set to 70 unique positions before SCORING assigns `p.emp_adj`. Significance threshold: `p.emp_adj < 0.1` (`scoring_p_emp_thr` default, `conf/scoring.config:49`).

| position (maize) | ref>alt | tier | in candidate set (n=70) | p.emp_adj | rank | sig p.emp_adj<0.1 |
|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | yes | 0.01658 | 1/70 | **yes** |
| 665 | H→N | mutagenesis | yes | 0.01658 | 3/70 | **yes** |
| 540 | P→T | selection | yes | 0.01658 | 2/70 | **yes** |
| 572 | E→Q | selection | yes | 0.04523 | 7/70 | **yes** |
| 733 | F→V | parallel | **no (filtered out upstream)** | — | — | — |
| 761 | S→A | parallel | yes | 0.11549 | 19/70 | no |
| 749 | L→T | weak | yes | 0.11549 | 21/70 | no |
| 505 | F→L | weak | yes | 0.08687 | 15/70 | **yes** |
| 573 | A→N | weak | yes | 0.11549 | 22/70 | no |
| 731 | I→V | weak | yes | 0.03554 | 5/70 | **yes** |

**Overlap: 6/10 significant** (780, 665, 540, 572, 505, 731); **1/10 absent from the candidate set entirely** (733, parallel tier — never reaches SCORING, filtered upstream in CT_DISAMBIGUATION/postproc, not just non-significant); **3/10 present but not significant** (761, 749, 573 — all rank in the bottom third of the 70-position candidate set, 19–22/70). Both mutagenesis-tier sites (780, 665) rank #1 and #3 of 70, the strongest possible showing for the highest-confidence truth tier.

**Beyond the truth set:** 11 additional positions clear `p.emp_adj < 0.1`:

| position (maize) | derived_residues | side | CAAS_score | p.emp_adj |
|---|---|---|---|---|
| 588 | IKL/R | top | 0.700 | 0.02764 |
| 751 | F/Y | top | 0.595 | 0.04523 |
| 862 | A/AS | bottom | 0.499 | 0.04523 |
| 852 | E/D | top | 0.636 | 0.05365 |
| 611 | L/F | top | 0.619 | 0.05365 |
| 620 | AC/S | top | 0.401 | 0.05365 |
| 579 | ET/A | top | 0.394 | 0.05365 |
| 627 | I/V | top | 0.436 | 0.07370 |
| 839 | K/G | top | 0.400 | 0.08687 |
| 584 | MST/I | top | 0.281 | 0.08687 |
| 625 | AV/IV | top | 0.230 | 0.09373 |

Position 751 sits immediately adjacent to truth position 749/750 (L→T, weak tier) but is a distinct candidate (F→Y) not annotated in the truth set.

## Cross-method comparison

| | FUBAR (OC) | PhyloPhere (CAAS/SCORING) |
|---|---|---|
| Truth positions recovered as significant | 0/10 | 6/10 |
| Truth positions scored but non-significant | 10/10 | 3/10 |
| Truth positions never reaching the method's own candidate/output set | 0/10 | 1/10 (733) |
| Both mutagenesis-tier sites (780, 665) significant | no | yes (rank 1, 3) |
| Novel significant sites beyond truth set | 2 | 11 |

## Caveat: the trait itself is defined by position 780 (circularity)

This fixture's `c4` trait (`validation/tier1/input/pepc/my_traits.tsv`) is
ConDor's **"genotypic"** annotation: `c4 = 1` iff the tip's ppc-1 sequence
carries the A780S substitution at maize position 780. Checked directly
against every tip's actual residue: **22/23 c4=1 tips carry S at 780, 53/54
c4=0 tips carry A — a 100% match among the 76/78 tips where the residue is
observable** (the 2 exceptions are gaps, not mismatches). This is not
approximate correlation; by construction it is very close to identity.

Besnard et al. (2009) themselves define "C4 ppc" lineages in their own gene
tree by this residue (*"Branches leading to genes encoding C4 PEPC (with a
serine at position 780) are in bold,"* Fig. 3 caption). Morel et al. (2024)
state explicitly, for this exact 78-tip dataset, that the genotypic
annotation is *"grounded in the fact that... the A780S mutation... has been
experimentally demonstrated to be a major determinant of C4-specific
characteristics"* and that they *"predicted the metabolism associated with
the sedge PEPC sequences according to the presence or absence of the A780S
mutation"* — i.e. the trait literally is the residue, by the source paper's
own account, not an independent phenotype that happens to correlate with it.

Morel et al. (2024) also demonstrate the empirical cost of this: their PCOC
benchmark recovers 7/11 true positives under the genotypic annotation but
**0/11 under an independent, phenotype-based annotation** (Bruhl & Wilson
2007's anatomical/physiological C3/C4 survey) on the same underlying data.

**Consequence for this report**: recovering position 780 (both methods'
top-ranked or near-top-ranked hit) is close to guaranteed by construction
under the genotypic trait and should not be read as a genuine test of either
method's power. The other 9 truth-set positions are unaffected — none of
them are the trait-definition site — and remain valid recovery tests.
A non-circular counterpart fixture using Bruhl & Wilson's independent
phenotypic classification (dropping the trait's dependency on 780 entirely)
is built at `validation/tier1/input/pepc_phenotypic/` for a follow-up run.

## Caveat: some "independent" contrast pairs are not phylogenetically adjacent

CAAS/CT_DISAMBIGUATION's contrast-selection step draws its foreground/
background pairs from PSS-selected extreme-divergence species, not from
sister taxa. Checked directly on `contrast_hypotheses_pairs.tsv` (135
distinct pairs across the 100 hypotheses used in this run): **no pair is a
literal sister-cherry** (minimum MRCA span across all 400 pair-instances is 4
tips), **31/135 (23%) pair members are already different genera**, and
**63/135 (47%) have a genus other than either member's nested somewhere
between them** in the tree. This reproduces on `tree.nwk` and independently
on Besnard et al. (2009)'s own original PhyML tree
(`besnard2009/pepc.phyml_tree.txt`) — not an artifact of this fixture's tree
construction.

Most of this traces to one documented case: *Cyperus* as classically
circumscribed is paraphyletic, with several segregate genera (*Kyllinga*,
*Pycreus*, *Remirea*, *Volkiella*) nested within it — a real feature of
Cyperaceae systematics, not a fixture error. Separately, two of the fixture's
multi-accession species (*Eleocharis baldwinii*, *E. vivipara*) have their
C4-associated ppc-1 accessions cluster with *each other* rather than with
their own species' C3 accessions; Besnard et al. (2009) attribute this to a
documented horizontal-transfer/hybridization event between the two species
(*"These two unrelated taxa seem to have acquired their C4 ppc from the same
source... through horizontal gene transfer or hybridization,"* p. 1916),
independent of the general genus-paraphyly pattern above.

**Consequence**: `min_contrasts=3` (of 4 pairs per hypothesis) is intended to
read as "detected in ≥3 independent lineages." With zero true sister-pairs in
use and a documented non-monophyletic backbone, some of that count may draw
from within the same local radiation rather than genuinely separate origins.
Worth caveating wherever `min_contrasts` is cited as a replication count.

## References

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009.
  Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple
  appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919.
  doi:10.1093/molbev/msp103.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4
  photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148.
  doi:10.5642/aliso.20072301.11.
- Morel B, Schade P, Le Douarin M, Lartillot N, Villoutreix R, Chevalier T,
  Bosseur F, Gascuel O, Guindon S. 2024. ConDor: accurate detection of
  convergent evolution with amino acid substitutions. Genome Biol Evol
  16(4):evae040. doi:10.1093/gbe/evae040.
- Rey C, Guéguen L, Sémon M, Boussau B. 2018. Accurate detection of
  convergent amino-acid evolution with PCOC. Mol Biol Evol 35(9):2296–2306.
  doi:10.1093/molbev/msy114.
