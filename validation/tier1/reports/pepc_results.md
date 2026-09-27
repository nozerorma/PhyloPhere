# Tier 1 PEPC: Results Report

Truth set: `validation/truthsets/tier1/pepc_c4.sites.tsv`, 10 positions in maize PEPC1 (CAA33317) numbering, 1:1 with the 970-column fixture alignment (`input/pepc/align/PEPC.fasta`, no reference row). Two methods: OC's FUBAR (site-level dN/dS, trait-independent) and PhyloPhere (CAAS → CT_DISAMBIGUATION → FADE → SCORING), run under two trait definitions on the same alignment and tree.

| PhyloPhere run | trait | definition | tips (in-group) | C4 / C3 |
|---|---|---|---|---|
| `results/c4_complete/` | `c4` | genotypic: tip's ppc-1 carries A780S (Besnard et al. 2009; Morel et al. 2024) | 77 | 23 / 54 |
| `results/c4_phenotypic_complete/` | `c4_phenotypic` | phenotypic: species-level C3/C4 call from Bruhl & Wilson (2007); 7 C3/C4-intermediate *Eleocharis* accessions pruned | 70 | 20 / 50 |

Both traits are columns of `input/pepc/my_traits.tsv`; sourcing of the phenotypic column is documented in `input/.pepc_phenotypic/README.md`. The outgroup *Chrysitrix dodii* is in the alignment and tree but carries no trait value in either run. Cross-trait interpretation is in `pepc_genotypic_vs_phenotypic.md`.

## Numbering

`scoring/position_scores.tsv` `Position` is the 0-based alignment column, i.e. maize position minus 1, in **both** runs. Verified per truth position by recomputing C4/C3 residue counts at alignment column *p* (1-based) from the raw FASTA under each run's tip set and trait, and matching them to `top_species_residues` / `bottom_species_residues` at `Position = p − 1`: exact match for all 9 truth positions present in each run (e.g. maize 780: genotypic `S:22` / `A:53`, phenotypic `S:16,A:2` / `A:49,S:1`). All positions below are in maize numbering.

## `position_scores.tsv` schema

Columns: `Gene, Position, n_schemes, scheme_set, n_hypotheses, participating_hypotheses, top_species_residues, bottom_species_residues, n_top_species, n_bottom_species, CAAS_score, side, caas, p.emp, p.emp_adj`. One row per (Position, side); a position detected on both sides has two rows with identical `p.emp`/`p.emp_adj` (verified: ≤ 1 distinct value per position in both runs).

`p.emp` is the pooled "detects AND exceeds" permulation p (`subworkflows/SCORING/local/src/scoring_compute.R` §2f-ter): `(k_emp + 1)/(N + 1)`, where `k_emp` counts null cycles that re-detect the position on any side with max-over-sides CAAS ≥ the observed max. `p.emp_adj` is BH over position-rows. Significance: `p.emp_adj < 0.1` (`scoring_p_emp_thr`). Raw `p.emp` is reported alongside throughout. Raw p-values in this candidate set are not uniform under the null (they are capped by each position's null detection rate), so raw-p counts overstate evidence; see `pepc_genotypic_vs_phenotypic.md` §6 for the choice of correction family.

Note for ad hoc pandas reads: the `caas` value `N/A` (maize 573, both runs) is parsed as missing under pandas defaults; read with `keep_default_na=False`.

| | genotypic | phenotypic |
|---|---|---|
| position-rows | 72 | 66 |
| unique positions (candidate set) | 61 | 51 |
| null cycles `N` | 999 | 1000 |
| minimum attainable `p.emp` | 0.001 | 0.000999 |
| positions with raw `p.emp < 0.05` | 9 | 10 |
| positions with raw `p.emp < 0.1` | 18 | 16 |
| positions with `p.emp_adj < 0.1` | 3 | 1 |

## Method 1: OC / FUBAR

HyPhy FUBAR on 78 sequences × 970 codons (`input/pepc/oc_run/PSEL/PEPC.FUBAR.json`, values re-extracted). Trait-independent, so a single run applies to both trait definitions. Threshold: posterior P[β > α] > 0.9.

| position | ref>alt | tier | P[β>α] | rank /970 | sig |
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

**0/10.** FUBAR flags 2 sites genome-wide (codon 630, P = 0.941; codon 474, P = 0.924), neither in the truth set.

## Method 2: PhyloPhere

Rank = competition rank of `p.emp_adj` among unique candidate positions (ties share the lowest rank; tie-group size in parentheses). Residues: C4 group | C3 group, from `position_scores.tsv`.

### Genotypic trait (`c4`)

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | n_hyp /100 | p.emp | p.emp_adj | rank /61 | sig |
|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:22 \| A:53 | 0.673 | 63 | 0.004 | 0.0576 | 1 (3) | **yes** |
| 665 | H→N | mutagenesis | N:23 \| H:52,N:2 | 0.683 | 100 | 0.004 | 0.0576 | 1 (3) | **yes** |
| 540 | P→T | selection | T:23 \| P:53,S:1 | 0.683 | 100 | 0.004 | 0.0576 | 1 (3) | **yes** |
| 572 | E→Q | selection | Q:18,K:5 \| E:52,Q:2 | 0.437 | 100 | 0.036 | 0.263 | 5 (14) | no |
| 733 | F→V | parallel | not in candidate set | | | | | | |
| 761 | S→A | parallel | A:14,S:8 \| S:53 | 0.226 | 63 | 0.169 | 0.329 | 25 (10) | no |
| 749 | L→T | weak | L:9,M:9,T:4 \| L:52,P:1 | 0.160 | 48 | 0.146 | 0.329 | 25 (10) | no |
| 505 | F→L | weak | L:17,F:6 \| F:54 | 0.343 | 100 | 0.170 | 0.329 | 25 (10) | no |
| 573 | A→N | weak | N:15,A:8 \| A:52,G:2 | 0.260 | 100 | 0.063 | 0.263 | 5 (14) | no |
| 731 | I→V | weak | V:17,Y:5 \| I:51,V:3 | 0.442 | 78 | 0.031 | 0.263 | 5 (14) | no |

**3/10 significant** (780, 665, 540), tied at rank 1; **6/10 present, not significant**; **1/10 absent** (733).

### Phenotypic trait (`c4_phenotypic`)

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | n_hyp /100 | p.emp | p.emp_adj | rank /51 | sig |
|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:16,A:2 \| A:49,S:1 | 0.715 | 68 | 0.0210 | 0.191 | 2 (6) | no |
| 665 | H→N | mutagenesis | N:17,H:3 \| H:47,N:3 | 0.713 | 89 | 0.0280 | 0.191 | 2 (6) | no |
| 540 | P→T | selection | T:17,P:3 \| P:48,S:1,T:1 | 0.710 | 92 | 0.0180 | 0.191 | 2 (6) | no |
| 572 | E→Q | selection | Q:17,E:3 \| E:47,Q:3 | 0.704 | 92 | 0.0280 | 0.191 | 2 (6) | no |
| 733 | F→V | parallel | not in candidate set | | | | | | |
| 761 | S→A | parallel | A:13,S:5 \| S:49,A:1 | 0.367 | 32 | 0.151 | 0.321 | 26 (1) | no |
| 749 | L→T | weak | M:9,L:6,T:3 \| L:48,P:1,T:1 | 0.114 | 24 | 0.266 | 0.379 | 36 (4) | no |
| 505 | F→L | weak | L:16,F:4 \| F:49,L:1 | 0.631 | 81 | 0.146 | 0.321 | 14 (12) | no |
| 573 | A→N | weak | N:14,A:6 \| A:49,N:1 | 0.355 | 39 | 0.202 | 0.349 | 27 (8) | no |
| 731 | I→V | weak | V:16,I:3 \| I:46,V:4 | 0.696 | 78 | 0.0360 | 0.216 | 8 (1) | no |

**0/10 significant; 9/10 present; 1/10 absent** (733). The four strongest truth positions (780, 665, 540, 572) share rank 2, behind one non-truth position.

### Position 733

Absent from the candidate set under both traits because the fixture carries no C4-specific residue there: C4 tips are `F:19, V:2, M:1, gap:1` (genotypic) and C3 tips `F:54`. The F→V change reported for grasses and sedges (Besnard et al. 2009, Table 2) is not a C4-group-wide state in this sequence sample. Its absence reflects the input, not an upstream filter.

## Significant non-truth positions

| trait | position | residues C4 \| C3 | CAAS_score | n_hyp /100 | sides | p.emp | p.emp_adj |
|---|---|---|---|---|---|---|---|
| genotypic | none | | | | | | |
| phenotypic | 859 | K:13,G:4 \| K:46,R:4 | 0.803 | 1 (H51) | bottom+top | 0.002 | 0.0659 |

Closest non-significant non-truth position under the genotypic trait: 818 (`E:21` \| `E:47,G:3,A:2,Q:1`, CAAS_score 0.077, 2 hypotheses, p.emp 0.009, p.emp_adj 0.108). Both are detected in ≤ 2 of 100 hypotheses. Hypothesis recurrence is descriptor-only in `scoring_compute.R` (§2g: it "never multiplies CAAS_score"), so neither `CAAS_score` nor `p.emp` penalises narrow detection; 859 reaches the highest `CAAS_score` in the phenotypic run from H51 alone. See `pepc_genotypic_vs_phenotypic.md` §5.

## Cross-method summary

| | FUBAR | PhyloPhere, genotypic | PhyloPhere, phenotypic |
|---|---|---|---|
| Truth positions with raw p < 0.05 (`p.emp`; FUBAR has no frequentist p) | n/a | 5/10 (780, 665, 540, 731, 572) | 5/10 (540, 780, 665, 572, 731) |
| Non-truth positions with raw p < 0.05 | n/a | 4 (818, 620, 588, 584) | 5 (859, 625, 751, 626, 460) |
| Truth positions significant (`p.emp_adj < 0.1`) | 0/10 | 3/10 | 0/10 |
| Present, not significant | 10/10 | 6/10 | 9/10 |
| Absent from method's output | 0/10 | 1/10 (733) | 1/10 (733) |
| Both mutagenesis sites (780, 665) significant | no | yes | no |
| Truth positions among 10 lowest `p.emp` / P[β>α]-ranked | 0 (best: 749, rank 33) | 5 | 5 |
| Significant non-truth sites | 2 | 0 | 1 |

## Caveat: the genotypic trait is the residue at 780

Under `c4`, the 23 C4 tips carry S at 780 (22) or a gap (1); the 54 C3 tips carry A (53) or a gap (1). Among observable residues the trait and the residue coincide exactly. Besnard et al. (2009) define "C4 ppc" by this residue (Fig. 3 caption: branches to "genes encoding C4 PEPC (with a serine at position 780) are in bold"), and Morel et al. (2024) state that the 78-sequence genotypic annotation was predicted "according to the presence or absence of the A780S mutation". Recovery of 780 under `c4` is therefore close to definitional. 540 and 665 co-vary almost perfectly with 780 across ppc-1 copies in this sample, and 572 largely (see the residue columns above; Besnard et al. 2009, Fig. 3), so their genotypic-run recovery inherits most of the same circularity. The phenotypic run is the non-circular test.

## Caveat: the truth set is itself genotype-derived

The selection-tier sites (540, 572) are Besnard et al. (2009) branch-site positive-selection codons on "C4 ppc" branches, which are defined by Ser780, and were confirmed by Morel et al. (2024) under the genotypic annotation. The weak tier is the same test without corroboration. Only the mutagenesis tier (780, 665; Bläsing et al. 2000; Svensson et al. 2003, via Morel et al. 2024) and the parallel tier (Christin et al. 2007) carry evidence independent of the A780S labelling.

## Caveat: contrast pairs are not sister pairs

`data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv`, both runs:

| | genotypic | phenotypic |
|---|---|---|
| hypotheses × pairs per hypothesis | 100 × 4 | 100 × 3 |
| distinct pairs | 135 | 106 |
| cross-genus pair-instances | 67/400 (17 %) | 147/300 (49 %) |
| same-species pair-instances (two ppc-1 copies of one species) | 41/400, in 35 hypotheses | 0/300 |

With `min_contrasts = 3`, a position must diverge in 3 of 4 pairs (genotypic) but in all 3 of 3 pairs (phenotypic). Under the genotypic trait, 41 pair-instances contrast the C4-type and non-C4 ppc-1 copy of the same *Eleocharis* or *Fimbristylis* species: 8 distinct pairs, all within *E. baldwinii*, *E. vivipara*, *F. dichotoma*, *F. ferruginea*, *F. littoralis*. These are paralog contrasts inside one genome, not lineage contrasts. *Cyperus* s.l. is paraphyletic with respect to *Kyllinga*, *Pycreus* and *Volkiella*, so a genus boundary between pair members does not imply separate lineages either; `min_contrasts` is not a count of independent C4 origins in either run.

## References

- Svensson P, et al. 2003. Cited via Morel et al. (2024) for the 665 functional claim; not available locally.

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bläsing OE, Westhoff P, Svensson P. 2000. Evolution of C4 phosphoenolpyruvate carboxylase in *Flaveria*: a conserved serine residue in the carboxyl-terminal part of the enzyme is a major determinant for C4-specific characteristics. J Biol Chem 275:27917–27923.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Christin PA, Salamin N, Savolainen V, Duvall MR, Besnard G. 2007. C4 photosynthesis evolved in grasses via parallel adaptive genetic changes. Curr Biol 17:1241–1247.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
