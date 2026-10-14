# Tier 1 PEPC: Results Report

Truth set: `validation/truthsets/tier1/pepc_c4.sites.tsv`, 10 positions in maize PEPC1 (CAA33317) numbering, 1:1 with the 970-column fixture alignment (`input/pepc/align/PEPC.fasta`, no reference row). Three methods are compared: OC's FUBAR (site-level dN/dS, trait-independent), PhyloPhere CAAS (CAAS discovery → CT_DISAMBIGUATION → SCORING, with the permulation null) and PhyloPhere FADE (directional selection on foreground-clade branches). The two PhyloPhere methods are run under two trait definitions on the same alignment and tree.

| PhyloPhere run | trait | definition | tips (in-group) | C4 / C3 |
|---|---|---|---|---|
| `results/c4_complete/` | `c4` | genotypic: tip's ppc-1 carries A780S (Besnard et al. 2009; Morel et al. 2024) | 77 | 23 / 54 |
| `results/c4_phenotypic_complete/` | `c4_phenotypic` | phenotypic: species-level C3/C4 call from Bruhl & Wilson (2007); 7 C3/C4-intermediate *Eleocharis* accessions pruned | 70 | 20 / 50 |

Both traits are columns of `input/pepc/my_traits.tsv`; sourcing of the phenotypic column is documented in `input/.pepc_phenotypic/README.md`. The outgroup *Chrysitrix dodii* is in the alignment and tree but carries no trait value in either run. Species are mapped to NCBI tax ids with the fixture's `input/pepc/taxid.tsv`. Cross-trait interpretation is in `pepc_genotypic_vs_phenotypic.md`; runtime and run configuration in `pepc_resources.md`.

## Numbering

`scoring/position_scores.tsv` `Position` and FADE's `fade_site_bf_top.tsv` `position` are 0-based alignment columns, i.e. maize position minus 1, in both runs. Verified per truth position by recomputing C4/C3 residue counts at alignment column *p* (1-based) from the raw FASTA under each run's tip set and trait and matching them to `top_species_residues` / `bottom_species_residues` at `Position = p − 1`. All positions below are in maize numbering.

## `position_scores.tsv` schema and statistics

Columns: `Gene, Position, n_schemes, scheme_set, GS_schemes, is_biochem, n_hypotheses, participating_hypotheses, top_species_residues, bottom_species_residues, n_top_species, n_bottom_species, n_conserved_pairs, CAAS_score, side, caas, amino_encoded, ancestral_aa, derived_aa, p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact, nd, null_unseen`. One row per (Position, side); a position detected on both sides has two rows sharing one `p.emp`, `p.adj_bh`, `p.emp_fact` and `p.adj_bh_fact`. `top_species_residues`, `bottom_species_residues`, `n_top_species`, `n_bottom_species` and `n_conserved_pairs` tally the species of the contrast pairs; the residue columns of the tables below count every tip with a trait value, read from the alignment (gaps omitted). Definitions are in `docs/scoring_v2_p_emp.md`.

- `p.emp`: `(k_emp + 1)/(N + 1)`, where `k_emp` counts null cycles that re-detect the position on any side with max-over-sides CAAS ≥ the observed max-over-sides `CAAS_score`; `N` = 1000. A position with a score of 0 has `p.emp` = 1.
- `p.adj_bh`: BH with one test per position over the null universe, i.e. every position detected in ≥ 1 null cycle plus every observed position, with p = 1 for null-detected positions the observed data did not detect.
- `p.emp_fact`: `(nd + 1)/(N + 1) × (1 + #{detections of the class with score ≥ s})/(1 + #{detections of the class})`, where `nd` is the number of null cycles that score the position and the class is one of 20 percentile classes of the null detections. It is not bounded below by `1/(N + 1)`.
- `p.adj_bh_fact`: BH of `p.emp_fact` over the same null universe as `p.adj_bh`.
- Null design: every null cycle is matched to the observed design, 100 FOP hypotheses of K pairs each (K = 4 genotypic, 3 phenotypic); `caas_permulation/resample_perms.tab` holds 100 hypothesis labellings for each of the 1000 cycles.
- The four p-values (`p.emp`, `p.adj_bh`, `p.emp_fact`, `p.adj_bh_fact`) are informative: the pipeline derives no threshold or significance column from them. This report marks values below 0.05 in bold and counts positions below 0.05 and 0.1 only to describe them; `p.emp` at its floor `1/(N + 1)` is written `<0.001`. `CAAS_score` is the exact-identity (US) score of a side. Each domain of the design (one per canonical pair) is scored by the species that represent it across the hypotheses: a species takes the mean, over the hypotheses in which it represents the domain, of the noisy-OR of the independence of that domain from its partner domains (the probability that the derived residue was not already present at their last common ancestor), with 0 where the hypothesis does not detect the position; the domain takes the mean over its species weighted by the PSS of their pairs; the score is the mean over the K domains. Every species counts once whatever the number of hypotheses it appears in, and the score lies between 0 and 1. The grouping schemes (GS1 to GS4) do not enter the score; `scheme_set`, `GS_schemes` and `is_biochem` describe the biochemical profile of each position. `nd` is the number of null cycles that score the position.

Note for ad hoc pandas reads: the `caas` value `N/A` (maize 573) is parsed as missing under pandas defaults; read with `keep_default_na=False`.

| | genotypic | phenotypic |
|---|---|---|
| position-rows | 36 | 34 |
| unique positions (candidate set) | 35 | 34 |
| BH family (null universe) | 80 (45 null-only) | 71 (37 null-only) |
| null cycles `N` | 1000 | 1000 |
| minimum attainable `p.emp` | 0.000999 | 0.000999 |
| positions with `p.emp < 0.05` / `p.emp_fact < 0.05` | 28 / 27 | 26 / 24 |
| positions with `p.adj_bh < 0.05` / `< 0.1` | 26 / 27 | 19 / 23 |
| positions with `p.adj_bh_fact < 0.05` / `< 0.1` | 20 / 27 | 16 / 19 |
| minimum `p.adj_bh` / `p.adj_bh_fact` | 0.0054 / 0.0044 | 0.0053 / 0.0049 |
| observed positions vs null per cycle | 35; null median 25, q95 35; P(null ≥ obs) = 0.065 | 34; null median 19, q95 35; P = 0.064 |

## Method 1: OC / FUBAR

HyPhy FUBAR on 78 sequences × 970 codons (`input/pepc/oc_run/PSEL/PEPC.FUBAR.json`). Trait-independent, so one run applies to both trait definitions. Threshold: posterior P[β > α] > 0.9.

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

**0/10.** FUBAR flags 2 sites (codon 630, P = 0.941; codon 474, P = 0.924), neither in the truth set.

## Method 2: PhyloPhere CAAS

Score rank = competition rank of `CAAS_score` among unique candidate positions (best side). BH rank = competition rank of `p.adj_bh` (tie-group size in parentheses). Residues: C4 group | C3 group, over every tip with a trait value. Values below 0.05 in bold.

### Genotypic trait (`c4`)

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | score rank /35 | n_hyp /100 | p.emp | p.adj_bh | BH rank | p.emp_fact | p.adj_bh_fact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:22 \| A:53 | 0.594 | 6 | 61 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0044** |
| 665 | H→N | mutagenesis | N:23 \| H:52,N:2 | 1.000 | 1 | 100 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0044** |
| 540 | P→T | selection | T:23 \| P:53,S:1 | 1.000 | 1 | 100 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0065** |
| 572 | E→Q | selection | Q:18,K:5 \| E:52,Q:2 | not in candidate set | | | | | | | |
| 733 | F→V | parallel | not in candidate set | | | | | | | | |
| 761 | S→A | parallel | A:14,S:8 \| S:53 | 0.288 | 12 | 61 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0044** |
| 749 | L→T | weak | L:9,M:9,T:4 \| L:52,P:1 | 0.025 | 32 | 37 | **0.016** | **0.048** | 25 (1) | 0.083 | 0.20 |
| 505 | F→L | weak | L:17,F:6 \| F:54 | 0.667 | 5 | 100 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0044** |
| 573 | A→N | weak | N:15,A:8 \| A:52,G:2 | not in candidate set | | | | | | | |
| 731 | I→V | weak | V:17,Y:5 \| I:51,V:3 | 0.560 | 7 | 77 | **<0.001** | **0.0054** | 1 (14) | **<0.001** | **0.0044** |

**Adjusted p below 0.05:** `p.adj_bh` 7/10 (505, 540, 665, 731, 749, 761, 780), `p.adj_bh_fact` 6/10 (505, 540, 665, 731, 761, 780; 749 at 0.20); **7/10 present, 3/10 absent** (572, 573, 733). Six of the seven present positions (780, 665, 540, 761, 505, 731) are exceeded by none of the 1000 null cycles and share the floor value 0.0054 with eight other positions, so their order rests on the score. 540 and 665 hold the maximum (1.000, rank 1); the non-truth 611 and 620 rank 3 and 4 (0.826, 0.690), and 505, 780 and 731 rank 5 to 7 (0.667 to 0.560). In the candidate set 26 of 35 positions are below 0.05 under `p.adj_bh`, 19 of them outside the truth set.

### Phenotypic trait (`c4_phenotypic`)

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | score rank /34 | n_hyp /100 | p.emp | p.adj_bh | BH rank | p.emp_fact | p.adj_bh_fact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:16,A:2 \| A:49,S:1 | 0.609 | 8 | 73 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0049** |
| 665 | H→N | mutagenesis | N:17,H:3 \| H:47,N:3 | 0.781 | 4 | 92 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0052** |
| 540 | P→T | selection | T:17,P:3 \| P:48,T:1,S:1 | 0.807 | 3 | 97 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0049** |
| 572 | E→Q | selection | Q:17,E:3 \| E:47,Q:3 | 0.814 | 2 | 97 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0049** |
| 733 | F→V | parallel | not in candidate set | | | | | | | | |
| 761 | S→A | parallel | A:13,S:5 \| S:49,A:1 | 0.285 | 14 | 47 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0052** |
| 749 | L→T | weak | M:9,L:6,T:3 \| L:48,P:1,T:1 | 0.007 | 33 | 29 | 0.109 | 0.24 | 31 (1) | 0.132 | 0.29 |
| 505 | F→L | weak | L:16,F:4 \| F:49,L:1 | 0.665 | 7 | 84 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0052** |
| 573 | A→N | weak | N:14,A:6 \| A:49,N:1 | 0.363 | 11 | 58 | **<0.001** | **0.0053** | 1 (13) | **0.0035** | **0.016** |
| 731 | I→V | weak | V:16,I:3 \| I:46,V:4 | 0.693 | 5 | 85 | **<0.001** | **0.0053** | 1 (13) | **<0.001** | **0.0049** |

**Adjusted p below 0.05:** `p.adj_bh` 8/10 (780, 665, 540, 572, 761, 505, 573, 731), `p.adj_bh_fact` 8/10 (the same eight; 573 at 0.016); **9/10 present, 1/10 absent** (733). The eight positions are exceeded by none of the 1000 null cycles and share the floor value 0.0053 with five other positions; 749 (score 0.007, 29 hypotheses) is the exception (`p.emp` 0.109). The non-truth 518 ranks first (0.873); 572, 540, 665, 731 follow at ranks 2 to 5 (0.814 to 0.693), then the non-truth 611 (0.689), 505 (0.665) and 780 (rank 8, 0.609). In the candidate set 19 of 34 positions are below 0.05 under `p.adj_bh`, 11 of them outside the truth set.

### Position 733

Absent from the candidate set under both traits because the fixture carries no C4-specific residue there: C4 tips are `F:19, V:2, M:1, gap:1` (genotypic) and C3 tips `F:54`. The F→V change reported for grasses and sedges (Besnard et al. 2009, Table 2) is not a C4-group-wide state in this sequence sample.

### Positions 572 and 573 under the genotypic trait

Absent from the candidate set of the genotypic run: the C4 tips carry two residues at each (572: Q 18, K 5; 573: N 15, A 8), and only the exact-identity scheme detects positions. They are detected under the phenotypic trait (572 and 573 in 97 and 58 hypotheses).

### Null of each truth position

`k_emp` = null cycles that re-detect the position and reach ≥ the observed max-over-sides `CAAS_score`; quantiles over the cycles that detect it (`pepc_pvalue_tables.py`).

| pos | trait | observed CAAS | n_hyp /100 | null detects | null q50 | null q90 | null q99 | k_emp | p.emp |
|---|---|---|---|---|---|---|---|---|---|
| 780 | geno | 0.594 | 61 | 868 | 0.035 | 0.119 | 0.226 | 0 | <0.001 |
| 780 | pheno | 0.609 | 73 | 543 | 0.054 | 0.152 | 0.298 | 0 | <0.001 |
| 665 | geno | 1.000 | 100 | 956 | 0.052 | 0.169 | 0.313 | 0 | <0.001 |
| 665 | pheno | 0.781 | 92 | 670 | 0.067 | 0.197 | 0.442 | 0 | <0.001 |
| 540 | geno | 1.000 | 100 | 824 | 0.048 | 0.145 | 0.295 | 0 | <0.001 |
| 540 | pheno | 0.807 | 97 | 572 | 0.067 | 0.196 | 0.425 | 0 | <0.001 |
| 572 | pheno | 0.814 | 97 | 576 | 0.071 | 0.203 | 0.438 | 0 | <0.001 |
| 761 | geno | 0.288 | 61 | 521 | 0.023 | 0.068 | 0.153 | 0 | <0.001 |
| 761 | pheno | 0.285 | 47 | 379 | 0.025 | 0.074 | 0.146 | 0 | <0.001 |
| 749 | geno | 0.025 | 37 | 211 | 0.008 | 0.018 | 0.050 | 15 | 0.016 |
| 749 | pheno | 0.007 | 29 | 143 | 0.010 | 0.023 | 0.036 | 108 | 0.109 |
| 505 | geno | 0.667 | 100 | 719 | 0.055 | 0.157 | 0.268 | 0 | <0.001 |
| 505 | pheno | 0.665 | 84 | 567 | 0.070 | 0.193 | 0.362 | 0 | <0.001 |
| 573 | pheno | 0.363 | 58 | 423 | 0.040 | 0.123 | 0.221 | 0 | <0.001 |
| 731 | geno | 0.560 | 77 | 685 | 0.030 | 0.091 | 0.201 | 0 | <0.001 |
| 731 | pheno | 0.693 | 85 | 563 | 0.060 | 0.162 | 0.334 | 0 | <0.001 |

Every truth position that is a candidate, except 749, has an observed score above the null q99 of its own position (by a factor of 1.6 to 3.4) and is exceeded by none of the 1000 cycles. The pooled null quartiles over all detected positions are 0.011 / 0.026 / 0.056 (genotypic) and 0.013 / 0.029 / 0.065 (phenotypic) against observed 0.056 / 0.117 / 0.421 and 0.031 / 0.122 / 0.447. `CAAS_score` is an average over the species of the design and its scale depends on the design; the p-values, calibrated within each run, are what compare across runs.

### Non-truth positions: highest scores and positions discussed in the text

Positions with an adjusted p below 0.05 outside the truth set number 19 (`p.adj_bh`) and 14 (`p.adj_bh_fact`) in the genotypic run, and 11 and 8 in the phenotypic run; the table lists the highest-scoring ones and those discussed below. Residues are the C4 and C3 tallies over every tip with a trait value, read from the alignment.

| trait | position | residues C4 \| C3 | CAAS_score | score rank | n_hyp /100 | p.emp | p.adj_bh | p.adj_bh_fact | FADE BF |
|---|---|---|---|---|---|---|---|---|---|
| genotypic | 611 | L:22,F:1 \| F:44,L:10 | 0.826 | 3 | 90 | **<0.001** | **0.0054** | **0.0044** | 7 644 |
| genotypic | 620 | C:17,A:3,T:1,F:1,S:1 \| S:51,T:3 | 0.690 | 4 | 100 | **<0.001** | **0.0054** | **0.0044** |  |
| genotypic | 511 | I:15,V:8 \| I:54 | 0.500 | 8 | 100 | **<0.001** | **0.0054** | **0.0049** |  |
| genotypic | 751 | F:11,Y:11 \| Y:53 | 0.445 | 9 | 61 | **<0.001** | **0.0054** | **0.0044** |  |
| genotypic | 517 | T:14,A:9 \| T:53,A:1 | 0.397 | 10 | 69 | **0.0030** | **0.012** | **0.0086** | 104 |
| genotypic | 509 | D:12,E:11 \| E:51,D:2,G:1 | 0.387 | 11 | 64 | **<0.001** | **0.0054** | **0.0044** | 138 |
| genotypic | 584 | T:8,M:7,S:5,I:2,K:1 \| I:52,V:2 | 0.257 | 13 | 100 | **<0.001** | **0.0054** | **0.019** |  |
| genotypic | 839 | G:13,K:8 \| G:53 | 0.201 | 14 | 41 | **0.0020** | **0.0088** | **0.0082** | 924 |
| genotypic | 518 | F:18,L:5 \| L:42,F:12 | 0.117 | 18 | 24 | **0.0090** | **0.032** | **0.049** |  |
| phenotypic | 518 | F:20 \| L:42,F:8 | 0.873 | 1 | 100 | **<0.001** | **0.0053** | **0.0049** | 68 828 |
| phenotypic | 611 | L:17,F:3 \| F:41,L:9 | 0.689 | 6 | 84 | **<0.001** | **0.0053** | **0.0052** | 1 196 |
| phenotypic | 620 | C:13,S:4,A:2,T:1 \| S:46,T:3,A:1 | 0.471 | 9 | 86 | **<0.001** | **0.0053** | **0.0049** | 1 208 |
| phenotypic | 579 | T:13,A:5,E:2 \| A:49,T:1 | 0.373 | 10 | 75 | **0.0020** | **0.0098** | **0.0052** |  |
| phenotypic | 751 | Y:12,F:6 \| Y:50 | 0.302 | 13 | 50 | **<0.001** | **0.0053** | **0.0052** |  |
| phenotypic | 501 | I:20 \| F:39,I:11 | 0.282 | 15 | 42 | **<0.001** | **0.0053** | **0.0049** |  |
| phenotypic | 588 | I:15,L:4,Y:1 \| R:34,K:13,I:2,L:1 | 0.259 | 16 | 94 | **0.0040** | **0.016** | **0.0052** |  |
| phenotypic | 744 | R:11,L:7,S:1 \| L:48,R:1,S:1 | 0.104 | 18 | 39 | **0.0030** | **0.013** | **0.033** |  |
| phenotypic | 474 | G:10,E:2 \| E:16,D:4,G:4,K:1 | 0.005 | 34 | 1 | **0.011** | **0.040** | 0.13 | 428 |

The two adjustments rank positions differently. `p.adj_bh` follows each position's own null (`p.emp`); `p.adj_bh_fact` follows the observed score against the scores of the null detections of positions with a similar detection frequency. A position that the null rarely detects reaches a small `p.emp` with a low score: 474 (detected in 1 of 100 hypotheses, score 0.005) has `p.emp` 0.011 and `p.adj_bh` 0.040 but `p.adj_bh_fact` 0.129. A call that passes only one adjustment is the least secure.

- **518** is present in all 100 phenotypic hypotheses and ranks first (0.873); every phenotypic-C4 tip carries F. Under the genotypic trait the C4 group is `F:18, L:5` and the position is weak (`CAAS_score` 0.117, 24 hypotheses, `p.adj_bh` 0.032, `p.adj_bh_fact` 0.049). The five L tips are the Ser780-carrying *E. baldwinii* (FM208014/015/016) and *E. vivipara* (AB085948, FM208029) accessions, which the phenotypic run prunes as C3/C4 intermediates. F therefore marks the C4 lineages other than *Eleocharis*; it is not associated with the C4-type ppc-1 copy in *Eleocharis*. Within C4 *Fimbristylis* both ppc-1 copies carry F (the non-C4 paralogs FM208032/034/036 included), so at 518 F is a lineage state rather than a property of the C4-recruited copy. CT_DISAMBIGUATION's ancestral reconstruction (`ct_disambiguation/caas_convergence_master.csv`, US scheme, top side) places L at the MRCA of the three C4 domains carrying F in 100, 93 and 74 of the 100 hypotheses (posterior 0.94 to 1.0; in the other hypotheses the third domain's MRCA is F, with posterior 0.56 to 0.57), and FADE maps three L→F substitutions on foreground branches, so the association reflects repeated gains of F rather than its retention. F also occurs in the outgroup and in early-diverging C3 lineages (*Coleochloa*, *Microdracoides*, *Carpha*), consistent with an early shift to L in the in-group followed by returns to F in the C4 lineages.
- **460** and **474** are detected in 2 and 1 of 100 hypotheses and have low coverage (460: 12 of 20 C4 and 25 of 50 C3 tips carry a residue). Their scores are 0.016 and 0.005, because a hypothesis that does not detect the position contributes zero; what lets 474 reach a `p.emp` of 0.011 is its low null detection rate.

## Method 3: PhyloPhere FADE

`selection/fade/top/fade_site_bf_top.tsv` and `json/PEPC.top.FADE.json`, top (C4) direction, Bayes factor threshold 100 (`fade_bf_threshold`). FADE tests, per site and target residue, for a substitution bias towards that residue on foreground branches relative to the background branches.

Foreground (`fade_internal_nodes = all_descendants`): every C4 tip's terminal branch plus every internal branch whose descendant tips are all C4, stem branches of C4 clades included (HyPhy LabelTrees' default "All descendants" strategy, and the foreground used for FADE by Morel et al. 2024 on the same data). Background (`fade_background_scope = all`): every other branch of the tree, the outgroup included in the genotypic run; the phenotypic tree carries no outgroup because pruning removes species without a trait value.

| | genotypic | phenotypic |
|---|---|---|
| foreground branches (terminal + internal) / total | 41 (23 + 18) / 154 | 35 (20 + 15) / 138 |

### Truth set

BF towards the truth set's derived residue; substitutions are FADE's reconstructed history on foreground branches (`Substitutions` site annotation).

| position | tier | target | genotypic BF | subs (geno) | phenotypic BF | subs (pheno) |
|---|---|---|---|---|---|---|
| 780 | mutagenesis | S | **1.3 × 10⁶** | A→S ×5 | **760** | A→S ×3 |
| 665 | mutagenesis | N | **5.9 × 10⁷** | H→N ×5 | **5 100** | H→N ×3 |
| 540 | selection | T | **1.3 × 10⁹** | P→T ×5 | **48 900** | P→T ×3 |
| 572 | selection | Q | **8 510** | E→Q ×4, E→K ×1 | **1 760** | E→Q ×3 |
| 733 | parallel | V | 2.3 | F→V, F→M | 2.3 | F→V, F→M |
| 761 | parallel | A | 64 | S→A ×3 | 12 | S→A ×2 |
| 749 | weak | T | **168** | L/M→T ×3 | 39 | L/M→T ×2 |
| 505 | weak | L | **547** | F→L ×4 | **174** | F→L ×3 |
| 573 | weak | N | **17 400** | A→N ×3 | **360** | A→N ×2 |
| 731 | weak | V | **281** | I→V ×4 | **122** | I→V ×3 |

**8/10 genotypic, 7/10 phenotypic.** Both mutagenesis sites and both selection-tier sites are significant under both traits, each with 3 to 5 reconstructed substitutions to the derived residue on foreground branches, i.e. repeated changes across C4 lineages. Under the phenotypic trait 780 keeps BF 760 despite the two A-carrying phenotypic-C4 paralogs and the S-carrying C3 *C. distichus*.

### Non-truth sites at BF ≥ 100

| trait | position | target | BF | residues C4 \| C3 | subs | CAAS `p.adj_bh` / `p.adj_bh_fact` |
|---|---|---|---|---|---|---|
| genotypic | 611 | L | 7 640 | L:22,F:1 \| F:44,L:10 | F→L ×5 | **0.0054** / **0.0044** |
| genotypic | 579 | T | 6 180 | T:19,A:2,E:2 \| A:53,T:1 | A→T ×4 | not detected |
| genotypic | 839 | K | 924 | G:13,K:8 \| G:53 | G→K ×3 | **0.0088** / **0.0082** |
| genotypic | 509 | D | 138 | D:12,E:11 \| E:51,D:2,G:1 | E→D ×4 | **0.0054** / **0.0044** |
| genotypic | 514 | C | 128 | V:21,C:2 \| V:53,I:1 | V→C ×2 | not detected |
| genotypic | 630 | K | 114 | Q:12,K:11 \| Q:43,K:9 | Q→K ×3 | not detected |
| genotypic | 471 | T | 112 | T:9,K:5 \| T:26,Q:2 (35/77 gapped) | K→T ×2, Q→T, T→K | not detected |
| genotypic | 517 | A | 104 | T:14,A:9 \| T:53,A:1 | T→A ×3 | **0.012** / **0.0086** |
| phenotypic | 518 | F | 68 800 | F:20 \| L:42,F:8 | L→F ×3 | **0.0053** / **0.0049** |
| phenotypic | 620 | C | 1 210 | C:13,S:4,A:2,T:1 \| S:46,T:3,A:1 | S→C, A→C, others | **0.0053** / **0.0049** |
| phenotypic | 611 | L | 1 200 | L:17,F:3 \| F:41,L:9 | F→L ×5 | **0.0053** / **0.0052** |
| phenotypic | 474 | G | 428 | G:10,E:2 \| E:16,D:4,G:4,K:1 (33/70 gapped) | E→G ×2, G→E | **0.040** / 0.13 |
| phenotypic | 630 | K | 146 | K:11,Q:9 \| Q:39,K:9 | Q→K ×3 | not detected |
| phenotypic | 514 | C | 129 | V:18,C:2 \| V:49,I:1 | V→C ×2 | not detected |

- **579, 611, 620** combine several reconstructed substitutions with a residue shared by C4 lineages that are not each other's closest relatives (*Cyperus*, *Eleocharis*, *Rhynchospora*, *Fimbristylis*, *Bulbostylis*). 611's L also occurs in about ten C3 tips across unrelated clades, so it is the least C4-specific of the three.
- **630** is also one of the two sites FUBAR flags for positive selection (P = 0.941). Residues are close to evenly split in C4, and K is present in nine C3 tips.
- **471** and **474** are heavily gapped and rest on 2 to 4 reconstructed changes; **514** rests on two V→C changes in two *Cyperus* tips. These are the weakest non-truth calls.
- **880** is not significant under either trait (BF 17 and 18): its I is a clade state of C4 *Cyperus*, reconstructed as one change on the clade stem. It is a candidate only in the genotypic run (`CAAS_score` 0.217, 20 hypotheses, `p.adj_bh` 0.025), and the phenotypic run does not report it.

### CAAS vs FADE overlap

| | genotypic | phenotypic |
|---|---|---|
| candidate positions | 35 | 34 |
| `p.adj_bh < 0.05` | 26 | 19 |
| `p.adj_bh_fact < 0.05` | 20 | 16 |
| FADE BF ≥ 100 | 16 (10 of them candidates) | 13 (11 of them candidates) |
| `p.adj_bh < 0.05` and FADE | 10 (505, 509, 517, 540, 611, 665, 731, 749, 780, 839) | 11 (474, 505, 518, 540, 572, 573, 611, 620, 665, 731, 780) |
| `p.adj_bh < 0.05` without FADE | 16 | 8 (501, 503, 509, 579, 588, 744, 751, 761) |
| `p.adj_bh_fact < 0.05` and FADE | 9 (505, 509, 517, 540, 611, 665, 731, 780, 839) | 10 (505, 518, 540, 572, 573, 611, 620, 665, 731, 780) |
| `p.adj_bh_fact < 0.05` without FADE | 11 | 6 (501, 579, 588, 744, 751, 761) |

Under the genotypic trait FADE confirms 10 of the 26 positions below 0.05 under `p.adj_bh` and 9 of the 20 under `p.adj_bh_fact`; under the phenotypic trait, 11 of 19 and 10 of 16. FADE calls fewer candidates than CAAS in both runs: 10 of 35 and 11 of 34 candidates, and 6 and 2 FADE calls fall outside the candidate set (the two of the phenotypic run are 514 and 630). Among the truth positions FADE calls 8/10 (genotypic) and 7/10 (phenotypic); CAAS has 7/10 present with adjusted p below 0.05 in the genotypic run and 8/10 in the phenotypic run. The two methods share the alignment and the foreground definition, so their agreement is not independent corroboration; it reflects agreement between two models of the same labelled data (pair contrasts under a design-matched permulation null vs a branch-level substitution-bias model).

## Cross-method summary

| | FUBAR | CAAS, genotypic | CAAS, phenotypic | FADE, genotypic | FADE, phenotypic |
|---|---|---|---|---|---|
| Truth positions with adjusted p below 0.05 (FADE: BF ≥ 100) | 0/10 | 7/10 (`p.adj_bh`); 6/10 (`p.adj_bh_fact`) | 8/10; 8/10 | 8/10 | 7/10 |
| Both mutagenesis sites (780, 665) | no | yes | yes | yes | yes |
| Non-truth sites (FUBAR: P > 0.9; FADE: BF ≥ 100) | 2 | 19 (`p.adj_bh`); 14 (`p.adj_bh_fact`) | 11; 8 | 8 | 6 |
| Truth positions absent from output | 0/10 | 3/10 (572, 573, 733) | 1/10 (733) | | |

## Caveat: the genotypic trait is the residue at 780

Under `c4`, the 23 C4 tips carry S at 780 (22) or a gap (1); the 54 C3 tips carry A (53) or a gap (1). Among observable residues the trait and the residue coincide exactly. Besnard et al. (2009) define "C4 ppc" by this residue, and Morel et al. (2024) predicted the 78-sequence genotypic annotation from the presence or absence of A780S. Recovery of 780 under `c4` is therefore close to definitional for both CAAS and FADE. 540 and 665 co-vary almost perfectly with 780 across ppc-1 copies in this sample, and 572 largely, so their genotypic-run recovery inherits most of the same circularity. The phenotypic run is the non-circular test.

## Caveat: the truth set is itself genotype-derived

The selection-tier sites (540, 572) are Besnard et al. (2009) branch-site positive-selection codons on "C4 ppc" branches, which are defined by Ser780, and were confirmed by Morel et al. (2024) under the genotypic annotation. The weak tier is the same test without corroboration. Only the mutagenesis tier (780, 665; Bläsing et al. 2000; Svensson et al. 2003, via Morel et al. 2024) and the parallel tier (Christin et al. 2007) carry evidence independent of the A780S labelling. Positions outside the truth set are not thereby false positives: the set is not exhaustive.

## Caveat: contrast pairs are not sister pairs

`data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv`, both runs:

| | genotypic | phenotypic |
|---|---|---|
| hypotheses × pairs per hypothesis | 100 × 4 | 100 × 3 |
| distinct pairs | 126 | 130 |
| cross-genus pair-instances | 50/400 (13 %) | 138/300 (46 %) |
| same-species pair-instances (two ppc-1 copies of one species) | 23/400, in 21 hypotheses | 0/300 |

With `min_contrasts = 3`, a position must diverge in 3 of 4 pairs (genotypic) but in all 3 of 3 pairs (phenotypic). Under the genotypic trait, same-species pair-instances contrast the C4-type and non-C4 ppc-1 copy of one *Eleocharis* or *Fimbristylis* species: paralog contrasts inside one genome, not lineage contrasts. *Cyperus* s.l. is paraphyletic with respect to *Kyllinga*, *Pycreus* and *Volkiella*, so a genus boundary between pair members does not imply separate lineages either; `min_contrasts` is not a count of independent C4 origins in either run.

## References

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bläsing OE, Westhoff P, Svensson P. 2000. Evolution of C4 phosphoenolpyruvate carboxylase in *Flaveria*: a conserved serine residue in the carboxyl-terminal part of the enzyme is a major determinant for C4-specific characteristics. J Biol Chem 275:27917–27923.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Christin PA, Salamin N, Savolainen V, Duvall MR, Besnard G. 2007. C4 photosynthesis evolved in grasses via parallel adaptive genetic changes. Curr Biol 17:1241–1247.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
- Svensson P, et al. 2003. Cited via Morel et al. (2024) for the 665 functional claim; not available locally.
