# Tier 1 PEPC: Genotypic vs Phenotypic Trait

Two PhyloPhere runs on the same 970-column ppc-1 alignment and tree, differing in the trait. Per-run tables, the FUBAR baseline, FADE (including its agreement with CAAS) and numbering verification are in `pepc_results.md`; runtime and run configuration in `pepc_resources.md`.

| | genotypic | phenotypic |
|---|---|---|
| results dir | `output/pepc/results/c4_complete/` | `output/pepc/results/c4_phenotypic_complete/` |
| trait column (`input/pepc/my_traits.tsv`) | `c4` | `c4_phenotypic` |
| definition | ppc-1 copy carries Ser780 (Besnard et al. 2009) | species-level anatomical/physiological call (Bruhl & Wilson 2007), assigned to every ppc-1 copy of that species (as in Morel et al. 2024) |
| tips (in-group) | 77 (23 C4 / 54 C3) | 70 (20 C4 / 50 C3); 7 *E. baldwinii* / *E. vivipara* C3/C4-intermediate accessions pruned |
| evolutionary model (AIC) | BM | BM |
| candidate positions (unique) | 35 lenient, 10 strict | 34 lenient, 17 strict |
| hypotheses × pairs | 100 × 4 | 100 × 3 |
| null cycles, hypotheses per cycle | 1000, 100 each | 1000, 100 each |

Each trait is run in two CAAS discovery modes (`pepc_results.md`): **lenient** (`min_divergent_fraction` 0.5, up to 2 of 4 or 1 of 3 pairs may keep the ancestral state) and **strict** (every pair diverges). Both modes use the same hypotheses and null labelings. Tables are lenient unless a section says strict.

Position numbers are maize PEPC1 (`position_scores.tsv` `Position` + 1). The four position p-values are informative and carry no threshold in the pipeline: `p.emp`, its BH adjustment over the null universe `p.adj_bh` (one test per position), the factorized `p.emp_fact` and its BH adjustment `p.adj_bh_fact` over the same family (§6c). This report marks values below 0.05 and counts positions below 0.05 and 0.1 to describe them. `CAAS_score` is the exact-identity (US) score of a side: each domain is scored by the species that represent it across the hypotheses (a species takes the mean of its domain score over the hypotheses in which it represents the domain, 0 where the hypothesis does not detect the position; species weighted by the PSS of their pairs), and the score is the mean over the K domains; it lies between 0 and 1 and does not include the grouping schemes. `pepc_pvalue_tables.py` prints the tables of this report from the results directory.

## 1. Truth positions against their nulls

`k_emp` = null cycles that re-detect the position and reach ≥ the observed max-over-sides `CAAS_score`. Null quantiles are over the cycles that detect the position (`pepc_results.md`, "Null of each truth position", has the full table). `k_emp` reproduces the reported `p.emp` in every case.

| pos | trait | observed CAAS | n_hyp /100 | null detects | null q50 | null q90 | null q99 | k_emp | p.emp | p.adj_bh |
|---|---|---|---|---|---|---|---|---|---|---|
| 780 | geno | 0.594 | 61 | 868 | 0.035 | 0.119 | 0.226 | 0 | <0.001 | **0.0054** |
| 780 | pheno | 0.609 | 73 | 543 | 0.054 | 0.152 | 0.298 | 0 | <0.001 | **0.0053** |
| 665 | geno | 1.000 | 100 | 956 | 0.052 | 0.169 | 0.313 | 0 | <0.001 | **0.0054** |
| 665 | pheno | 0.781 | 92 | 670 | 0.067 | 0.197 | 0.442 | 0 | <0.001 | **0.0053** |
| 540 | geno | 1.000 | 100 | 824 | 0.048 | 0.145 | 0.295 | 0 | <0.001 | **0.0054** |
| 540 | pheno | 0.807 | 97 | 572 | 0.067 | 0.196 | 0.425 | 0 | <0.001 | **0.0053** |
| 572 | geno | not in candidate set | | | | | | | | |
| 572 | pheno | 0.814 | 97 | 576 | 0.071 | 0.203 | 0.438 | 0 | <0.001 | **0.0053** |
| 731 | geno | 0.560 | 77 | 685 | 0.030 | 0.091 | 0.201 | 0 | <0.001 | **0.0054** |
| 731 | pheno | 0.693 | 85 | 563 | 0.060 | 0.162 | 0.334 | 0 | <0.001 | **0.0053** |
| 505 | geno | 0.667 | 100 | 719 | 0.055 | 0.157 | 0.268 | 0 | <0.001 | **0.0054** |
| 505 | pheno | 0.665 | 84 | 567 | 0.070 | 0.193 | 0.362 | 0 | <0.001 | **0.0053** |

- **Both traits.** The positions in the table are exceeded by none of the 1000 null cycles (`p.emp` at its floor), and their observed score is 2.5 to 3.4 times the null q99 of the position (genotypic) and 1.8 to 2.1 times (phenotypic).
- **Null scale.** Null medians at the truth positions are 0.008 to 0.055 (genotypic) and 0.010 to 0.071 (phenotypic); the pooled null quartiles over all detected positions are 0.011 / 0.026 / 0.056 and 0.013 / 0.029 / 0.065, and the observed quartiles 0.056 / 0.117 / 0.421 and 0.031 / 0.122 / 0.447. `CAAS_score` is an average over the species of the design and its scale depends on the design; the p-values, calibrated within each run, are what compare across runs.
- The phenotypic null detects each truth position less often (543 to 670 cycles vs 685 to 956 for the five positions detected in both runs), consistent with 3-of-3 all-or-nothing detection per hypothesis (§4).

## 2. Truth set, side by side

| pos | tier | geno score rank /35 | geno p.adj_bh | geno p.adj_bh_fact | pheno score rank /34 | pheno p.adj_bh | pheno p.adj_bh_fact |
|---|---|---|---|---|---|---|---|
| 780 | mutagenesis | 6 | **0.0054** | **0.0044** | 8 | **0.0053** | **0.0049** |
| 665 | mutagenesis | 1 | **0.0054** | **0.0044** | 4 | **0.0053** | **0.0052** |
| 540 | selection | 1 | **0.0054** | **0.0065** | 3 | **0.0053** | **0.0049** |
| 572 | selection | absent |  |  | 2 | **0.0053** | **0.0049** |
| 733 | parallel | absent |  |  | absent |  |  |
| 761 | parallel | 12 | **0.0054** | **0.0044** | 14 | **0.0053** | **0.0052** |
| 749 | weak | 32 | **0.048** | 0.20 | 33 | 0.24 | 0.29 |
| 505 | weak | 5 | **0.0054** | **0.0044** | 7 | **0.0053** | **0.0052** |
| 573 | weak | absent |  |  | 11 | **0.0053** | **0.016** |
| 731 | weak | 7 | **0.0054** | **0.0044** | 5 | **0.0053** | **0.0049** |

| criterion | genotypic | phenotypic |
|---|---|---|
| truth positions present in the candidate set | 7/10 | 9/10 |
| truth positions with `p.adj_bh` below 0.05 | 7/10 (505, 540, 665, 731, 749, 761, 780) | 8/10 (all but 749 and 733) |
| truth positions with `p.adj_bh_fact` below 0.05 | 6/10 (505, 540, 665, 731, 761, 780) | 8/10 (same eight) |
| truth positions with `p.emp` below 0.05 | 7/10 | 8/10 |
| truth positions among the 10 highest `CAAS_score` | 5 | 6 |
| positions at the floor of `p.emp` (one shared `p.adj_bh`) | 14 | 13 |
| non-truth positions with `p.adj_bh` below 0.05 | 19 of 28 candidates | 11 of 25 candidates |
| non-truth positions with `p.adj_bh_fact` below 0.05 | 14 | 8 |

- **Recovery under the non-circular trait.** Eight of the nine phenotypic truth positions that enter the candidate set are exceeded by none of the 1000 null cycles and have `p.adj_bh` 0.0053 and `p.adj_bh_fact` 0.0049 to 0.016 (573 the largest); 749 does not (0.24 and 0.29). The scores rank 572, 540, 665 and 731 at 2 to 5 (0.814 to 0.693), 505 at 7 and 780 at 8 (0.609); the non-truth 518 ranks first (0.873) and 611 sixth.
- **The genotypic hits are partly circular.** 780 defines the genotypic label; 540 and 665 carry the derived residue in 23/23 genotypic-C4 tips (`pepc_results.md`, caveats). Their genotypic p-values are close to definitional and are not evidence for the method. 505, 731 and 761 are not definitional, but 505 and 731 belong to the weak tier of the truth set, which comes from the same genotype-defined branches.
- **The weak and parallel tiers** separate from the background in both runs for 505, 731 and 761 (at the floor in both), and for 573 under the phenotypic trait; 749 does not separate under the phenotypic trait (0.24) and sits at the threshold under the genotypic one (0.048).
- **Specificity.** Most candidates pass: in the genotypic run 19 of the 28 non-truth candidates are below 0.05 under `p.adj_bh` and 11 of 25 in the phenotypic run, so a call is a statement about signal in the labelled data and not about membership of the truth set (§5).
- **Limit.** The top-10 cutoff is post hoc and descriptive. With 1000 cycles the floor of `p.emp` is 0.001; the positions at the floor share one `p.adj_bh` value, so their order rests on `CAAS_score` and on the factorized p.

## 2b. Lenient and strict: where the truth set falls

Rank among the candidates of each run (rank / N, score); `no` = not a candidate.

| position | change | tier | geno, lenient (N=35) | geno, strict (N=10) | pheno, lenient (N=34) | pheno, strict (N=17) |
|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | 6 · 0.59 | 4 · 0.59 | 8 · 0.61 | 5 · 0.44 |
| 665 | H→N | mutagenesis | 1 · 1.00 | 1 · 1.00 | 4 · 0.78 | 4 · 0.51 |
| 540 | P→T | selection | 1 · 1.00 | 1 · 1.00 | 3 · 0.81 | 2 · 0.55 |
| 572 | E→Q | selection | no (cluster filter) | 3 · 0.75 | 2 · 0.81 | 3 · 0.55 |
| 733 | F→V | parallel | no | no | no | no |
| 761 | S→A | parallel | 12 · 0.29 | no | 14 · 0.29 | no |
| 749 | L→T | weak | 32 · 0.02 | no | 33 · 0.01 | no |
| 505 | F→L | weak | 5 · 0.67 | no | 7 · 0.67 | 8 · 0.35 |
| 573 | A→N | weak | no (cluster filter) | no | 11 · 0.36 | no |
| 731 | I→V | weak | 7 · 0.56 | 5 · 0.56 | 5 · 0.69 | 5 · 0.44 |

| | geno, lenient | geno, strict | pheno, lenient | pheno, strict |
|---|---|---|---|---|
| truth positions present | 7 | 5 | 9 | 6 |
| ranks | 1, 1, 5, 6, 7, 12, 32 | 1, 1, 3, 4, 5 | 2, 3, 4, 5, 7, 8, 11, 14, 33 | 2, 3, 4, 5, 5, 8 |
| median rank (percentile of N) | 6 (17 %) | 3 (30 %) | 7 (21 %) | 4.5 (26 %) |
| in the top 5 / top 10 | 3 / 5 | 5 / 5 | 4 / 6 | 5 / 6 |
| `p.adj_bh` below 0.05 | 7 | 5 | 8 | 6 |

- **Why a truth position is missing.** 733 carries no C4-specific residue. 572 and 573 under the genotypic trait are discovered in lenient mode and removed by the cluster filter (`filter_minlen` 3, `filter_maxcaas` 0.7: they sit in a group of 6 CAAS in 8 columns); strict mode forms no cluster (0 positions discarded) and keeps 572. 505, 749 and 761 (genotypic) and 573, 749 and 761 (phenotypic) keep 1 to 2 conserved pairs in lenient mode and are not called in strict mode.
- **Rank and percentile.** The strict runs have fewer candidates, so the rank improves without the positions moving up: 780 goes from the 17th to the 40th percentile in the genotypic run. The comparison of ranks across modes is not a comparison of support.
- **Mutagenesis and selection tiers** are present in the four runs except 572 (genotypic, lenient). The weak tier is the one that loses members in strict mode: 505, 573 and 749 disappear and 731 stays.

## 3. The four discordant tips

Residues at the truth positions, read from the alignment:

| tip | `c4` | `c4_phenotypic` | 780 | 665 | 540 | 572 | 761 | 749 | 505 | 573 | 731 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| *F. dichotoma* FM208033 | 1 | 1 | – | N | T | Q | – | – | L | N | – |
| *F. dichotoma* **FM208032** | 0 | **1** | – | H | P | E | – | – | F | A | I |
| *F. ferruginea* FM208035 | 1 | 1 | S | N | T | Q | A | T | L | N | V |
| *F. ferruginea* **FM208034** | 0 | **1** | A | H | P | E | S | L | F | A | I |
| *F. littoralis* FM208037 | 1 | 1 | S | N | T | Q | A | T | L | N | V |
| *F. littoralis* **FM208036** | 0 | **1** | A | H | P | E | S | L | F | A | I |
| ***Cyperus distichus*** | 1 | **0** | S | N | T | Q | A | T | L | N | V |

Each C4 *Fimbristylis* species contributes two accessions: one with the C4-type residue at every truth site, one with the ancestral residue at every truth site. Besnard et al. (2009, Results) report two distinct ppc-1 gene clusters in C4 *Fimbristylis*, i.e. a duplication predating the clade, and state that at least one Ser780-encoding ppc-1 copy was isolated from every C4 species analysed. On this evidence FM208032/034/036 are the non-C4 paralog of a C4 species (inferred from residues and the duplication statement; copy identity was not checked against their gene tree). The phenotypic annotation assigns the species phenotype to every gene copy (Morel et al. 2024), so these paralogs are labelled C4. *Cyperus distichus* (*Volkiella disticha* in Bruhl & Wilson 2007) is the reverse case: a full C4-type sequence in a species classified C3. Its cause (misclassification, lookup, or an unexpressed C4-type copy) was not tested.

At 780 these tips account for every non-concordant residue under the phenotypic trait: the only phenotypic-C3 tip with S is *C. distichus*, and the phenotypic-C4 tips with A are FM208034/036 (FM208032/033 are gapped).

### Hypotheses that do not detect a truth position

| pos | trait | hypotheses not detecting | of which contain a member gapped at the position | of which contain a gapped member or a discordant tip |
|---|---|---|---|---|
| 780 | geno | 39 | 39 | 39 |
| 780 | pheno | 27 | 26 | 27 |
| 665 | pheno | 8 | 0 | 8 (discordant) |
| 540 | pheno | 3 | 0 | 3 (discordant) |
| 572 | pheno | 3 | 0 | 3 (discordant) |
| 731 | geno | 23 | 23 | 23 |
| 731 | pheno | 15 | 12 | 15 |

Under the genotypic trait, 665, 540 and 505 are detected in all 100 hypotheses; every non-detection of 780 and 731 is explained by a gapped pair member, and no detecting hypothesis contains one. Under the phenotypic trait every non-detection of these positions is explained by a gap or a discordant tip. 42/100 phenotypic hypotheses contain ≥ 1 discordant tip (*C. distichus* 7, FM208032 14, FM208034 12, FM208036 12). With 3 pairs per hypothesis and `min_contrasts = 3`, one non-divergent pair removes the position from that hypothesis.

## 4. Contrast-pair structure changes with the trait

| | genotypic | phenotypic |
|---|---|---|
| pairs per hypothesis | 4 | 3 |
| distinct pairs | 129 | 130 |
| same-species pair-instances (C4 vs non-C4 ppc-1 copy of one species) | 27/400, in 24 hypotheses | 0/300 |
| cross-genus pair-instances | 55/400 (14 %) | 138/300 (46 %) |

Under the genotypic trait, within-species paralog pairs are available and selected; each differs at every C4-associated site by construction of the label. Under the phenotypic trait the *Eleocharis* accessions are pruned and the *Fimbristylis* paralogs become C4–C4, so the selector draws on more distant pairs. The null scales of the two runs are similar (§1), so the difference in pair structure does not show up as a difference in the null.

## 5. Non-truth positions

| | genotypic | phenotypic |
|---|---|---|
| `p.adj_bh < 0.05` | 19 positions (491, 501, 503, 509, 511, 517, 518, 527, 584, 590, 611, 616, 620, 744, 751, 839, 862, 879, 880) | 11 (474, 501, 503, 509, 518, 579, 588, 611, 620, 744, 751) |
| `p.adj_bh_fact < 0.05` | 14 (501, 509, 511, 517, 518, 527, 584, 611, 616, 620, 744, 751, 839, 879) | 8 (501, 518, 579, 588, 611, 620, 744, 751) |
| FADE BF ≥ 100 as well (either adjustment below 0.05) | 509, 517, 611, 839 | 474, 518, 611, 620 |

518 ranks first in the phenotypic run (0.873, all 100 hypotheses, `p.emp` at the floor, `p.adj_bh` 0.0053, `p.adj_bh_fact` 0.0049, FADE BF 6.9 × 10⁴); under the genotypic trait its score is 0.117 in 24 hypotheses (0.032 and 0.049). The difference comes from pruning the five L-carrying *Eleocharis* accessions (`pepc_results.md`, Method 2). 611 ranks third in the genotypic run (0.826) and sixth in the phenotypic run (0.689), at the floor in both and with a FADE call in both (BF 7.6 × 10³ and 1.2 × 10³). 620 passes both adjustments in both runs and has a FADE call in the phenotypic run (BF 1.2 × 10³). 751 is at the floor in both runs and has no FADE call. 474 and 460 are detected in 1 and 2 hypotheses and have scores of 0.005 and 0.016 (474: `p.emp` 0.011, `p.adj_bh` 0.040, `p.adj_bh_fact` 0.129, FADE BF 428). A hypothesis that does not detect a position contributes zero to the species it represents, so positions detected in one or two hypotheses score near 0; a small `p.emp` still follows when the null rarely detects them.

In strict mode the non-truth candidates are 5 (genotypic: 579, 584, 611, 620, 623) and 11 (phenotypic: 460, 509, 518, 579, 584, 588, 611, 620, 744, 813, 855), and 5 and 10 are below 0.05 under `p.adj_bh` (all but 813). FADE calls 579 and 611 (genotypic) and 518, 611 and 620 (phenotypic) among them. 518 ranks first in the phenotypic run in both modes (0.873 lenient, 0.584 strict; 65 hypotheses in strict mode against 100). The strict runs keep 620 and 611 in the genotypic top 7 with lower scores (0.535 and 0.324, against 0.690 and 0.826).

## 6. Multiple-testing family

### 6a. Why candidate-set p-values are sub-uniform

`p.emp = (k + 1)/(N + 1)` counts null cycles that **detect** the position **and** reach the observed score. Its implicit test statistic is T = `CAAS_score` if detected, −∞ otherwise. Over all positions this is a valid permutation p-value; a position never detected in the observed data has T = −∞ and p = 1. A detected position with a score of 0 takes p = 1 by rule.

The candidate set keeps only positions the observed data detected, i.e. it filters on T. For a truly null position that passes this filter and has a positive score, `p.emp` ≈ d × U, where d is the fraction of null cycles that detect the position, so `p.emp` cannot exceed ≈ d. Across candidate positions d has median 0.32 (genotypic, range 0.05 to 0.96) and 0.33 (phenotypic, 0.05 to 0.83). This is the non-independent filtering case of Bourgon et al. (2010): any FDR procedure applied only to the survivors is anti-conservative.

### 6b. Corrections compared

Counts of positions below 0.1, except the raw count (below 0.05). "Candidate set" = observed-detected positions (35 / 34). "Null universe" = candidate set ∪ every position detected in ≥ 1 null cycle (80 / 71), with p = 1 for positions the observed data did not detect.

| procedure | family | genotypic | phenotypic |
|---|---|---|---|
| raw `p.emp` < 0.05 | none | 28 | 26 |
| BH of `p.emp` | candidate set | 31 (min 0.002) | 30 (min 0.003) |
| Storey, λ = 0.5 | candidate set | 35 of 35 (π̂₀ = 0.00) | 34 of 34 (π̂₀ = 0.00) |
| **BH of `p.emp` (`p.adj_bh`, pipeline)** | **null universe** | **27** (min 0.0054) | **23** (min 0.0053) |
| **BH of `p.emp_fact` (`p.adj_bh_fact`, pipeline)** | **null universe** | **27** (min 0.0044) | **19** (min 0.0049) |

Strict mode (null universe 25 and 47 positions, 15 and 30 of them null-only):

| procedure | family | genotypic | phenotypic |
|---|---|---|---|
| raw `p.emp` < 0.05 | none | 10 | 17 |
| BH of `p.emp` | candidate set | 10 (min 0.001) | 17 (min 0.002) |
| Storey, λ = 0.5 | candidate set | 10 of 10 (π̂₀ = 0.00) | 17 of 17 (π̂₀ = 0.00) |
| **BH of `p.emp` (`p.adj_bh`, pipeline)** | **null universe** | **10** (min 0.0027) | **17** (min 0.0051) |
| **BH of `p.emp_fact` (`p.adj_bh_fact`, pipeline)** | **null universe** | **10** (min 0.0027) | **17** (min 0.0051) |

In strict mode the candidates are almost all at the floor and the families have a handful of positions, so the procedures coincide.

- **Storey on the candidate set is invalid, not merely unstable.** π̂₀ = #{p > λ}/(m(1 − λ)) assumes null p-values fill [0, 1]; §6a caps them near d, so π̂₀ is biased low (0.00 in both runs: no candidate has p > 0.5, the largest `p.emp` being 0.23 and 0.48) and Storey calls every candidate (35 and 34) where BH over the null universe calls 27 and 23.
- **Candidate-set BH** counts 31 and 30 against 27 and 23, with smaller adjusted values because m is 35 and 34 instead of 80 and 71.
- **On the null universe**, restoring the p = 1 mass makes π̂₀ ≈ 1, and Storey reduces to BH.

### 6c. Factorized p and calibration

`p.emp_fact = (nd + 1)/(N + 1) × (1 + #{detections of the class with score ≥ s})/(1 + #{detections of the class})`, where `nd` is the number of null cycles that score the position and the class is one of 20 percentile classes of the null detections (`FACT_PROP_CLASSES`). The evaluated position counts as one more detection of its own and is classed by `nd + 1`. It is not bounded below by 1/(N + 1). `p.adj_bh_fact` is its BH over the null universe.

Calibration: each of the 1000 null cycles in turn is treated as observed and scored against the remaining 999, with the `.fact_*` functions of `scoring_compute.R` (`pepc_null_calibration.R`). A position the cycle does not score has p = 1. Under exchangeability the share of (position, cycle) pairs with p ≤ α cannot exceed α. For `p.emp` the share divided by α is 1.00, 0.97 and 0.92 (genotypic) and 0.00, 0.88 and 0.89 (phenotypic) at α = 0.001, 0.01 and 0.05 (0 at 0.001 in the phenotypic run because the smallest p of a pseudo-observed cycle, 1/999, exceeds 0.001). For `p.emp_fact` the ratios are 1.13, 0.99 and 0.91 (genotypic) and 1.10, 1.07 and 0.94 (phenotypic). The share of null cycles, taken as observed, with at least one position below 0.05 and 0.1 after BH of `p.emp` is 0.006 and 0.077 (genotypic) and 0.015 and 0.064 (phenotypic); after BH of `p.emp_fact` it is 0.053 and 0.100 (genotypic) and 0.047 and 0.096 (phenotypic). `pepc_null_calibration.R` writes the share tables and the per-cycle rates.

- **Both p-values are at or below nominal.** `p.emp` is at or under the nominal level at every α (ratio 0.88 to 1.00), as expected for a statistic that counts ties and detection at the same time, and BH of `p.emp` calls in 0.6 % and 1.5 % of the null cycles at 0.05. `p.emp_fact` is within 14 % of nominal at the three levels, and BH of `p.emp_fact` calls in 5.3 % and 4.7 % of the cycles at 0.05 and in 10.0 % and 9.6 % at 0.1, i.e. at the nominal levels.
- **Strict mode** is conservative: `p.emp` is 1.00, 0.73 and 0.40 of nominal at α = 0.001, 0.01 and 0.05 (genotypic) and 1.00, 0.90 and 0.70 (phenotypic), and `p.emp_fact` 0.92, 0.66 and 0.40 and 0.94, 0.88 and 0.69. Taken as observed, a null cycle has at least one position with BH < 0.05 / < 0.1 in 2.7 % / 6.7 % (genotypic) and 4.4 % / 7.8 % (phenotypic) of the cycles under `p.emp`, and 3.6 % / 7.5 % and 4.1 % / 8.7 % under `p.emp_fact`. The null of strict mode holds 513 (genotypic) and 3 200 (phenotypic) position-cycle detections over 1000 cycles, so most cycles detect nothing and the p-value of a position rests on few cycles.
- **Scope.** Null cycles are draws from one generator, so this checks the estimators under exchangeability, not whether the observed labelling is exchangeable with the null ones (§6d). Pairs of the same position or the same cycle are not independent, and no standard error is given for the shares. Behavior under partial nulls and at genome scale is not tested.

### 6d. Is the observed labelling exchangeable with the null?

| check | genotypic | phenotypic |
|---|---|---|
| observed detections vs null per-cycle detections | 35; null median 25, q95 35; P(null ≥ obs) = 0.065 | 34; null median 19, q95 35; P = 0.064 |
| observed vs pooled-null `CAAS_score` quartiles | 0.056 / 0.117 / 0.421 vs 0.011 / 0.026 / 0.056 | 0.031 / 0.122 / 0.447 vs 0.013 / 0.029 / 0.065 |
| strict: observed detections vs null per-cycle detections | 10; null median 0, q95 2; P(null ≥ obs) = 0.005 | 17; null median 1, q95 14; P = 0.020 |
| strict: observed vs pooled-null `CAAS_score` quartiles | 0.266 / 0.548 / 0.711 vs 0.009 / 0.013 / 0.022 | 0.070 / 0.250 / 0.438 vs 0.009 / 0.015 / 0.030 |
| null cycles with the observed design (100 hypotheses) | 1000/1000 | 1000/1000 |

The observed data detect as many positions as the upper tail of the null cycles (P(null ≥ observed) 0.064 to 0.065), and score above their median; this is what signal would produce and is not by itself diagnostic. Because every null cycle carries the same 100-hypothesis design as the observed data, the null is conditioned on design size by construction. The real labelling, routed through the null path as cycle b_0, reproduces every position score of the observed run (maximum difference 0 over 36 and 34 position-side scores in lenient mode and 0 and 1 × 10⁻¹⁶ over 10 and 17 in strict mode), so the scoring is shared between the two paths. Whether the full observed path is exchangeable with the null under no association is tested by the negative controls (§6e).

### 6e. Negative controls

Negative controls (control traits drawn from the permutation generator and analyzed through the complete observed workflow) have not been run with the current pipeline. The check available in this report is the null-cycle calibration of §6c, which tests the estimators under exchangeability and not whether the observed labelling is exchangeable with the null ones.

## 7. Sanity checks on the phenotypic input

| check | result |
|---|---|
| trait counts (in-group) | 20 C4 / 50 C3, 7 pruned (NA); matches `input/.pepc_phenotypic/README.md` |
| pruned *Eleocharis* tips in any contrast pair | 0 (the genotypic run uses all 7) |
| pair members with missing trait | 0/300 phenotypic, 0/400 genotypic |
| pair label structure | 300/300 C4–C3 (`abs_diff = 1`); genotypic 400/400 |
| DAG | 33/33 (genotypic) and 34/34 (phenotypic) tasks succeeded; the runs differ only in the pruning path (`pepc_resources.md`) |

The discordant tips behave as their residues predict and account for every phenotypic non-detection of the truth positions not explained by gaps (§3).

The 2-tip count discrepancy (this fixture 20 C4 vs Morel et al.'s 22 on the same 71 tips) cannot be resolved here. One bound follows from the alignment: *C. distichus* is the only phenotypic-C3 tip carrying Ser780, so any two further tips relabelled C4 would carry Ala or a gap at 780 and add label noise in the same direction as the *Fimbristylis* paralogs.

## 8. Relation to Morel et al. (2024)

On the same sedge PEPC data, Morel et al. (2024, Table 2) report PCOC recovering 7 of 12 convergent mutations under the genotypic annotation and none of 11 under the phenotypic one, while ConDor retains the best phenotypic precision (0.57). Under the phenotypic annotation PhyloPhere CAAS gives both mutagenesis sites and both selection-tier sites an adjusted p of 0.0053 or lower under both adjustments. Morel et al. report FADE, with all branches of convergent clades as foreground, as the best method on both annotations (F₁ 0.81 genotypic, 0.57 phenotypic); PhyloPhere's FADE, run with the same foreground definition, recovers 8/10 and 7/10 of this truth set. The comparison is qualitative only: Morel et al. use a 458-column alignment with their own convergent-mutation set (12 genotypic / 11 phenotypic mutations; e.g. 749 recorded as M→T, where this truth set has L→T), so counts are not transferable.

## 9. Limits

- One run per trait, one seed (1998); `p.emp` floor 0.001.
- The truth set is not phenotype-independent: selection- and weak-tier sites come from Besnard et al. (2009) tests on branches defined by Ser780. Only the mutagenesis tier (780, 665) and the parallel tier (733, 761) carry independent evidence, and 733 is uninformative in this sample.
- The phenotypic trait is not clean either: at gene level it mislabels non-C4 paralogs of C4 species. A gene-copy-aware phenotypic label would require copy assignment, which in this dataset is made by the Ser780 residue itself; the two concerns cannot be fully separated with this fixture.
- n = 10 truth positions in one gene.
- Negative controls with the current pipeline are pending (§6e).
- The strict runs have 10 and 17 candidates and 5 and 6 reference positions; their percentiles and shares (50 % and 35 % of the candidates are reference positions) describe one gene and carry no interval.
- The value of `min_divergent_fraction` of the strict runs is not recorded; the strict mode is identified by `n_conserved_pairs` = 0 in every position.
- The positions at the floor of `p.emp` (14 and 13) share one adjusted value, so the order among them rests on `CAAS_score` and on the factorized p; with 1000 cycles `p.emp` cannot separate them.

## References

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. Proc Natl Acad Sci USA 100(16):9440–9445. doi:10.1073/pnas.1530509100.
