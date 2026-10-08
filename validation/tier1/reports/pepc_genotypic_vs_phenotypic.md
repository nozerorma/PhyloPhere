# Tier 1 PEPC: Genotypic vs Phenotypic Trait

Two PhyloPhere runs on the same 970-column ppc-1 alignment and tree, differing in the trait. Per-run tables, the FUBAR baseline, FADE (including its agreement with CAAS) and numbering verification are in `pepc_results.md`; runtime and run configuration in `pepc_resources.md`.

| | genotypic | phenotypic |
|---|---|---|
| results dir | `output/pepc/results/c4_complete/` | `output/pepc/results/c4_phenotypic_complete/` |
| trait column (`input/pepc/my_traits.tsv`) | `c4` | `c4_phenotypic` |
| definition | ppc-1 copy carries Ser780 (Besnard et al. 2009) | species-level anatomical/physiological call (Bruhl & Wilson 2007), assigned to every ppc-1 copy of that species (as in Morel et al. 2024) |
| tips (in-group) | 77 (23 C4 / 54 C3) | 70 (20 C4 / 50 C3); 7 *E. baldwinii* / *E. vivipara* C3/C4-intermediate accessions pruned |
| evolutionary model (AIC) | BM | BM |
| candidate positions (unique) | 59 | 57 |
| hypotheses × pairs | 100 × 4 | 100 × 3 |
| null cycles, hypotheses per cycle | 1000, 100 each | 1000, 100 each |

Position numbers are maize PEPC1 (`position_scores.tsv` `Position` + 1). `p.adj_bh` is BH of `p.emp` over the null universe (one test per position) and `p.adj_bh_fact` the BH of the factorized `p.emp_fact` over the same family (§6c); the pipeline threshold is 0.05 for both. `CAAS_score` is the sum of the per-scheme scores over the five schemes divided by five; a scheme that does not detect the position counts 0. `pepc_pvalue_tables.py` prints the tables of this report from the results directory.

## 1. Truth positions against their nulls

`k_emp` = null cycles that re-detect the position and reach ≥ the observed max-over-sides `CAAS_score`. Null quantiles are over the cycles that detect the position. Recomputed from `caas_permulation/perm_pos_cycle_caas.tsv.gz` with the §2f-ter definition; `k_emp` reproduces the reported `p.emp` in every case.

| pos | trait | observed CAAS | n_hyp /100 | null detects | null q50 | null q90 | null q99 | k_emp | p.emp | p.adj_bh |
|---|---|---|---|---|---|---|---|---|---|---|
| 780 | geno | 0.600 | 49 | 862 | 0.232 | 0.342 | 0.416 | 0 | 0.001 | **0.026** |
| 780 | pheno | 0.549 | 73 | 632 | 0.316 | 0.415 | 0.459 | 2 | 0.003 | **0.041** |
| 665 | geno | 0.800 | 100 | 914 | 0.305 | 0.458 | 0.575 | 0 | 0.001 | **0.026** |
| 665 | pheno | 0.730 | 92 | 702 | 0.420 | 0.555 | 0.620 | 1 | 0.002 | **0.038** |
| 540 | geno | 0.800 | 100 | 913 | 0.276 | 0.440 | 0.549 | 0 | 0.001 | **0.026** |
| 540 | pheno | 0.725 | 97 | 817 | 0.349 | 0.524 | 0.601 | 0 | 0.001 | **0.038** |
| 572 | geno | 0.283 | 100 | 900 | 0.121 | 0.274 | 0.351 | 79 | 0.080 | 0.339 |
| 572 | pheno | 0.543 | 97 | 671 | 0.313 | 0.416 | 0.462 | 1 | 0.002 | **0.038** |
| 731 | geno | 0.422 | 76 | 733 | 0.143 | 0.286 | 0.359 | 2 | 0.003 | 0.053 |
| 731 | pheno | 0.536 | 85 | 655 | 0.312 | 0.414 | 0.457 | 1 | 0.002 | **0.038** |
| 505 | geno | 0.366 | 100 | 741 | 0.212 | 0.320 | 0.380 | 19 | 0.020 | 0.235 |
| 505 | pheno | 0.479 | 84 | 667 | 0.325 | 0.422 | 0.487 | 10 | 0.011 | 0.095 |

- **Under the phenotypic trait the five strongest truth positions (780, 665, 540, 572, 731) sit above their own null q99**, with `k_emp` ≤ 2 of 1000. Under the genotypic trait four do (780, 665, 540, 731); 572 (0.283) is below its null q99 (0.351) and has `k_emp` = 79. At 540 no null cycle reaches the observed score under either trait.
- **Null scale differs between runs.** Null medians at the truth positions are 0.12 to 0.31 (genotypic) and 0.31 to 0.42 (phenotypic); the pooled null quartiles over all detected positions are 0 / 0.115 / 0.234 and 0 / 0.183 / 0.352, and the observed quartiles 0.078 / 0.217 / 0.326 and 0.062 / 0.277 / 0.469. `CAAS_score` is therefore not comparable across the two runs; `p.emp`, calibrated within each run, is.
- The phenotypic null detects each truth position less often (632 to 817 cycles vs 733 to 914), consistent with 3-of-3 all-or-nothing detection per hypothesis (§4).

## 2. Truth set, side by side

| pos | tier | geno score rank /59 | geno p.adj_bh | geno p.adj_bh_fact | pheno score rank /57 | pheno p.adj_bh | pheno p.adj_bh_fact |
|---|---|---|---|---|---|---|---|
| 780 | mutagenesis | 4 | **0.026** | **0.035** | 7 | **0.041** | **0.029** |
| 665 | mutagenesis | 1 | **0.026** | **0.027** | 1 | **0.038** | **0.029** |
| 540 | selection | 1 | **0.026** | **0.027** | 2 | **0.038** | **0.029** |
| 572 | selection | 20 | 0.339 | 0.470 | 8 | **0.038** | **0.029** |
| 733 | parallel | absent | | | absent | | |
| 761 | parallel | 27 | 0.402 | 0.309 | 23 | 0.480 | 0.276 |
| 749 | weak | 18 | 0.284 | 0.368 | 25 | 0.427 | 0.635 |
| 505 | weak | 11 | 0.235 | 0.354 | 14 | 0.095 | 0.089 |
| 573 | weak | 17 | 0.482 | 0.354 | 12 | 0.553 | 0.285 |
| 731 | weak | 10 | 0.053 | 0.309 | 10 | **0.038** | 0.198 |

| criterion | genotypic | phenotypic |
|---|---|---|
| truth positions with `p.adj_bh < 0.05` | 3/10 (780, 665, 540) | 5/10 (780, 665, 540, 572, 731) |
| truth positions with `p.adj_bh_fact < 0.05` | 3/10 (780, 665, 540) | 4/10 (780, 665, 540, 572) |
| truth positions among the 10 highest `CAAS_score` | 4 | 5 |
| non-truth positions with `p.adj_bh < 0.05` | 1 (751) | 2 (474, 518) |
| non-truth positions with `p.adj_bh_fact < 0.05` | 2 (611, 620) | 1 (518) |

- **Calibrated recovery holds under the non-circular trait.** Both mutagenesis sites and both selection-tier sites are significant under the phenotypic trait with both adjustments. Their score ranks are 1, 2, 7 and 8: positions 662 and 620 (outside the truth set) rank 3 and 4 and are not significant after adjustment.
- **The genotypic hits are partly circular.** 780 defines the genotypic label; 540 and 665 carry the derived residue in 23/23 genotypic-C4 tips (`pepc_results.md`, caveats). Their genotypic significance is close to definitional and is not evidence for the method. Under the genotypic trait 572 does not pass (`p.adj_bh` 0.339).
- **The weak tier does not separate from background** in either run, except 731 under the phenotypic trait with `p.adj_bh` (0.038; `p.adj_bh_fact` 0.198) and 505 under the phenotypic trait at 0.095 (0.089).
- **Limit.** The top-10 cutoff is post hoc and descriptive. Ranks among the top positions rest on `p.emp` values a few null cycles apart; with 1000 cycles the floor is 0.001.

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
| 780 | geno | 51 | 51 | 51 |
| 780 | pheno | 27 | 26 | 27 |
| 665 | pheno | 8 | | 8 (discordant) |
| 540 | pheno | 3 | | 3 (discordant) |
| 572 | pheno | 3 | | 3 (discordant) |
| 731 | geno | 24 | 24 | 24 |
| 731 | pheno | 15 | 12 | 15 |

Under the genotypic trait, 665, 540 and 572 are detected in all 100 hypotheses; every non-detection of 780 and 731 is explained by a gapped pair member (`max_*gaps_fraction = 0`), and no detecting hypothesis contains one. Under the phenotypic trait every non-detection is explained by a gap or a discordant tip. 42/100 phenotypic hypotheses contain ≥ 1 discordant tip (*C. distichus* 7, FM208032 14, FM208034 12, FM208036 12). With 3 pairs per hypothesis and `min_contrasts = 3`, one non-divergent pair removes the position from that hypothesis.

## 4. Contrast-pair structure changes with the trait

| | genotypic | phenotypic |
|---|---|---|
| pairs per hypothesis | 4 | 3 |
| distinct pairs | 126 | 130 |
| same-species pair-instances (C4 vs non-C4 ppc-1 copy of one species) | 23/400, in 21 hypotheses | 0/300 |
| cross-genus pair-instances | 50/400 (13 %) | 138/300 (46 %) |

Under the genotypic trait, within-species paralog pairs are available and selected; each differs at every C4-associated site by construction of the label. Under the phenotypic trait the *Eleocharis* accessions are pruned and the *Fimbristylis* paralogs become C4–C4, so the selector draws on more distant pairs. This coincides with the higher null scale in §1: 3-pair all-or-nothing detection and pairs spanning more genus boundaries could each raise the per-cycle score of whatever a random relabelling detects. Neither was tested in isolation.

## 5. Non-truth positions

| | genotypic | phenotypic |
|---|---|---|
| `p.adj_bh < 0.05` | 751 | 474, 518 |
| `0.05 ≤ p.adj_bh < 0.1` | 611, 620 | 460, 501, 620 |
| `p.adj_bh_fact < 0.05` | 611, 620 | 518 |
| `0.05 ≤ p.adj_bh_fact < 0.1` | none | 474, 620 |
| FADE BF ≥ 100 as well (`p.adj_bh < 0.1`) | 611 | 474, 518, 620 |

611 is present in both runs (genotypic `p.adj_bh` 0.061, `p.adj_bh_fact` 0.035; phenotypic 0.293 and 0.507) and is FADE-significant in both (BF 7.6 × 10³ and 1.2 × 10³). 620 ranks 3 (genotypic) and 4 (phenotypic) by score, and 751 is the only genotypic call of `p.adj_bh` that no FADE call and no `p.adj_bh_fact` supports (0.216). 518 passes both adjustments only under the phenotypic trait (0.038 and 0.029); the difference comes from pruning the five L-carrying *Eleocharis* accessions (`pepc_results.md`, Method 2). 460 and 474 are detected in 2 and 1 hypotheses with low tip coverage; hypothesis recurrence and coverage do not enter `CAAS_score` or `p.emp`, so such positions can rank alongside positions detected in every hypothesis.

## 6. Multiple-testing family

### 6a. Why candidate-set p-values are sub-uniform

`p.emp = (k + 1)/(N + 1)` counts null cycles that **detect** the position **and** reach the observed score. Its implicit test statistic is T = `CAAS_score` if detected, −∞ otherwise. Over all positions this is a valid permutation p-value; a position never detected in the observed data has T = −∞ and p = 1. A detected position with a score of 0 takes p = 1 by rule.

The candidate set keeps only positions the observed data detected, i.e. it filters on T. For a truly null position that passes this filter and has a positive score, `p.emp` ≈ d × U, where d is the fraction of null cycles that detect the position, so `p.emp` cannot exceed ≈ d. Across candidate positions d has median 0.51 (genotypic, range 0.08 to 0.96) and 0.40 (phenotypic, 0.04 to 0.88). This is the non-independent filtering case of Bourgon et al. (2010): any FDR procedure applied only to the survivors is anti-conservative.

### 6b. Corrections compared

Counts of positions below 0.1, except the raw count. "Candidate set" = observed-detected positions (59 / 57). "Null universe" = candidate set ∪ every position detected in ≥ 1 null cycle (106 / 95), with p = 1 for positions the observed data did not detect.

| procedure | family | genotypic | phenotypic |
|---|---|---|---|
| raw `p.emp` < 0.05 | none | 18 | 18 |
| BH of `p.emp` | candidate set | 7 (min 0.015) | 11 (min 0.023) |
| Storey, λ = 0.5 | candidate set | 34 of 59 (π̂₀ = 0.44) | 24 of 57 (π̂₀ = 0.49) |
| **BH of `p.emp` (`p.adj_bh`, pipeline)** | **null universe** | **7** (min 0.026) | **11** (min 0.038) |
| **BH of `p.emp_fact` (`p.adj_bh_fact`, pipeline)** | **null universe** | **5** (min 0.027) | **8** (min 0.029) |

- **Storey on the candidate set is invalid, not merely unstable.** π̂₀ = #{p > λ}/(m(1 − λ)) assumes null p-values fill [0, 1]; §6a caps them near d, so π̂₀ is biased low (0.44 and 0.49) and Storey calls 34 and 24 positions where BH over the null universe calls 7 and 11.
- **Candidate-set BH** gives the same counts as the null universe here, with smaller adjusted values because m is 59 and 57 instead of 106 and 95. The agreement is a property of these two runs: the family only matters when positions near the floor are few.
- **On the null universe**, restoring the p = 1 mass makes π̂₀ ≈ 1, and Storey reduces to BH.

### 6c. Factorized p and calibration

`p.emp_fact = (nd + 1)/(N + 1) × (1 + #{detections of the class with score ≥ s})/(1 + #{detections of the class})`, where `nd` is the number of null cycles that score the position and the class is one of 20 percentile classes of the null detections (`FACT_PROP_CLASSES`). It is not bounded below by 1/(N + 1). `p.adj_bh_fact` is its BH over the null universe.

Calibration: each of the 1000 null cycles in turn is treated as observed and scored against the remaining 999, with the `.fact_*` functions of `scoring_compute.R` (`pepc_null_calibration.R`). A position the cycle does not score has p = 1. Under exchangeability the share of (position, cycle) pairs with p ≤ α cannot exceed α, overall or within a class. Entries are that share divided by α (≤ 1 is calibrated or conservative):

| p | trait | α = 0.001 | α = 0.01 | α = 0.05 |
|---|---|---|---|---|
| `p.emp` | genotypic | 0.84 | 0.83 | 0.79 |
| `p.emp_fact` | genotypic | 1.09 | 0.95 | 0.91 |
| `p.emp` | phenotypic | 0.77 | 0.79 | 0.76 |
| `p.emp_fact` | phenotypic | 0.92 | 0.94 | 0.89 |

By propensity class of the position (cycles that score it; the model's own 20 classes are finer), `p.emp_fact`, share / α at 0.001 / 0.01 / 0.05:

| nd | positions geno / pheno | genotypic | phenotypic |
|---|---|---|---|
| ≤ 5 | 1 / 2 | 3.0 / 0.30 / 0.06 | 1.5 / 0.15 / 0.03 |
| 6 to 20 | 6 / 6 | 3.7 / 0.85 / 0.17 | 3.8 / 0.72 / 0.14 |
| 21 to 100 | 12 / 12 | 1.2 / 1.10 / 0.82 | 0.08 / 0.83 / 0.80 |
| > 100 | 87 / 75 | 0.87 / 0.95 / 0.98 | 0.80 / 1.00 / 0.99 |

Share of null cycles, taken as observed, with at least one position at adjusted p < 0.05 / < 0.1: BH of `p.emp` 0.001 / 0.011 (genotypic) and 0.012 / 0.066 (phenotypic); BH of `p.emp_fact` 0.059 / 0.107 and 0.053 / 0.092 (Monte Carlo SE ≈ 0.007 at 0.05).

- **`p.emp` is conservative** at every α and in every class, as expected for a statistic that counts ties and detection at the same time.
- **`p.emp_fact` is calibrated overall, with violations in the extreme tail of small classes.** At α = 0.001 the classes nd ≤ 5 and 6 to 20 exceed α (1.5 to 3.8 times; 3 to 23 pairs out of 1000 to 6000), and the overall genotypic share is 1.09 times α (115 pairs where 106 are expected). At α = 0.01 and 0.05 the only entry above 1 is the genotypic class 21 to 100 at 0.01 (1.10). The classes with violations hold 1 to 6 positions, so the check has little resolution there.
- **BH of `p.emp_fact` is close to nominal** (0.059 and 0.053 at 0.05, 0.107 and 0.092 at 0.1), whereas BH of `p.emp` spends well under its budget, most visibly in the genotypic run where the 1/N floor binds.
- **Scope.** Null cycles are draws from one generator, so this checks the estimator under exchangeability, not whether the observed labelling is exchangeable with the null ones (§6d). Pairs of the same position or the same cycle are not independent, and no standard error is given for the shares. Behavior under partial nulls and at genome scale is not tested.

### 6d. Is the observed labelling exchangeable with the null?

| check | genotypic | phenotypic |
|---|---|---|
| observed detections vs null per-cycle detections | 59; null median 42, q95 57; P(null ≥ obs) = 0.036 | 57; null median 34, q95 50; P = 0.018 |
| observed vs pooled-null `CAAS_score` quartiles | 0.078 / 0.217 / 0.326 vs 0 / 0.115 / 0.234 | 0.062 / 0.277 / 0.469 vs 0 / 0.183 / 0.352 |
| null cycles with the observed design (100 hypotheses) | 1000/1000 | 1000/1000 |

The observed data detect more positions, at higher scores, than a typical null cycle; this is what signal would produce and is not by itself diagnostic. Because every null cycle carries the same 100-hypothesis design as the observed data, the null is conditioned on design size by construction. The real labelling, routed through the null path as cycle b_0, reproduces every position score of the observed run (maximum difference 1 × 10⁻¹⁶ over 67 and 71 position-side scores), so the scoring is shared between the two paths. Whether the full observed path is exchangeable with the null under no association is tested by the negative controls (§6e).

### 6e. Negative controls

These controls were run with the pipeline before the unified CAAS core and have not been repeated with it (mean aggregation of the scheme scores, no factorized p); their inputs and outputs are in `previous_work/`. Twenty control traits (`previous_work/input/pepc_negctrl/my_traits.tsv`, columns nc01 to nc20) were drawn by `previous_work/input/pepc_negctrl/build_negctrl_traits.R` from the permulation null's own generator fitted on the genotypic trait (BM, 23 foreground tips, Tier-1 Dunn acceptance at 4 pairs, seed 2026). Each was run through the observed path with the genotypic run's settings (`previous_work/input/pepc_negctrl/run_negctrl_local.sh`; FADE, enrichment and reports off). Overlap with the real C4 set is low (Jaccard ≤ 0.10; φ from −0.43 to 0.38). Summary: `previous_work/output/pepc_negctrl/negctrl_summary.tsv`, from `previous_work/input/pepc_negctrl/analyze_negctrl.py`.

| quantity | value |
|---|---|
| controls with zero observed discoveries (no candidate positions, hence no calls) | 6/20 |
| controls with candidate positions | 14/20; 1 to 18 positions each (median 2.5) |
| observed detections vs null per-cycle detections | above the null median in 8/14; median P(null ≥ obs) = 0.40 (real traits: 0.036, 0.018) |
| family p-values (null universe), pooled: P(p ≤ 0.01 / 0.05 / 0.10) | 0.028 / 0.069 / 0.093 |
| observed-detected p-values only, pooled: P(p ≤ 0.01 / 0.05 / 0.10) | 0.24 / 0.58 / 0.78 |
| controls with ≥ 1 call, BH `p.adj_bh < 0.05` | 1/20 |
| controls with ≥ 1 call, BH `p.adj_bh < 0.1` | 4/20 |

- **Detection counts behave like null draws.** Control labellings detect about as many positions as their own null cycles (median P = 0.40), while both real traits sit in the upper tail of theirs. This is the most direct check that the observed and null paths are exchangeable, and it holds.
- **Run-level false calls are near nominal.** Under a complete null, P(≥ 1 call) is the realised FDR. At 0.05 it is 1/20; at 0.1 it is 4/20 (binomial P(X ≥ 4 | n = 20, p = 0.1) = 0.13, so not distinguishable from nominal with 20 controls).
- **Filtering on detection makes candidate-set p-values anti-conservative, as §6a predicts.** Restricted to observed-detected positions, 58 % of p-values fall below 0.05.
- **The family p-values are anti-conservative in the extreme tail** (2.8 × nominal at 0.01, 1.4 × at 0.05, nominal at 0.10). The 19 control positions with p ≤ 0.01 are positions the null rarely detects (median 28 of 1000 cycles, against 195 for all observed control positions). For such positions the add-one p has little resolution, and a position the null never detects enters at p = 1/1001 while every other never-detected column stays outside the family. nc13 shows the case: its single observed position (maize 723) is detected in 0 of 1000 null cycles and its family has 17 positions, so BH gives 0.017.
- **Small candidate sets and the coordinate guard.** `scoring_compute.R` sets `p.emp` to NA when fewer than half of the observed positions are re-detected by the null, a check for coordinate mismatches between the observed and null tables; it applies only when ≥ 10 positions are observed, and below that an unmatched position is scored at k = 0. `analyze_negctrl.py` applies the same rule, so nc13's position enters the counts above at p = 1/1001.
- **The tail excess is not explained by rare positions or by the null model.** Restricting the family to positions detected in ≥ 10 null cycles removes nc13's call but does not restore calibration (P(p ≤ 0.01 / 0.05 / 0.10) = 0.029 / 0.079 / 0.108), and controls whose null selects OU and BM show the same excess (0.027 vs 0.030 at 0.01). The excess is spread across controls (11 of 14 have ≥ 1 position at p ≤ 0.01). Run-level calls stay near nominal because, with families of 17 to 65 positions, a BH call at 0.05 needs a `p.emp` at or near the 1/1001 floor.
- **Generator mismatch.** Controls were drawn from the BM fit on the genotypic trait, but each control's own null refits the model on that control and selects OU in 7 of the 14 controls with candidate positions. Control and null generators therefore differ in those runs. Where OU is rejected despite a lower AIC (e.g. nc01), `select_model` has fallen back to BM because the OU α optimised to its upper bound.

- **Limit.** The null universe is defined by this run's 1000 cycles; positions with very low detection rates may be missing from it, which understates *m*; nc13 in §6e is a case where this decides a call. Columns never detected in the null sample are excluded unless the observed data detect them.

## 7. Sanity checks on the phenotypic input

| check | result |
|---|---|
| trait counts (in-group) | 20 C4 / 50 C3, 7 pruned (NA); matches `input/.pepc_phenotypic/README.md` |
| pruned *Eleocharis* tips in any contrast pair | 0 (the genotypic run uses all 7) |
| pair members with missing trait | 0/300 phenotypic, 0/400 genotypic |
| pair label structure | 300/300 C4–C3 (`abs_diff = 1`); genotypic 400/400 |
| DAG | 37/37 tasks succeeded in both runs; the runs differ only in the pruning path (`pepc_resources.md`) |

The discordant tips behave as their residues predict and account for every phenotypic non-detection of the truth positions not explained by gaps (§3).

The 2-tip count discrepancy (this fixture 20 C4 vs Morel et al.'s 22 on the same 71 tips) cannot be resolved here. One bound follows from the alignment: *C. distichus* is the only phenotypic-C3 tip carrying Ser780, so any two further tips relabelled C4 would carry Ala or a gap at 780 and add label noise in the same direction as the *Fimbristylis* paralogs.

## 8. Relation to Morel et al. (2024)

On the same sedge PEPC data, Morel et al. (2024, Table 2) report PCOC recovering 7 of 12 convergent mutations under the genotypic annotation and none of 11 under the phenotypic one, while ConDor retains the best phenotypic precision (0.57). PhyloPhere CAAS keeps both mutagenesis sites and both selection-tier sites significant under the phenotypic annotation. Morel et al. report FADE, with all branches of convergent clades as foreground, as the best method on both annotations (F₁ 0.81 genotypic, 0.57 phenotypic); PhyloPhere's FADE, run with the same foreground definition, recovers 8/10 and 7/10 of this truth set. The comparison is qualitative only: Morel et al. use a 458-column alignment with their own convergent-mutation set (12 genotypic / 11 phenotypic mutations; e.g. 749 recorded as M→T, where this truth set has L→T), so counts are not transferable.

## 9. Limits

- One run per trait, one seed (1998); `p.emp` resolution ≈ 0.001, and ranks among the top positions turn on a few null cycles.
- The truth set is not phenotype-independent: selection- and weak-tier sites come from Besnard et al. (2009) tests on branches defined by Ser780. Only the mutagenesis tier (780, 665) and the parallel tier (733, 761) carry independent evidence, and 733 is uninformative in this sample.
- The phenotypic trait is not clean either: at gene level it mislabels non-C4 paralogs of C4 species. A gene-copy-aware phenotypic label would require copy assignment, which in this dataset is made by the Ser780 residue itself; the two concerns cannot be fully separated with this fixture.
- The mechanism proposed for the null-scale difference (§4) is a hypothesis.
- n = 10 truth positions in one gene.
- The calibration of `p.emp_fact` (§6c) rests on 95 and 106 scored positions; the classes where it exceeds α hold 1 to 6 positions.
- The negative controls (§6e) were not repeated with the factorized p.

## References

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. Proc Natl Acad Sci USA 100(16):9440–9445. doi:10.1073/pnas.1530509100.
