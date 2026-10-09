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

Position numbers are maize PEPC1 (`position_scores.tsv` `Position` + 1). The four position p-values are informative and carry no threshold in the pipeline: `p.emp`, its BH adjustment over the null universe `p.adj_bh` (one test per position), the factorized `p.emp_fact` and its BH adjustment `p.adj_bh_fact` over the same family (§6c). This report marks values below 0.05 and counts positions below 0.05 and 0.1 to describe them. `CAAS_score` is `(US + mean(GS)) / 2`: half the exact-identity score (US) plus half the mean of the grouping-scheme scores (GS1 to GS4) of the schemes that scored the position, each term 0 when absent; it lies between 0 and 1. `pepc_pvalue_tables.py` prints the tables of this report from the results directory.

## 1. Truth positions against their nulls

`k_emp` = null cycles that re-detect the position and reach ≥ the observed max-over-sides `CAAS_score`. Null quantiles are over the cycles that detect the position. Recomputed from `caas_permulation/perm_pos_cycle_caas.tsv.gz` with the §2f-ter definition; `k_emp` reproduces the reported `p.emp` in every case.

| pos | trait | observed CAAS | n_hyp /100 | null detects | null q50 | null q90 | null q99 | k_emp | p.emp | p.adj_bh |
|---|---|---|---|---|---|---|---|---|---|---|
| 780 | geno | 1.000 | 49 | 867 | 0.455 | 0.622 | 0.745 | 0 | <0.001 | **0.026** |
| 780 | pheno | 0.915 | 73 | 543 | 0.592 | 0.725 | 0.807 | 0 | <0.001 | **0.016** |
| 665 | geno | 1.000 | 100 | 936 | 0.451 | 0.621 | 0.713 | 0 | <0.001 | **0.026** |
| 665 | pheno | 0.912 | 92 | 670 | 0.601 | 0.726 | 0.797 | 0 | <0.001 | **0.016** |
| 540 | geno | 1.000 | 100 | 931 | 0.415 | 0.595 | 0.709 | 0 | <0.001 | **0.026** |
| 540 | pheno | 0.906 | 97 | 821 | 0.460 | 0.693 | 0.779 | 0 | <0.001 | **0.016** |
| 572 | geno | 0.353 | 100 | 907 | 0.236 | 0.504 | 0.633 | 275 | 0.276 | 0.75 |
| 572 | pheno | 0.904 | 97 | 576 | 0.597 | 0.725 | 0.785 | 0 | <0.001 | **0.016** |
| 731 | geno | 0.704 | 76 | 875 | 0.277 | 0.507 | 0.627 | 1 | 0.0020 | **0.042** |
| 731 | pheno | 0.893 | 85 | 563 | 0.595 | 0.725 | 0.780 | 0 | <0.001 | **0.016** |
| 505 | geno | 0.611 | 100 | 762 | 0.404 | 0.569 | 0.674 | 42 | 0.043 | 0.44 |
| 505 | pheno | 0.798 | 84 | 577 | 0.599 | 0.752 | 0.822 | 15 | 0.016 | 0.17 |

- **Under the phenotypic trait the five strongest truth positions (540, 572, 665, 731, 780) exceed all 1000 null cycles** (`k_emp` = 0, `p.emp` at its floor). Under the genotypic trait 780, 665, 540 do; 731 is matched by 1 cycle and 572 (0.353) by 275, so 572 lies below its own null q99 (0.633).
- **Null scale differs between runs.** Null medians at the truth positions are 0.24 to 0.45 (genotypic) and 0.46 to 0.60 (phenotypic); the pooled null quartiles over all detected positions are 0 / 0.197 / 0.374 and 0 / 0.295 / 0.541, and the observed quartiles 0.112 / 0.253 / 0.470 and 0.090 / 0.369 / 0.624. `CAAS_score` is therefore not comparable across the two runs; the p-values, calibrated within each run, are.
- The phenotypic null detects each of the five strongest truth positions less often (543 to 821 cycles vs 867 to 936), consistent with 3-of-3 all-or-nothing detection per hypothesis (§4).

## 2. Truth set, side by side

| pos | tier | geno score rank /59 | geno p.adj_bh | geno p.adj_bh_fact | pheno score rank /57 | pheno p.adj_bh | pheno p.adj_bh_fact |
|---|---|---|---|---|---|---|---|
| 780 | mutagenesis | 1 | **0.026** | **0.021** | 1 | **0.016** | **0.015** |
| 665 | mutagenesis | 1 | **0.026** | **0.013** | 2 | **0.016** | **0.016** |
| 540 | selection | 1 | **0.026** | **0.013** | 3 | **0.016** | 0.088 |
| 572 | selection | 20 | 0.75 | 1.00 | 4 | **0.016** | **0.045** |
| 733 | parallel | absent |  |  | absent |  |  |
| 761 | parallel | 18 | 0.58 | 0.58 | 19 | 0.62 | 0.34 |
| 749 | weak | 31 | 0.44 | 0.54 | 39 | 0.77 | 0.98 |
| 505 | weak | 8 | 0.44 | 0.24 | 8 | 0.17 | 0.16 |
| 573 | weak | 39 | 1.00 | 1.00 | 20 | 0.68 | 0.47 |
| 731 | weak | 6 | **0.042** | 0.053 | 6 | **0.016** | **0.031** |

| criterion | genotypic | phenotypic |
|---|---|---|
| truth positions with `p.adj_bh` below 0.05 | 4/10 (540, 665, 731, 780) | 5/10 (540, 572, 665, 731, 780) |
| truth positions with `p.adj_bh_fact` below 0.05 | 3/10 (540, 665, 780) | 4/10 (572, 665, 731, 780) |
| truth positions among the 10 highest `CAAS_score` | 5 | 6 |
| non-truth positions with `p.adj_bh` below 0.05 | 1 (751) | 1 (518) |
| non-truth positions with `p.adj_bh_fact` below 0.05 | 2 (611, 751) | 1 (518) |

- **Recovery under the non-circular trait.** Both mutagenesis sites and both selection-tier sites have `p.adj_bh` below 0.05 under the phenotypic trait, and 780, 665 and 572 also under `p.adj_bh_fact` (540: 0.088). They hold score ranks 1 to 4 (0.915 to 0.904); 518 (outside the truth set) and 731 follow at ranks 5 and 6, and the next position scores 0.812.
- **The genotypic hits are partly circular.** 780 defines the genotypic label; 540 and 665 carry the derived residue in 23/23 genotypic-C4 tips (`pepc_results.md`, caveats). All three share the maximum of the score scale (1.000) and rank 1; their genotypic p-values are close to definitional and are not evidence for the method. Under the genotypic trait 572 is not below 0.05 (`p.adj_bh` 0.75).
- **The weak tier does not separate from the background** in either run, except 731: `p.adj_bh` 0.042 genotypic and 0.016 phenotypic, `p.adj_bh_fact` 0.053 and 0.031. 505 under the phenotypic trait is the nearest of the others (0.17; 0.16).
- **Limit.** The top-10 cutoff is post hoc and descriptive. With 1000 cycles the floor of `p.emp` is 0.001; positions at the floor share one `p.adj_bh` value, so their order rests on `CAAS_score` and on the factorized p.

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
| `p.adj_bh < 0.05` | 751 | 518 |
| `0.05 ≤ p.adj_bh < 0.1` | 620 | 460, 474 |
| `p.adj_bh_fact < 0.05` | 611, 751 | 518 |
| `0.05 ≤ p.adj_bh_fact < 0.1` | none | 620 |
| FADE BF ≥ 100 as well (`p.adj_bh < 0.1`) | none | 474, 518 |

518 is below 0.05 under both adjustments only under the phenotypic trait (0.0157 and 0.0154; `p.emp` at its floor); under the genotypic trait it ranks 29 of 59 (`p.adj_bh` 0.61). The difference comes from pruning the five L-carrying *Eleocharis* accessions (`pepc_results.md`, Method 2). In the genotypic run 751 is below 0.05 under both adjustments (0.026 and 0.013) and has no FADE call; 611, a FADE call in both runs (BF 7.6 × 10³ and 1.2 × 10³), is below 0.05 only under `p.adj_bh_fact` (0.013; `p.adj_bh` 0.11), and 620 has `p.adj_bh` 0.088 and `p.adj_bh_fact` 0.20. In the phenotypic run 620 has 0.17 and 0.100. 460 and 474 are detected in 2 and 1 hypotheses with low tip coverage; hypothesis recurrence and coverage do not enter `CAAS_score` or `p.emp`, so such positions can rank alongside positions detected in every hypothesis.

## 6. Multiple-testing family

### 6a. Why candidate-set p-values are sub-uniform

`p.emp = (k + 1)/(N + 1)` counts null cycles that **detect** the position **and** reach the observed score. Its implicit test statistic is T = `CAAS_score` if detected, −∞ otherwise. Over all positions this is a valid permutation p-value; a position never detected in the observed data has T = −∞ and p = 1. A detected position with a score of 0 takes p = 1 by rule.

The candidate set keeps only positions the observed data detected, i.e. it filters on T. For a truly null position that passes this filter and has a positive score, `p.emp` ≈ d × U, where d is the fraction of null cycles that detect the position, so `p.emp` cannot exceed ≈ d. Across candidate positions d has median 0.53 (genotypic, range 0.10 to 0.96) and 0.45 (phenotypic, 0.07 to 0.84). This is the non-independent filtering case of Bourgon et al. (2010): any FDR procedure applied only to the survivors is anti-conservative.

### 6b. Corrections compared

Counts of positions below 0.1, except the raw count (below 0.05). "Candidate set" = observed-detected positions (59 / 57). "Null universe" = candidate set ∪ every position detected in ≥ 1 null cycle (106 / 94), with p = 1 for positions the observed data did not detect.

| procedure | family | genotypic | phenotypic |
|---|---|---|---|
| raw `p.emp` < 0.05 | none | 10 | 14 |
| BH of `p.emp` | candidate set | 7 (min 0.015) | 8 (min 0.009) |
| Storey, λ = 0.5 | candidate set | 7 of 59 (π̂₀ = 0.54) | 20 of 57 (π̂₀ = 0.49) |
| **BH of `p.emp` (`p.adj_bh`, pipeline)** | **null universe** | **6** (min 0.026) | **8** (min 0.016) |
| **BH of `p.emp_fact` (`p.adj_bh_fact`, pipeline)** | **null universe** | **6** (min 0.013) | **7** (min 0.015) |

- **Storey on the candidate set is invalid, not merely unstable.** π̂₀ = #{p > λ}/(m(1 − λ)) assumes null p-values fill [0, 1]; §6a caps them near d, so π̂₀ is biased low (0.54 and 0.49) and Storey counts 7 and 20 positions where BH over the null universe counts 6 and 8.
- **Candidate-set BH** counts 7 and 8 against 6 and 8, with smaller adjusted values because m is 59 and 57 instead of 106 and 94.
- **On the null universe**, restoring the p = 1 mass makes π̂₀ ≈ 1, and Storey reduces to BH.

### 6c. Factorized p and calibration

`p.emp_fact = (nd + 1)/(N + 1) × (1 + #{detections of the class with score ≥ s})/(1 + #{detections of the class})`, where `nd` is the number of null cycles that score the position and the class is one of 20 percentile classes of the null detections (`FACT_PROP_CLASSES`). The evaluated position counts as one more detection of its own and is classed by `nd + 1`. It is not bounded below by 1/(N + 1). `p.adj_bh_fact` is its BH over the null universe.

Calibration: each of the 1000 null cycles in turn is treated as observed and scored against the remaining 999, with the `.fact_*` functions of `scoring_compute.R` (`pepc_null_calibration.R`). A position the cycle does not score has p = 1. Under exchangeability the share of (position, cycle) pairs with p ≤ α cannot exceed α. For `p.emp` the share divided by α is 0.84, 0.82 and 0.79 (genotypic) and 0.74, 0.81 and 0.77 (phenotypic) at α = 0.001, 0.01 and 0.05. The share of null cycles, taken as observed, with at least one position below 0.05 and 0.1 after BH of `p.emp` is 0.001 and 0.010 (genotypic) and 0.009 and 0.065 (phenotypic). The script writes the same tables for `p.emp_fact` and `p.adj_bh_fact` to the `calibration/` folder of the run.

- **`p.emp` is conservative** at every α, as expected for a statistic that counts ties and detection at the same time; BH of `p.emp` spends well under its budget, most visibly in the genotypic run where the 1/N floor binds.
- **Scope.** Null cycles are draws from one generator, so this checks the estimators under exchangeability, not whether the observed labelling is exchangeable with the null ones (§6d). Pairs of the same position or the same cycle are not independent, and no standard error is given for the shares. Behavior under partial nulls and at genome scale is not tested.

### 6d. Is the observed labelling exchangeable with the null?

| check | genotypic | phenotypic |
|---|---|---|
| observed detections vs null per-cycle detections | 59; null median 44, q95 58; P(null ≥ obs) = 0.049 | 57; null median 32, q95 57; P = 0.056 |
| observed vs pooled-null `CAAS_score` quartiles | 0.112 / 0.253 / 0.470 vs 0 / 0.197 / 0.374 | 0.090 / 0.369 / 0.624 vs 0 / 0.295 / 0.541 |
| null cycles with the observed design (100 hypotheses) | 1000/1000 | 1000/1000 |

The observed data detect as many positions as the upper tail of the null cycles, and score above their median; this is what signal would produce and is not by itself diagnostic. Because every null cycle carries the same 100-hypothesis design as the observed data, the null is conditioned on design size by construction. The real labelling, routed through the null path as cycle b_0, reproduces every position score of the observed run (maximum difference 1 × 10⁻¹⁶ over 67 and 71 position-side scores), so the scoring is shared between the two paths. Whether the full observed path is exchangeable with the null under no association is tested by the negative controls (§6e).

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

On the same sedge PEPC data, Morel et al. (2024, Table 2) report PCOC recovering 7 of 12 convergent mutations under the genotypic annotation and none of 11 under the phenotypic one, while ConDor retains the best phenotypic precision (0.57). PhyloPhere CAAS gives both mutagenesis sites and both selection-tier sites an adjusted p below 0.05 under the phenotypic annotation. Morel et al. report FADE, with all branches of convergent clades as foreground, as the best method on both annotations (F₁ 0.81 genotypic, 0.57 phenotypic); PhyloPhere's FADE, run with the same foreground definition, recovers 8/10 and 7/10 of this truth set. The comparison is qualitative only: Morel et al. use a 458-column alignment with their own convergent-mutation set (12 genotypic / 11 phenotypic mutations; e.g. 749 recorded as M→T, where this truth set has L→T), so counts are not transferable.

## 9. Limits

- One run per trait, one seed (1998); `p.emp` resolution ≈ 0.001, and ranks among the top positions turn on a few null cycles.
- The truth set is not phenotype-independent: selection- and weak-tier sites come from Besnard et al. (2009) tests on branches defined by Ser780. Only the mutagenesis tier (780, 665) and the parallel tier (733, 761) carry independent evidence, and 733 is uninformative in this sample.
- The phenotypic trait is not clean either: at gene level it mislabels non-C4 paralogs of C4 species. A gene-copy-aware phenotypic label would require copy assignment, which in this dataset is made by the Ser780 residue itself; the two concerns cannot be fully separated with this fixture.
- The mechanism proposed for the null-scale difference (§4) is a hypothesis.
- n = 10 truth positions in one gene.
- The negative controls (§6e) were not repeated with the factorized p.

## References

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. Proc Natl Acad Sci USA 100(16):9440–9445. doi:10.1073/pnas.1530509100.
