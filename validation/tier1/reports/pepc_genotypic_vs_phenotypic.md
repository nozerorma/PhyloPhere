# Tier 1 PEPC: Genotypic vs Phenotypic Trait

Two PhyloPhere runs (CAAS → CT_DISAMBIGUATION → FADE → SCORING → ENRICHMENT) on the same 970-column ppc-1 alignment and tree, differing only in the trait. Per-run tables, FUBAR baseline and numbering verification are in `pepc_results.md`; runtime in `pepc_resources.md`.

| | genotypic | phenotypic |
|---|---|---|
| results dir | `output/pepc/results/c4_complete/` | `output/pepc/results/c4_phenotypic_complete/` |
| trait column (`input/pepc/my_traits.tsv`) | `c4` | `c4_phenotypic` |
| definition | ppc-1 copy carries Ser780 (Besnard et al. 2009) | species-level anatomical/physiological call (Bruhl & Wilson 2007), assigned to every ppc-1 copy of that species (as in Morel et al. 2024) |
| tips (in-group) | 77 (23 C4 / 54 C3) | 70 (20 C4 / 50 C3); 7 *E. baldwinii* / *E. vivipara* C3/C4-intermediate accessions pruned |
| candidate positions (unique) | 61 | 51 |
| hypotheses × pairs | 100 × 4 | 100 × 3 |
| null cycles | 999 | 1000 |

`params.json` differs only in run-scoped paths, `traitname`/`secondary_trait`, and the prune flags. Position numbers below are maize PEPC1 (`position_scores.tsv` `Position` + 1; offset verified against the raw alignment for both runs). Significance: `p.emp_adj < 0.1`.

## 1. Position 780 and its co-circular block

`k_emp` = null cycles that re-detect the position and reach ≥ the observed max-over-sides `CAAS_score`. Null quantiles are over the cycles that detect the position. Recomputed from `caas_permulation/perm_pos_cycle_caas.tsv.gz` with the `scoring_compute.R` §2f-ter definition; `k_emp` reproduces the reported `p.emp` exactly in all cases.

| pos | trait | observed CAAS | n_hyp /100 | null detects | null q50 | null q90 | null q99 | k_emp | p.emp | p.emp_adj |
|---|---|---|---|---|---|---|---|---|---|---|
| 780 | geno | 0.673 | 63 | 869/999 | 0.276 | 0.468 | 0.615 | 3 | 0.004 | **0.058** |
| 780 | pheno | 0.715 | 68 | 601/1000 | 0.463 | 0.663 | 0.745 | 20 | 0.021 | 0.191 |
| 665 | geno | 0.683 | 100 | 907/999 | 0.243 | 0.453 | 0.615 | 3 | 0.004 | **0.058** |
| 665 | pheno | 0.713 | 89 | 664/1000 | 0.477 | 0.672 | 0.746 | 27 | 0.028 | 0.191 |
| 540 | geno | 0.683 | 100 | 900/999 | 0.251 | 0.452 | 0.592 | 3 | 0.004 | **0.058** |
| 540 | pheno | 0.710 | 92 | 788/1000 | 0.415 | 0.635 | 0.731 | 17 | 0.018 | 0.191 |
| 572 | geno | 0.437 | 100 | 885/999 | 0.094 | 0.363 | 0.493 | 35 | 0.036 | 0.263 |
| 572 | pheno | 0.704 | 92 | 637/1000 | 0.474 | 0.662 | 0.740 | 27 | 0.028 | 0.191 |
| 731 | geno | 0.442 | 78 | 742/999 | 0.131 | 0.387 | 0.522 | 30 | 0.031 | 0.263 |
| 731 | pheno | 0.696 | 78 | 616/1000 | 0.478 | 0.662 | 0.740 | 35 | 0.036 | 0.216 |

Magnitude at 780: `k_emp` 3 → 20, `p.emp` × 5.2, `p.emp_adj` 0.058 → 0.191 (× 3.3), crossing the threshold. The difference in `k_emp` is well outside Monte Carlo noise at N ≈ 1000 (Poisson SD ≈ 1.7 vs ≈ 4.5). 665 and 540 move in lockstep (× 7.0 and × 4.5 in `p.emp`).

The drop is carried by the null, not the observed score. Observed `CAAS_score` at 780 *rises* (0.673 → 0.715) while the null q99 rises further (0.615 → 0.745). The same upward null shift appears at every position inspected, including non-signal ones (818: null median 0 → 0.173), and the run-wide median observed `CAAS_score` goes from 0.117 to 0.355. `CAAS_score` values are therefore not comparable between the two runs; only `p.emp`, calibrated within each run, is.

## 2. Truth set, side by side

Rank = competition rank of `p.emp_adj` among unique candidate positions (tie-group size).

| pos | tier | geno p.emp | geno p.emp_adj | geno rank /61 | pheno p.emp | pheno p.emp_adj | pheno rank /51 | pheno/geno p.emp |
|---|---|---|---|---|---|---|---|---|
| 780 | mutagenesis | 0.004 | **0.058** | 1 (3) | 0.021 | 0.191 | 2 (6) | 5.2 |
| 665 | mutagenesis | 0.004 | **0.058** | 1 (3) | 0.028 | 0.191 | 2 (6) | 7.0 |
| 540 | selection | 0.004 | **0.058** | 1 (3) | 0.018 | 0.191 | 2 (6) | 4.5 |
| 572 | selection | 0.036 | 0.263 | 5 (14) | 0.028 | 0.191 | 2 (6) | 0.78 |
| 733 | parallel | absent | | | absent | | | |
| 761 | parallel | 0.169 | 0.329 | 25 (10) | 0.151 | 0.321 | 26 (1) | 0.89 |
| 749 | weak | 0.146 | 0.329 | 25 (10) | 0.266 | 0.379 | 36 (4) | 1.8 |
| 505 | weak | 0.170 | 0.329 | 25 (10) | 0.146 | 0.321 | 14 (12) | 0.86 |
| 573 | weak | 0.063 | 0.263 | 5 (14) | 0.202 | 0.349 | 27 (8) | 3.2 |
| 731 | weak | 0.031 | 0.263 | 5 (14) | 0.036 | 0.216 | 8 (1) | 1.2 |

733 is absent from both candidate sets because the fixture has no C4-group-wide residue there (C4 `F:19, V:2, M:1, gap:1`; C3 `F:54`).

## 3. How much recovery survives

| criterion | genotypic | phenotypic |
|---|---|---|
| truth positions with `p.emp_adj < 0.1` | 3/10 (780, 665, 540) | 0/10 |
| truth positions among the 10 lowest `p.emp` | 5 (540, 572, 665, 731, 780) | 5 (same five) |
| hypergeometric P(≥ 5 of 10), 9 truth positions in candidate set | 0.004 (of 61) | 0.009 (of 51) |
| Spearman ρ of `p.emp`, 44 shared candidate positions | 0.54 | |

- **Calibrated significance does not survive.** All three genotypic hits are the block that co-varies with Ser780 across ppc-1 copies (780 itself; 540 and 665 carry the derived residue in 23/23 genotypic-C4 tips). No truth position, and no non-780-block position, is significant under the genotypic trait either, so the genotypic run's significant recovery is entirely the circular block.
- **Relative ordering survives.** The same five truth positions occupy the 10 lowest `p.emp` values under both traits. Two non-circular-block sites improve under the phenotypic trait (572: rank 5 → 2, `p.emp` 0.036 → 0.028; 505: rank 25 → 14). 749 and 573 degrade.
- **Limit.** The top-10 cutoff is post hoc; the hypergeometric P is descriptive, not a pre-registered test. Ranks rest on `p.emp` values separated by a few null cycles.

## 4. Why significance is lost

### 4a. The four discordant tips are gene-copy mislabels

Residues at the truth positions, re-read from the alignment:

| tip | `c4` | `c4_phenotypic` | 780 | 665 | 540 | 572 | 761 | 749 | 505 | 573 | 731 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| *F. dichotoma* FM208033 | 1 | 1 | – | N | T | Q | – | – | L | N | – |
| *F. dichotoma* **FM208032** | 0 | **1** | – | H | P | E | – | – | F | A | I |
| *F. ferruginea* FM208035 | 1 | 1 | S | N | T | Q | A | T | L | N | V |
| *F. ferruginea* **FM208034** | 0 | **1** | A | H | P | E | S | L | F | A | I |
| *F. littoralis* FM208037 | 1 | 1 | S | N | T | Q | A | T | L | N | V |
| *F. littoralis* **FM208036** | 0 | **1** | A | H | P | E | S | L | F | A | I |
| ***Cyperus distichus*** | 1 | **0** | S | N | T | Q | A | T | L | N | V |

Each C4 *Fimbristylis* species contributes two accessions: one with the C4-type residue at every truth site, one with the ancestral residue at every truth site. Besnard et al. (2009, Results) report "two distinct ppc-1 gene clusters in C4 Fimbristylis", i.e. a duplication predating the clade, and that "at least one C4 ppc-1 sequence (i.e., encoding PEPC with Ser780) was isolated from all C4 species analyzed". On this evidence the discordant FM208032/034/036 are the non-C4 paralog of a C4 species (inferred from residues plus Besnard et al.'s duplication statement; copy identity was not checked against their gene tree). The phenotypic annotation assigns the species phenotype to every gene (Morel et al. 2024: "we annotated each gene using the annotation of the plant species in which it was sequenced"), so these paralogs are labelled C4. *Cyperus distichus* (*Volkiella disticha* in Bruhl & Wilson 2007) is the reverse case: a full C4-type sequence in a species classified C3. Its cause is unresolved (misclassification, lookup, or an unexpressed C4-type copy are all possible; none was tested).

At 780 these four tips account for every non-concordant residue: the only phenotypic-C3 tip with S is *C. distichus*, and the phenotypic-C4 tips with A are FM208034/036 (FM208032/033 are gapped). At 665 they account for all C4-side exceptions (H in FM208032/034/036); the C3 side also carries N in *Carpha glomerata* and *Hellmuthia membranacea* FM208041, which is the same under both traits.

### 4b. Hypotheses that lose detection contain those tips

| pos | phenotypic hypotheses not detecting | of which contain a discordant tip | remainder |
|---|---|---|---|
| 780 | 32 | 25 | 7, all using *F. dichotoma* FM208033, gapped at 780 (`max_*gaps_fraction = 0`) |
| 665 | 11 | 11 | 0 |
| 540 | 8 | 8 | 0 |
| 572 | 8 | 8 | 0 |

61/100 phenotypic hypotheses contain ≥ 1 discordant tip. With 3 pairs per hypothesis and `min_contrasts = 3`, one non-divergent pair removes the position from that hypothesis. Under the genotypic trait, each of the 37 hypotheses not detecting 780 contains exactly one pair with a gapped member at 780, and none of the 63 detecting hypotheses contains any: non-detection there is entirely gap-driven, as expected when label and residue coincide by construction.

### 4c. Contrast-pair structure changes with the trait

| | genotypic | phenotypic |
|---|---|---|
| pairs per hypothesis | 4 | 3 |
| distinct pairs | 135 | 106 |
| same-species pair-instances (C4 vs non-C4 ppc-1 copy of one species) | 41/400, in 35 hypotheses | 0/300 |
| cross-genus pair-instances | 67/400 (17 %) | 147/300 (49 %) |

Under the genotypic trait, 8 distinct within-species paralog pairs (*E. baldwinii*, *E. vivipara*, *F. dichotoma*, *F. ferruginea*, *F. littoralis*) are available and selected; each differs at every C4-associated site by the definition of the label. Under the phenotypic trait the *Eleocharis* pairs are pruned and the *Fimbristylis* pairs become C4–C4, so the selector falls back on more distant pairs (the three discordant *Fimbristylis* tips pair mostly with *Actinoschoenus thouarsii*). The loss of one pair per hypothesis follows from losing these within-species contrasts; the pair-count rule itself was not traced in `CONTRAST_ALGORITHM`.

### 4d. Null inflation: candidate causes, not isolated

The upward shift of the null `CAAS_score` distribution (§1) coincides with two structural changes: 3-pair all-or-nothing detection, and pairs spanning twice the proportion of genus boundaries. Either could raise the per-cycle score of whatever a random relabelling detects. Neither was tested in isolation (that would need, e.g., a phenotypic run with the 4 discordant tips relabelled, or pruned).

## 5. Novel (non-truth) positions

| | genotypic | phenotypic |
|---|---|---|
| significant (`p.emp_adj < 0.1`) | none | 859 |
| non-truth among 10 lowest `p.emp` | 584, 588, 620, 818 | 460, 625, 626, 751, 859 |
| overlap | none | |

The two traits share their truth-set signal but none of their non-truth candidates.

**859** (phenotypic only): residues C4 `K:13, G:4` vs C3 `K:46, R:4`; detected in 1/100 hypotheses (H51), on both sides. H51's pairs at 859 are K vs R, K vs R, and G (*F. littoralis* FM208037) vs R (*A. thouarsii*). G occurs only in the four *F. ferruginea* / *F. littoralis* accessions, i.e. a *Fimbristylis* clade state, and K is the majority residue in both groups. Its `CAAS_score` (0.803) is the highest in the run and `k_emp = 1`. Two properties of the scoring make this possible:

1. Hypothesis recurrence is descriptor-only (`scoring_compute.R` §2g: it "never multiplies CAAS_score"), so a 1/100 detection is not penalised relative to a 100/100 one.
2. `p.emp_adj` is BH over position-rows, and 859 has two rows (both sides) with the same pooled `p.emp`. BH over the 66 rows gives 0.0659; BH over the 51 unique positions gives 0.102, not significant. The three genotypic hits are unaffected by this choice (0.058 → 0.081).

Since `p.emp` is defined per (Gene, Position), BH over unique positions is the minimum consistent family; under it the phenotypic run has **0** significant positions of any kind (see §6 for the wider family).

**818** (genotypic, `p.emp_adj` 0.108, just above threshold): C4 `E:21` vs C3 `E:47, G:3, A:2, Q:1`, 2 hypotheses, observed `CAAS_score` 0.077. The null detects it in 275/999 cycles with q90 = 0, so a small positive observed score suffices for `p.emp = 0.009`. No convergent pattern is visible in the residues.

## 6. Multiple-testing family and correction

### 6a. Why the p-values are sub-uniform

`p.emp = (k + 1)/(N + 1)` counts null cycles that **detect** the position **and** reach the observed score. Its implicit test statistic is T = `CAAS_score` if detected, −∞ otherwise. Over *all* positions this is a valid permutation p-value; a position never detected in the observed data has T = −∞ and p = 1.

The candidate set, however, keeps only positions the observed data detected, i.e. it filters on T itself. For a truly null position that passes this filter, `p.emp` ≈ d × U, where d is the fraction of null cycles that detect the position. So `p.emp` cannot exceed ≈ d. Across the candidate positions d has median 0.49 (genotypic, range 0.12 to 0.94) and 0.47 (phenotypic, 0.09 to 0.85). The largest observed `p.emp` is 0.72 and 0.80, and no position sits near 1. This is the non-independent filtering case described by Bourgon et al. (2010): the filter is not independent of the p-value under the null, so any FDR procedure applied only to the survivors is anti-conservative.

### 6b. Corrections compared

Counts of positions below 0.1. "Candidate set" = observed-detected positions (61 / 51). "Null universe" = every position detected in ≥ 1 null cycle (107 / 95; the candidate set is a subset), with p = 1 for positions the observed data did not detect.

| procedure | family | genotypic | phenotypic |
|---|---|---|---|
| raw `p.emp` < 0.05 | none | 9 | 10 |
| raw `p.emp` < 0.1 | none | 18 | 16 |
| BH over position-rows (`p.emp_adj` as written in these runs' outputs) | candidate position-rows (72 / 66) | 3 (min 0.058) | 1 (859, 0.066) |
| BH | candidate unique positions (61 / 51) | 3 (min 0.081) | 0 (min 0.102) |
| Storey, λ = 0.5 | candidate unique positions | **61 of 61** (π̂₀ = 0.13) | **40 of 51** (π̂₀ = 0.24) |
| Storey, smoother / bootstrap π̂₀ | candidate unique positions | fails (π̂₀(λ) reaches 0 at λ ≥ 0.8) | fails |
| BH (current `scoring_compute.R` §2h; reproduced by re-running `SCORING_COMPUTE` on both task dirs) | null universe (107 / 95) | 0 (min 0.143) | 0 (min 0.190) |
| Storey, smoother | null universe | 0 (π̂₀ = 1) | 0 (π̂₀ = 1) |

### 6c. Reading

- **Rows vs unique positions.** `p.emp` is one test per (Gene, Position), but `scoring_compute.R` §2h runs BH over position-rows, so a position detected on both sides enters twice. Duplicating a p-value lowers both its own adjusted value (it gains a rank) and inflates *m*; the net effect here is anti-conservative for two-sided positions. It is the only thing making 859 significant.
- **Storey on the candidate set is invalid, not merely unstable.** π̂₀ = #{p > λ}/(m(1 − λ)) assumes null p-values fill [0, 1] uniformly. §6a caps them near d ≈ 0.5, so the region above λ is nearly empty and π̂₀ is driven towards 0. The resulting q-values would declare every genotypic candidate a discovery, including 818, which shows no convergent residue pattern. Estimating π₀ is precisely the step that relies most on the uniform-null assumption.
- **On the null universe, Storey reduces to BH.** Restoring the p = 1 mass gives π̂₀ = 1 (or 0.94), and q = BH.
- **Under the consistent family, nothing is significant under either trait.** The circular genotypic block (780, 665, 540) reaches BH 0.143, not below 0.1. The genotypic–phenotypic contrast in §1 (k_emp 3 vs 20 at 780) is unaffected: it is a within-position null comparison, independent of the correction.
### 6d. Permutation FDR (SAM-style) and calibration

Estimator: for a threshold t on the pooled max-over-sides score, FDR(t) = (mean number of positions per null cycle with score ≥ t) / (number of observed positions with score ≥ t), with π₀ = 1 and q made monotone. The expected null count is averaged over cycles, so it has no per-position 1/N floor.

Calibration: each of the 1000 null cycles in turn is scored as if it were the observed data, against the remaining 999. Under this complete null every call is false, so P(≥ 1 call) is the realised FDR.

| | genotypic | phenotypic |
|---|---|---|
| observed min q (top block) | 0.105 (780, 665, 540) | 0.333 (859, 780, 665, 540, 572, 731) |
| same, two random half-splits of the cycles | 0.095 / 0.113 | 0.308 / 0.318 |
| calibration, P(≥ 1 call) at q < 0.1 | 0.111 | 0.102 |
| calibration, P(≥ 1 call) at q < 0.2 / < 0.3 | 0.211 / 0.306 | 0.188 / 0.297 |
| same calibration for BH over the null universe, at 0.1 | 0.012 | 0.068 |

Scope of this check: null cycles are i.i.d. draws from one generator, so scoring one against the others tests the estimator's arithmetic under exchangeability, not whether the observed labelling is exchangeable with the null ones (§6e). Within that scope, permutation FDR is calibrated at nominal level under the complete null in both runs (Monte Carlo SE of the proportion ≈ 0.01), while BH over the null universe spends far less than its nominal budget, most visibly in the genotypic run, where the 1/N floor binds. Neither procedure calls anything in the observed data at 0.1; the genotypic top block sits on the threshold (half-splits 0.095 to 0.113), the phenotypic one well above it. Not tested: behaviour under partial nulls (π₀ < 1, where π₀ = 1 makes the estimate conservative), and at genome scale.

### 6e. Is `p.emp` itself calibrated?

Within the null, trivially yes: a leave-one-out rank among i.i.d. cycles is uniform. The substantive question is whether the observed labelling behaves like one more null draw when there is no association. Two checks are possible with the existing outputs; neither settles it.

| check | genotypic | phenotypic | reading |
|---|---|---|---|
| observed detections vs null per-cycle detections | 61; null median 44, q95 61; P(null ≥ obs) = 0.062 | 51; null median 34, q95 53; P = 0.075 | upper edge, inside the range; consistent with signal, not diagnostic |
| observed vs pooled-null `CAAS_score` quartiles | 0 / 0.117 / 0.290 vs 0 / 0.124 / 0.286 | 0.106 / 0.355 / 0.546 vs 0 / 0.284 / 0.524 | similar scale |
| null cycles with the full 100-hypothesis design (`resample_perms.tab`) | 795/1000 (min 1, q05 30) | 933/1000 (min 3, q05 80) | observed always has 100 |
| Spearman(hypotheses per cycle, detections per cycle) | 0.33 | 0.21 | fewer hypotheses, fewer detections |
| `p.emp` recomputed on full-design cycles only: share of positions that increase; median ratio | 0.90; 1.07 (top block 0.004 → 0.005) | 0.92; 1.02 | small anti-conservative bias |

Null cycles whose permuted labelling admits fewer than 100 valid FOP hypotheses detect fewer positions, while the observed design always has 100; comparing the observed design only with full-design cycles is the like-for-like conditional null. The bias this introduces here is a few percent of `p.emp`. Whether the observed path and the null path are otherwise exchangeable (contrast selection, discovery, disambiguation, postproc filters all mirrored) can only be tested with negative-control traits run through the observed path; that has not been done.

- **Limit.** The null universe is itself defined by this run's null (1000 cycles). Positions with very low detection rates may be missing from it, which would understate *m* slightly. Columns never detected under any labelling carry no permutation information and are excluded by construction.

## 7. Sanity checks on the phenotypic input

| check | result |
|---|---|
| trait counts (in-group) | 20 C4 / 50 C3, 7 pruned; matches `input/.pepc_phenotypic/README.md` |
| pruned *Eleocharis* tips in any contrast pair | 0 (genotypic run uses all 7) |
| pair members with missing trait | 0/600 in both runs |
| pair label structure | 300/300 C4–C3 (`abs_diff = 1`); genotypic 400/400 |
| discordant-tip usage (hypotheses) | *C. distichus* 25; FM208032 13; FM208034 19; FM208036 12 |
| DAG / params | only prune and trait parameters differ; 37/37 tasks succeeded |

Nothing pathological in pipeline mechanics. The discordant tips behave exactly as their residues predict, and they are the proximate cause of the lost significance (§4a, §4b).

The 2-tip count discrepancy (this fixture 20 C4 vs Morel et al.'s 22 on the same 71 tips) cannot be resolved here. One bound does follow from the alignment: *C. distichus* is the only phenotypic-C3 tip carrying Ser780, so any two further tips relabelled C4 would carry Ala or a gap at 780 and add label noise in the same direction as the *Fimbristylis* paralogs, not remove it.

## 8. Relation to Morel et al. (2024)

On the same sedge PEPC data, Morel et al. (2024, Table 2) report PCOC recovering 7 of 12 convergent mutations under the genotypic annotation and none of 11 under the phenotypic one, while ConDor retains the best phenotypic precision (0.57). The PhyloPhere result has the same shape: calibrated significance collapses, ranking of the known sites does not. The comparison is qualitative only. Morel et al. use a 458-column alignment with their own convergent-mutation set (12 genotypic / 11 phenotypic mutations; e.g. 749 recorded as M→T, where this truth set has L→T), so counts are not transferable.

## 9. Limits

- One run per trait, one seed (1998); `p.emp` resolution ≈ 0.001, and ranks among the top positions turn on a few null cycles.
- The truth set is not phenotype-independent: selection- and weak-tier sites come from Besnard et al. (2009) tests on branches defined by Ser780. Only the mutagenesis tier (780, 665) and the parallel tier (733, 761) carry independent evidence, and 733 is uninformative in this sample.
- The phenotypic trait is not "clean" either: at gene level it mislabels non-C4 paralogs of C4 species. A gene-copy-aware phenotypic label (species phenotype applied only to the C4-recruited copy) would require copy assignment, which in this dataset is made by the Ser780 residue itself; the two concerns cannot be fully separated with this fixture.
- The mechanism in §4d is a hypothesis.
- n = 10 truth positions in one gene.

## References

- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. Proc Natl Acad Sci USA 100(16):9440–9445. doi:10.1073/pnas.1530509100.
- Tusher VG, Tibshirani R, Chu G. 2001. Significance analysis of microarrays applied to the ionizing radiation response. Proc Natl Acad Sci USA 98(9):5116–5121. doi:10.1073/pnas.091062498.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.

- Besnard G, Muasya AM, Russier F, Roalson EH, Salamin N, Christin PA. 2009. Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence. Mol Biol Evol 26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4 photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148. doi:10.5642/aliso.20072301.11.
- Morel M, Zhukova A, Lemoine F, Gascuel O. 2024. Accurate detection of convergent mutations in large protein alignments with ConDor. Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
