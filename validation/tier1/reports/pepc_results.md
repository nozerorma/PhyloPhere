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

Columns: `Gene, Position, n_schemes, scheme_set, n_hypotheses, participating_hypotheses, top_species_residues, bottom_species_residues, n_top_species, n_bottom_species, n_conserved_pairs, CAAS_score, side, caas, amino_encoded, ancestral_aa, derived_aa, p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact`. One row per (Position, side); a position detected on both sides has two rows sharing one `p.emp`, `p.adj_bh`, `p.emp_fact` and `p.adj_bh_fact`. `top_species_residues`, `bottom_species_residues`, `n_top_species`, `n_bottom_species` and `n_conserved_pairs` tally the species of the contrast pairs; the residue columns of the tables below count every tip with a trait value, read from the alignment (gaps omitted). Definitions are in `docs/scoring_v2_p_emp.md`.

- `p.emp`: `(k_emp + 1)/(N + 1)`, where `k_emp` counts null cycles that re-detect the position on any side with max-over-sides CAAS ≥ the observed max-over-sides `CAAS_score`; `N` = 1000. A position with a score of 0 has `p.emp` = 1.
- `p.adj_bh`: BH with one test per position over the null universe, i.e. every position detected in ≥ 1 null cycle plus every observed position, with p = 1 for null-detected positions the observed data did not detect.
- `p.emp_fact`: `(nd + 1)/(N + 1) × (1 + #{detections of the class with score ≥ s})/(1 + #{detections of the class})`, where `nd` is the number of null cycles that score the position and the class is one of 20 percentile classes of the null detections. It is not bounded below by `1/(N + 1)`.
- `p.adj_bh_fact`: BH of `p.emp_fact` over the same null universe as `p.adj_bh`.
- Null design: every null cycle is matched to the observed design, 100 FOP hypotheses of K pairs each (K = 4 genotypic, 3 phenotypic); `caas_permulation/resample_perms.tab` holds 100 hypothesis labellings for each of the 1000 cycles.
- The four p-values (`p.emp`, `p.adj_bh`, `p.emp_fact`, `p.adj_bh_fact`) are informative: the pipeline derives no threshold or significance column from them. This report marks values below 0.05 in bold and counts positions below 0.05 and 0.1 only to describe them; `p.emp` at its floor `1/(N + 1)` is written `<0.001`. `CAAS_score` is `US + mean(GS)`: the exact-identity score (US) plus the mean of the grouping-scheme scores (GS1 to GS4) of the schemes that scored the position, each term 0 when absent; it lies between 0 and 2.

Note for ad hoc pandas reads: the `caas` value `N/A` (maize 573) is parsed as missing under pandas defaults; read with `keep_default_na=False`.

| | genotypic | phenotypic |
|---|---|---|
| position-rows | 67 | 71 |
| unique positions (candidate set) | 59 | 57 |
| BH family (null universe) | 106 (47 null-only) | 95 (38 null-only) |
| null cycles `N` | 1000 | 1000 |
| minimum attainable `p.emp` | 0.000999 | 0.000999 |
| positions with `p.emp < 0.05` / `p.emp_fact < 0.05` | 10 / 18 | 14 / 22 |
| positions with `p.adj_bh < 0.05` / `< 0.1` | 5 / 6 | 6 / 8 |
| positions with `p.adj_bh_fact < 0.05` / `< 0.1` | 5 / 6 | 5 / 7 |
| minimum `p.adj_bh` / `p.adj_bh_fact` | 0.0265 / 0.0132 | 0.0157 / 0.0154 |

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

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | score rank /59 | n_hyp /100 | p.emp | p.adj_bh | BH rank | p.emp_fact | p.adj_bh_fact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:22 \| A:53 | 2.000 | 1 | 49 | <0.001 | **0.026** | 1 (4) | **0.0010** | **0.021** |
| 665 | H→N | mutagenesis | N:23 \| H:52,N:2 | 2.000 | 1 | 100 | <0.001 | **0.026** | 1 (4) | **0.0005** | **0.013** |
| 540 | P→T | selection | T:23 \| P:53,S:1 | 2.000 | 1 | 100 | <0.001 | **0.026** | 1 (4) | **0.0005** | **0.013** |
| 572 | E→Q | selection | Q:18,K:5 \| E:52,Q:2 | 0.707 | 20 | 100 | 0.276 | 0.75 | 39 (1) | 0.43 | 1.00 |
| 733 | F→V | parallel | not in candidate set | | | | | | | | |
| 761 | S→A | parallel | A:14,S:8 \| S:53 | 0.779 | 18 | 49 | 0.182 | 0.58 | 32 (3) | 0.15 | 0.58 |
| 749 | L→T | weak | L:9,M:9,T:4 \| L:52,P:1 | 0.502 | 31 | 34 | 0.059 | 0.44 | 8 (10) | 0.13 | 0.54 |
| 505 | F→L | weak | L:17,F:6 \| F:54 | 1.221 | 8 | 100 | 0.043 | 0.44 | 8 (10) | **0.037** | 0.24 |
| 573 | A→N | weak | N:15,A:8 \| A:52,G:2 | 0.400 | 39 | 100 | 0.429 | 1.00 | 40 (20) | 0.64 | 1.00 |
| 731 | I→V | weak | V:17,Y:5 \| I:51,V:3 | 1.407 | 6 | 76 | 0.0020 | **0.042** | 5 (1) | **0.0030** | 0.053 |

**Adjusted p below 0.05:** `p.adj_bh` 4/10 (540, 665, 731, 780), `p.adj_bh_fact` 3/10 (540, 665, 780; 731 at 0.053); **9/10 present, 1/10 absent** (733). 540, 665 and 780 share the maximum score of the scale (2.000) and rank 1; non-truth 611 and 751 rank 4 and 5 (1.802 and 1.679), 731 ranks 6, and 572 (rank 20, `p.adj_bh` 0.75) is not below 0.05.

### Phenotypic trait (`c4_phenotypic`)

| position | ref>alt | tier | residues C4 \| C3 | CAAS_score | score rank /57 | n_hyp /100 | p.emp | p.adj_bh | BH rank | p.emp_fact | p.adj_bh_fact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 780 | A→S | mutagenesis | S:16,A:2 \| A:49,S:1 | 1.830 | 1 | 73 | <0.001 | **0.016** | 1 (6) | **0.0003** | **0.015** |
| 665 | H→N | mutagenesis | N:17,H:3 \| H:47,N:3 | 1.824 | 2 | 92 | <0.001 | **0.016** | 1 (6) | **0.0005** | **0.016** |
| 540 | P→T | selection | T:17,P:3 \| P:48,S:1,T:1 | 1.812 | 3 | 97 | <0.001 | **0.016** | 1 (6) | **0.0056** | 0.088 |
| 572 | E→Q | selection | Q:17,E:3 \| E:47,Q:3 | 1.809 | 4 | 97 | <0.001 | **0.016** | 1 (6) | **0.0024** | **0.045** |
| 733 | F→V | parallel | not in candidate set | | | | | | | | |
| 761 | S→A | parallel | A:13,S:5 \| S:49,A:1 | 1.006 | 19 | 47 | 0.211 | 0.62 | 31 (2) | 0.098 | 0.34 |
| 749 | L→T | weak | M:9,L:6,T:3 \| L:48,P:1,T:1 | 0.376 | 39 | 29 | 0.297 | 0.77 | 36 (1) | 0.43 | 0.98 |
| 505 | F→L | weak | L:16,F:4 \| F:49,L:1 | 1.597 | 8 | 84 | 0.016 | 0.17 | 9 (1) | **0.022** | 0.16 |
| 573 | A→N | weak | N:14,A:6 \| A:49,N:1 | 0.994 | 20 | 58 | 0.254 | 0.68 | 33 (3) | 0.15 | 0.47 |
| 731 | I→V | weak | V:16,I:3 \| I:46,V:4 | 1.786 | 6 | 85 | <0.001 | **0.016** | 1 (6) | **0.0013** | **0.031** |

**Adjusted p below 0.05:** `p.adj_bh` 5/10 (540, 572, 665, 731, 780), `p.adj_bh_fact` 4/10 (572, 665, 731, 780; 540 at 0.088); **9/10 present, 1/10 absent** (733). The four mutagenesis- and selection-tier sites hold score ranks 1 to 4 (780, 665, 540, 572: 1.830 to 1.809); 518 (non-truth) and 731 follow at ranks 5 and 6 (1.800 and 1.786), and the next position, 611, scores 1.624.

### Position 733

Absent from the candidate set under both traits because the fixture carries no C4-specific residue there: C4 tips are `F:19, V:2, M:1, gap:1` (genotypic) and C3 tips `F:54`. The F→V change reported for grasses and sedges (Besnard et al. 2009, Table 2) is not a C4-group-wide state in this sequence sample.

### Non-truth positions with `p.adj_bh < 0.1` or `p.adj_bh_fact < 0.1`

| trait | position | residues C4 \| C3 | schemes | CAAS_score | n_hyp /100 | p.emp | p.adj_bh | p.adj_bh_fact |
|---|---|---|---|---|---|---|---|---|
| genotypic | 611 | L:22,F:1 \| F:44,L:10 | GS3+GS4+US | 1.802 | 87 | 0.0070 | 0.11 | **0.013** |
| genotypic | 620 | C:17,A:3,F:1,S:1,T:1 \| S:51,T:3 | GS1–GS4+US | 1.266 | 100 | 0.0050 | 0.088 | 0.20 |
| genotypic | 751 | F:11,Y:11 \| Y:53 | GS2+GS4+US | 1.679 | 49 | <0.001 | **0.026** | **0.013** |
| phenotypic | 460 | E:11,D:1 \| D:20,E:3,X:2 | US | 0.901 | 2 | 0.0070 | 0.082 | 0.16 |
| phenotypic | 474 | G:10,E:2 \| E:16,D:4,G:4,K:1 | GS1–GS4+US | 0.844 | 1 | 0.0040 | 0.054 | 0.16 |
| phenotypic | 518 | F:20 \| L:42,F:8 | GS3+GS4+US | 1.800 | 100 | <0.001 | **0.016** | **0.015** |
| phenotypic | 620 | C:13,S:4,A:2,T:1 \| S:46,T:3,A:1 | GS1–GS4+US | 1.254 | 86 | 0.018 | 0.17 | 0.100 |

The two adjustments rank positions differently. `p.adj_bh` follows each position's own null (`p.emp`); `p.adj_bh_fact` follows the observed score against the scores of the null detections of positions with a similar detection frequency. In the genotypic run 751 is below 0.05 under both (0.026 and 0.013), 611 only under `p.adj_bh_fact` (0.013; `p.adj_bh` 0.11), and 620 only near 0.09 under `p.adj_bh` (0.088; `p.adj_bh_fact` 0.20). In the phenotypic run 518 is below 0.05 under both (0.0157 and 0.0154); 474 and 460 are above 0.05 under `p.adj_bh` (0.054 and 0.082) and above 0.1 under `p.adj_bh_fact`; 620 is nearest under `p.adj_bh_fact` (0.100). A position below 0.05 under only one adjustment is the least secure call.

- **518** is present in all 100 phenotypic hypotheses; every phenotypic-C4 tip carries F. Under the genotypic trait the C4 group is `F:18, L:5` and the position is weak (`CAAS_score` 0.525, 14 hypotheses, `p.adj_bh` 0.61). The five L tips are the Ser780-carrying *E. baldwinii* (FM208014/015/016) and *E. vivipara* (AB085948, FM208029) accessions, which the phenotypic run prunes as C3/C4 intermediates. F therefore marks the C4 lineages other than *Eleocharis*; it is not associated with the C4-type ppc-1 copy in *Eleocharis*. Within C4 *Fimbristylis* both ppc-1 copies carry F (the non-C4 paralogs FM208032/034/036 included), so at 518 F is a lineage state rather than a property of the C4-recruited copy. CT_DISAMBIGUATION's ancestral reconstruction (`ct_disambiguation/caas_convergence_master.csv`, US scheme, top side) places L at the MRCA of the three C4 domains carrying F in 100, 93 and 74 of the 100 hypotheses (posterior 0.94 to 1.0; in the other hypotheses the third domain's MRCA is F, with posterior 0.56 to 0.57), and FADE maps three L→F substitutions on foreground branches, so the association reflects repeated gains of F rather than its retention. F also occurs in the outgroup and in early-diverging C3 lineages (*Coleochloa*, *Microdracoides*, *Carpha*), consistent with an early shift to L in the in-group followed by returns to F in the C4 lineages.
- **460** and **474** are detected in 2 and 1 of 100 hypotheses and have low coverage (460: 12 of 20 C4 and 25 of 50 C3 tips carry a residue). Hypothesis recurrence is a descriptor in `scoring_compute.R` (§2g) and does not enter `CAAS_score` or `p.emp`, so narrow detection is not penalised. Their null detection rates are low, which is what lets a call in one or two hypotheses reach a `p.emp` of 0.004 to 0.007.

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
| genotypic | 611 | L | 7 640 | L:22,F:1 \| F:44,L:10 | F→L ×5 | 0.11 / **0.013** |
| genotypic | 579 | T | 6 180 | T:19,A:2,E:2 \| A:53,T:1 | A→T ×4 | 1.00 / 0.71 |
| genotypic | 839 | K | 924 | G:13,K:8 \| G:53 | G→K ×3 | 0.50 / 0.24 |
| genotypic | 509 | D | 138 | D:12,E:11 \| E:51,D:2,G:1 | E→D ×4 | 0.55 / 1.00 |
| genotypic | 514 | C | 128 | V:21,C:2 \| V:53,I:1 | V→C ×2 | not detected |
| genotypic | 630 | K | 114 | Q:12,K:11 \| Q:43,K:9 | Q→K ×3 | 1.00 / 1.00 |
| genotypic | 471 | T | 112 | T:9,K:5 \| T:26,Q:2 (35/77 gapped) | K→T ×2, Q→T, T→K | not detected |
| genotypic | 517 | A | 104 | T:14,A:9 \| T:53,A:1 | T→A ×3 | 1.00 / 1.00 |
| phenotypic | 518 | F | 68 800 | F:20 \| L:42,F:8 | L→F ×3 | **0.016** / **0.015** |
| phenotypic | 620 | C | 1 210 | C:13,S:4,A:2,T:1 \| S:46,T:3,A:1 | S→C, A→C, others | 0.17 / 0.100 |
| phenotypic | 611 | L | 1 200 | L:17,F:3 \| F:41,L:9 | F→L ×5 | 0.26 / 0.16 |
| phenotypic | 474 | G | 429 | G:10,E:2 \| E:16,D:4,G:4,K:1 (33/70 gapped) | E→G ×2, G→E | 0.054 / 0.16 |
| phenotypic | 630 | K | 146 | K:11,Q:9 \| Q:39,K:9 | Q→K ×3 | 1.00 / 1.00 |
| phenotypic | 514 | C | 129 | V:18,C:2 \| V:49,I:1 | V→C ×2 | not detected |

- **579, 611, 620** combine several reconstructed substitutions with a residue shared by C4 lineages that are not each other's closest relatives (*Cyperus*, *Eleocharis*, *Rhynchospora*, *Fimbristylis*, *Bulbostylis*). 611's L also occurs in about ten C3 tips across unrelated clades, so it is the least C4-specific of the three.
- **630** is also one of the two sites FUBAR flags for positive selection (P = 0.941). Residues are close to evenly split in C4, and K is present in nine C3 tips.
- **471** and **474** are heavily gapped and rest on 2 to 4 reconstructed changes; **514** rests on two V→C changes in two *Cyperus* tips. These are the weakest non-truth calls.
- **880** is not significant under either trait (BF 17 and 18): its I is a clade state of C4 *Cyperus*, reconstructed as one change on the clade stem. CAAS discovery calls it (US, 28/100 phenotypic hypotheses) but `CT_FILTER` discards it as part of a CAAS cluster (`postproc/filter_selected/*.minlen3.maxcaas70.tsv`, flag `Discarded`).

### CAAS vs FADE overlap

| | genotypic | phenotypic |
|---|---|---|
| `p.adj_bh < 0.05` | 5 (540, 665, 731, 751, 780) | 6 (518, 540, 572, 665, 731, 780) |
| `p.adj_bh_fact < 0.05` | 5 (540, 611, 665, 751, 780) | 5 (518, 572, 665, 731, 780) |
| FADE BF ≥ 100 | 16 | 13 |
| `p.adj_bh < 0.05` and FADE | 4 (540, 665, 731, 780) | 6 (518, 540, 572, 665, 731, 780) |
| `p.adj_bh < 0.05` without FADE | 751 (FADE BF 96) | none |
| `p.adj_bh_fact < 0.05` and FADE | 4 (540, 611, 665, 780) | 5 (518, 572, 665, 731, 780) |
| `p.adj_bh_fact < 0.05` without FADE | 1 (751) | none |

Under the phenotypic trait every position with an adjusted p below 0.05 is also a FADE call (6/6 and 5/5); under the genotypic trait the exception is 751 under both adjustments. FADE calls more sites than CAAS in both runs; its truth-set calls beyond the positions below 0.05 under `p.adj_bh` are 505, 572, 573, 749 (genotypic) and 505, 573 (phenotypic): weak-tier positions, and 572 under the genotypic trait, that CAAS does not separate from its null. The two methods share the alignment and the foreground definition, so their agreement is not independent corroboration; it reflects agreement between two models of the same labelled data (pair contrasts under a design-matched permulation null vs a branch-level substitution-bias model).

## Cross-method summary

| | FUBAR | CAAS, genotypic | CAAS, phenotypic | FADE, genotypic | FADE, phenotypic |
|---|---|---|---|---|---|
| Truth positions with adjusted p below 0.05 (FADE: BF ≥ 100) | 0/10 | 4/10 (`p.adj_bh`); 3/10 (`p.adj_bh_fact`) | 5/10; 4/10 | 8/10 | 7/10 |
| Both mutagenesis sites (780, 665) | no | yes | yes | yes | yes |
| Non-truth sites (FUBAR: P > 0.9; FADE: BF ≥ 100) | 2 | 751 (`p.adj_bh`); 611, 751 (`p.adj_bh_fact`) | 518; 518 | 8 | 6 |
| Truth positions absent from output | 0/10 | 1/10 (733) | 1/10 (733) | | |

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
