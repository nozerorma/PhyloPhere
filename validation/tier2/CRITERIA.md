# Tier 2 (marine mammals): input, truth tiers and outcome criteria

DRAFT for review. Items marked **[decide]** need a decision before any run. The criteria are fixed
before the pipeline is run on this case, so that no outcome is judged after the fact.

## Input

- **Species:** the 58 placental mammals of the RER set (UCSC 100-way hg19 assembly codes); the four
  non-placentals (`ornAna1`, `monDom5`, `sarHar1`, `macEug2`) are removed.
- **Gene trees:** the 19,149 trees of `UCSC62mammaltrees.rds` (RERconverge object), pruned to the 58 species.
- **Species tree:** the master tree of that object ("Meredith+"). Pruned to the 58 placentals it is
  identical (Robinson-Foulds 0, both rooted at the armadillo) to the consensus topology in the
  supplement of Meyer et al. 2018, which is the Meredith et al. 2011 topology with three local
  resolutions (star-nosed mole outside hedgehog and shrew; cow outside Tibetan antelope, sheep and
  goat; ursids outside mustelids and pinnipeds). It differs from `hg19.100way.nh` by a Robinson-Foulds
  distance of 22 (maximum 110 on 58 tips). The master tree is used because the branch lengths of the
  gene trees were estimated on it.
- **Alignments:** UCSC 100-way amino acid alignments of the canonical transcripts, stitched per gene,
  restricted to the 58 species. Human (hg19) is the numbering frame.
- **Marine foreground:** `turTru2`, `orcOrc1` (Cetacea), `odoRosDi`, `lepWed1` (Pinnipedia), `triMan1`
  (Sirenia): three independent origins with 2+2+1 tips.
- **Pinniped ancestral branch:** both codings are run, marine and not marine. **[decide]** which one is
  the primary coding reported in the main text.
- **Negative contrast:** terrestrial control foregrounds (aardvark, alpaca, camel, microbat, David's
  myotis) run through the same pipeline; marine signal is read against them.

## Step 0: presence check (data only, no pipeline run)

For each truth item, check that the gene is in the alignment and that the five marine tips carry
the reported residue or an unambiguous alignment column. Items that fail are listed and excluded
before any result is seen.

### Step 0 results

- **RH1 (rhodopsin), transcript `uc003emt.3`, human numbering equal to the bovine numbering of the
  assays** (human D83, A292). Marine tips at the two sites:

  | Tip | Lineage | 83 | 292 |
  |---|---|---|---|
  | `turTru2` dolphin | Cetacea | N | S |
  | `orcOrc1` killer whale | Cetacea | N | S |
  | `odoRosDiv1` walrus | Pinnipedia | N | A |
  | `lepWed1` Weddell seal | Pinnipedia | D | A |
  | `triMan1` manatee | Sirenia | D | A |

  Site 292 carries the derived residue (S) only in the two cetaceans. Site 83 carries N in the two
  cetaceans and the walrus, not in the Weddell seal or the manatee, and N83 also occurs in terrestrial
  placentals of the alignment (panda, guinea pig, pika, elephant, three bats). The two sites are
  therefore not three-origin convergent sites on this alignment; as T1 items they test the Cetacea
  lineage only.
- **APPL1** (`uc003dio.3`, 709 aa) and **NEIL1** (`uc002bae.4`, 477-aa isoform of the alignment). The
  numbers reported by Xu 2025 (APPL1 378, NEIL1 71) do not hold the derived residue in the human
  protein of this alignment, and Yuan 2021 reports the same APPL1 P to L at 120, so the reported
  numbers belong to other frames. Scanning each whole protein for a column with the human residue
  and the derived residue in the marine tips gives one column per gene:

  | Gene | Human column | Marine tips with the derived residue | Terrestrial placentals with it |
  |---|---|---|---|
  | APPL1 | P421, derived L | dolphin, killer whale, Weddell seal, manatee (walrus keeps P) | 4 |
  | NEIL1 | E157, derived G | dolphin, killer whale, walrus, Weddell seal (manatee keeps E) | 1 |

  The match is by pattern within the gene, not by a verified position map. The NEIL1 offset (157 against
  71) is consistent with the 87 extra residues of this isoform relative to the 390-aa canonical
  protein; the APPL1 offset (421 against 378) has no such explanation.
- **Decision (T1 scoring):** items are scored per lineage carrying the residue. Two-level convergence
  is part of what the analysis validates, so a site shared by two of the three origins is a valid T1
  item for those lineages.

## Truth tiers

| Tier | Content | Items | Use |
|---|---|---|---|
| T1 experimental | Positions with a measured effect of the exact substitution | RH1 83 and 292; APPL1 P378L; NEIL1 E71G | Rank of each site among the tested positions of its gene |
| T2 external replication | Gene and derived residue from lists computed on other species sets and orthology | Zhou 2015 S1, Yuan 2021 S24, Xu 2025 S7, Huang 2024 S15, Foote 2015 Table 14, CAAP S12 | Share of residue-matched positions detected, against a null of random positions of the same genes |
| T3 gene sets | Independent gene-level classes | GO classes shared by RER and TRACCER (keratinization, cornification, protein-glutamine gamma-glutamyltransferase, myosin and muscle, olfactory receptors); pseudogenes of Meyer 2018; losses of Huelsmann 2019 | Enrichment of each pillar's ranking; no per-gene significance claim |

Reference null: the count of unique convergent substitutions of Thomas et al. 2017 on the same five
marine tips (131 with 15 species, 1 with 54 background species) is the expected order of magnitude of a
naive convergence call on this species set.

## Outcome criteria

- **CAAS:** one table per truth item with the native outputs of `position_scores.tsv`, as in the Tier 1
  report: presence in the candidate set, `CAAS_score` and its rank among all candidate positions, and
  `p.emp`, `p.adj_bh`, `p.emp_fact`, `p.adj_bh_fact`. The four p-values are informative and no
  threshold or significance column is derived from them; values below 0.05 are marked only to describe
  them. Items absent from the candidate set are reported as absent. A small total count is an
  admissible result: it is the expected behavior under dense background sampling and is not a failure.
  The case supports the tool if T1 and T2 items rank high among the genome-wide candidate positions and
  marine exceeds the terrestrial control foregrounds. **[decide]** the genome-wide statistic that
  summarizes the truth items together (for example the share of items inside the top fraction of
  candidate positions against the fraction expected by chance) and the margin that counts as support.
  With one BH family over all positions of 19,149 genes and `p.emp` bounded below by 1/(N+1), the
  adjusted values are expected to lose resolution, so rank is the primary readout.
- **RER:** gene-set level only. The T3 classes must be enriched in the accelerated tail, and the
  enrichment is compared with the published GO tables (Partha 2019 S7 to S12, TRACCER S2).
  Per-gene FDR is reported, not used as a criterion: the best published per-gene FDR on five tips is
  0.0455 (one gene).
- **Pillars are reported side by side and never combined into a single p-value.**
- **Sensitivity:** every result is reported under both pinniped codings.

## Limits stated in advance

- Three origins and five marine tips bound the power of any contrast-based method on this input.
- T2 lists come from other species sets; matches are by gene and derived residue after remapping to
  the human frame, not by raw position number.
- The alignment is a 2014 projection of whole-genome alignments onto human exons; alignment errors can
  create spurious convergent substitutions, which the T1 tier and the terrestrial controls address
  only in part.

## Toy run

A 121-gene toy (81 genes of interest of the truth tiers and 40 random background genes, `input/marine/README.md`) runs the same readouts with 100 null cycles to check the pipeline before the full run. Its results are in `reports/marine_toy_results.md`. At 100 cycles the smallest `p.emp` is 0.0099 and no adjusted p-value is informative, so only the rank-based readouts, the marine against control comparison and the RER gene-level signal are read; the toy does not decide any **[decide]** item.
