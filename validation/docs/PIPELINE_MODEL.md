# PhyloPhere pipeline model — what the validation suite is actually testing

Built from reading the code (2026-08-30), not from assumptions. Cited by
`DESIGN.md`. If any of this drifts from the code, the code wins — fix this file.

> **NOTE**: sections 6.4 and 7 discuss "Implications for Tier 0". Tier 0 is
> demoted (D-DIR-01) — read those as background on the pipeline's
> operative-fg/min-contrasts behaviour, not as live Tier 0 design.

Sources read: `subworkflows/CT/local/modules/{pindex,caas_id,disco,hyper}.py`,
`subworkflows/CT/local/scripts/{permulations.R,lean_contrast_selector.R}`,
`subworkflows/CT/{ct_resample,ct_bootstrap,ct_check_min_contrasts,caas_permulation}.nf`,
`workflows/contrast_selection.nf`, `subworkflows/TRAIT_ANALYSIS/local/src/commons.R`,
`workflows/rerconverge.nf`, `subworkflows/RERCONVERGE/local/{build_rer_trait,continuous_rer}.R`,
`subworkflows/SELECTION/selection_{prep,utils}.nf`, `subworkflows/FADE/local/src/run_hyphy_fade_batch.sh`,
`main.nf` wiring, a real result at
`~/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer_final_noprune/neoplasia_prevalence/`.

---

## 1. The phenotype is one continuous (or binary) trait table

Everything starts from `--my_traits`: a table with a `species` column and one or
more value columns (`--traitname` picks the column; `--n_trait`/`--c_trait` are
optional count columns). There is **no** "hand me explicit fg/bg groups" entry
point in the GUI recipe (`--caas_config` exists but bypasses contrast selection).

The three downstream modules each turn that trait into a foreground **their own
way**. They do not share an fg/bg definition.

---

## 2. CAAS — contrast selection produces the groups

### 2.1 contrast selection (`CONTRAST_ALGORITHM`, `pair_sel.f` / lean port)

Input: the continuous trait + the species tree. Steps:

1. **Discretise** the trait into high / low extremes by `discrete_method`
   (`quartile` q25/q75, `quintile` q20/q80, `decile` q10/q90, `median_sd`
   med±sd, `parameterized` = `bottom_quantile`/`top_quantile`). For a **binary**
   trait the 1s land in "high", the 0s in "low".
   *(For count data — `n_trait`+`c_trait` present — it instead uses Jeffreys-CI
   non-overlap as the candidate rule.)*
2. **Candidate pairs** = every (high tip, low tip) combination.
3. **Seed** with the pair of smallest patristic distance (largest trait diff
   breaks ties).
4. **Greedily add** the candidate with the highest *modified Dunn index*
   (min inter-cluster distance ÷ own diameter) until no more can be added.
5. **Validate independence**: overall Dunn ≥ 1 ⇒ every pair is phylogenetically
   independent of every other. This is what stops two nested foreground origins
   from both being counted.

Output: `traitfile.tab` — `species \t 1|0 \t pair_id`, **N pairs**, one fg tip +
one bg tip per pair. Real example (neoplasia): 4 pairs, 8 species total. Also
`boot_traitfile.tab` (the annotated trait table, kept for permulations) and the
pruned tree.

**`CHECK_MIN_CONTRASTS`**: if fewer than `min_contrasts` (default 3) foreground
tips survive, it writes `low_contrasts.skip` and the **entire CT/signification/
disambiguation chain is silently skipped**, run exits 0. Tier 0 must count how
often this fires.

### 2.2 CAAS discovery (`disco.discovery` → `caas_id.fetch_caas`)

Per gene, per alignment column, per trait:

- Build `fg_string` / `bg_string` = concatenated residues of the ungapped fg
  tips / bg tips, **ordered by `pair_id`**.
- `iscaas(fg_string/bg_string, max_conserved)`:
  - `max_conserved = int(n_pairs · (1 − min_divergent_fraction))`.
    Default `min_divergent_fraction = 0.5`, 4 pairs ⇒ `max_conserved = 2`: up to
    2 pairs may share their residue and it is still a CAAS. This is the "declared
    conservation" tolerance.
  - CAAS iff `overlap ≤ max_conserved` **and** (≥2 non-overlapping residues on
    the fg side **or** ≥2 on the bg side).
  - `pattern`: 1 = fg uniform & bg uniform; 2 = fg uniform, bg varied; 3 = fg
    varied, bg uniform; 4 = both varied. `--patterns` default `1,2,3`.
- `pvalue` column = `calcpval_random` (`hyper.py`): a hypergeometric probability
  of the observed fg/bg symbol partition **given that column's own residue
  composition**. This is a per-site surprise measure, **not** a permulation p and
  **not** multiple-testing corrected.
- **CAAP mode** (`caap_mode=true`, default): the same test is also run on
  physicochemically-encoded strings — `caap_group` ∈ {US (unencoded / strict),
  GS1, GS2, GS3} — so one position yields up to 4 rows, one per encoding scheme.

Output: `discovery.tab` (columns: gene, mode, caap_group, trait, position, caas,
[amino_encoded], pvalue, pattern, ffgn, fbgn, gfg, gbg, mfg, mbg, ffg, fbg, ms,
[is_conserved_meta, conserved_pair]).

### 2.3 The permulation null (`ct_resample` → `permulations.R`)

This is the null the DESIGN gate cares about ("KS-uniform p.perm").

- Read `target_pairs = max(pair_id)` from `traitfile.tab` V3 — **the null is
  filtered to the exact observed pair count**, deliberately not a parameter.
- Generate a **pool** of `perm_pool_size` accepted null labelings. Each draw:
  - `--perm_strategy`:
    - `FGBG` — plain label shuffle, no phylogeny.
    - `BM` — RERconverge `simpermvec`: Brownian simulation on the real tree,
      then **rank-match** the simulated ordering onto the sorted observed values
      (marginal distribution preserved exactly; only the phylogenetic dependence
      structure is randomised).
    - `lambda` — same, but simulate on a tree rescaled by the ML Pagel's λ of the
      empirical trait (accepts far more designs than BM). λ̂ ≈ 0.74 in the
      committed `perms_lambda` work.
  - Run the **same** `evaluate_lean_contrast_selection` (Dunn-gated greedy
    selector) on the draw. Accept iff it yields `target_pairs` pairs:
    Tier 1 = all pairs Dunn ≥ 1; Tier 2 = exactly one pair below 1 (fallback).
  - Budget = `max_tries`, escalated ×1.5 up to twice, else the run **fails
    loudly** rather than ship a thin null.
- Output: `resample_*.tab` (`b_0` = the real labeling if `include_b0`; then
  `b_1..b_pool`), `permulation_manifest.tsv` (tier / Dunn per record).

### 2.4 How the pool becomes p-values

- **Classic bootstrap** (`ct_bootstrap` with `--discovery` + `-s resample`):
  for each real CAAS, re-test that position under the pool → `bootstrap.tab`,
  `pvalue_boot` ≈ fraction of pool labelings that also call a CAAS there.
- **Permulation-excess / FCS null** (`caas_permulation.nf`, gated on
  `caas_permulation_enrichment`): seeded random subset of `caas_full_perms` pool
  labelings → full-pool bootstrap with `--export_perm_discovery` → **replay
  through disambiguation once (load ASR once, replay N labelings)** →
  `caas_perms.rds` (genes×N null score matrices on the `asr_path_score` axes) +
  `perm_pos_pval.tsv` (leave-one-out position-level `null_pvalue_boot`). Feeds
  the FCS `p.perm` path exactly like RER.

---

## 3. CAAS disambiguation + ASR path score

`CT_DISAMBIGUATION` re-analyses every discovered CAAS with ancestral state
reconstruction (codeml, **model = LG** by default — `ct_disambig_asr_model`).
`convergence_mode` ∈ `focal_clade` | `mrca`. Produces `asr_path_score`
(MRCA→root hop-weighted signed isolation/conservation score; see
`docs/ASR_PATH_SCORE.md` and memory `asr_path_score_axis_redundancy`) plus the
convergence/divergence/parallel classification. ASR mode `compute` runs codeml;
`precomputed` needs a cache dir (won't exist for simulated data).

`CT_POSTPROC` then applies `filter_minlen` / `filter_maxcaas` (single) or a
sweep (`minlen_values` × `maxcaas_values`), and `gene_filter_mode`
(none|extreme|dubious|both) drops outlier genes.

---

## 4. RER — raw trait, its own null

- RER_MAIN is called with an **empty** traitfile channel ⇒ falls back to
  `params.my_traits` (the original headered table, **not** contrast selection's
  output). `build_rer_trait.R` auto-detects binary vs continuous from the value
  distribution (exactly 2 unique non-NA values ⇒ binary, recode min→0/max→1).
- Binary: foreground = species with value 1 (`foreground2Tree`). Continuous: the
  values are used directly.
- Statistic: correlation (Rho) of per-gene relative evolutionary rate with the
  trait across the tree.
- **Null**: `getPermsContinuous` — Brownian-motion null phenotypes on the
  master tree, `rer_perm_batches` × `rer_perms_per_batch` (default 10×100).
  `p.perm` = proportion of null |Rho| ≥ observed |Rho|; `p.perm.adj` =
  BH across genes. Raw null matrices saved for pathway-level permulation.

RER does **not** use contrast selection's pairs at all.

---

## 5. FADE — extreme species, no permulation

- `SELECTION_PREP` → `EXTRACT_EXTREME_SPECIES` derives top / bottom species
  lists from `trait_stats.csv` (quantile cut on the trait). FADE foreground =
  the extreme set (**top** = trait-accelerated direction).
- HyPhy FADE tests **directional selection / mutational bias** toward each amino
  acid along the foreground branches. Per-site **Bayes Factor**; `fade_bf_threshold`
  default 100. Model **LG**. `fade_mode` all | gene_set.
- No permulation null — the BF is the evidence measure. Runs `top` and `bottom`
  directions separately.

FADE does **not** use contrast selection's pairs either.

---

## 6. SCORING — where the rankings actually come from

`scoring_compute.R` (`--scoring`, requires `--ct_postproc` output) integrates
CT_POSTPROC + FADE + RER + CT_ACCUMULATION into the tables the pipeline's
headline claims rest on. **Turning it off means Tier 0 has no gene/position
ranking to score power against** — the raw module outputs are evidence columns,
not a prioritisation.

### 6.1 Position level (`position_scores.tsv`)

Per (Gene, Position, scheme) row — schemes are `US, GS4, GS3, GS2, GS1` (the
CAAP encodings), treated **symmetrically, no per-scheme weight**:

```
phen_score  = 1 - percent_rank(pvalue_boot)      # permutation confidence, 0..1
asr_score   = asr_path_score                     # unified ASR/convergence/parallel signal
caas_row    = phen_score * asr_score
```

Aggregate to (Gene, Position): `CAAS_score = mean(caas_row over detecting schemes)`
(mean, not max — schemes are correlated readings of the same position; how *many*
schemes fire is already evidence). `n_schemes` / `scheme_set` kept as descriptors.
FADE site-level BF joined per position. The hypergeometric `pvalue` is used only
as a **significance gate** (`gate_fdr` / `gate_pvalue`), never as a score input.

### 6.2 Gene level (`gene_scores.tsv`) — separate axes, not one composite

| axis | how |
|------|-----|
| `gene_caas_score` | **size-adjusted max** of `CAAS_score` over the gene's positions — `F(max)^n`, which removes the `Spearman(n_positions, raw max) ≈ +0.47` order-statistic bias (memory `gene_caas_score_size_adjusted`). Also `_top` / `_bottom` by `change_side`. |
| `accum_cct_p` | Cauchy combination (CCT/ACAT) of the per-scheme accumulation empirical p-values; BH FDR over genes with ≥1 observed CAAS. |
| `rer_min_pval`, `rer_rho`, `rer_acceleration` | from the RER summary (`p.perm` preferred), `rer_significant = p ≤ 0.05`. |
| `fade_max_bf_top/bottom`, `fade_significant_*` | BF ≥ 100 per direction. |

Joins are **`full_join`** (memory `bugfix_scoring_gene_scores_caas_truncation`):
each module has its own, larger gene universe. `gene_scores.tsv` is sorted by
`gene_caas_score` desc. Cross-module correlation is computed on `gene_caas_score`
only; RER/FADE/accum enter as their native significance, not as score axes.

### 6.3 Ranked outputs

- `gene_lists/slice_*.tsv` — 12 slices (top/bottom/global × 25/10/5/1%), each
  ranked on its own axis (`gene_caas_score`, `gene_caas_score_top/bottom`), used
  by STRING/DOMINO.
- `fcs_stats*.tsv` — per-module ranking files feeding FCS (`fastwilcoxGMT`
  rank-shift AUC = Wilcoxon; plus a label-permuted score-accumulation NES null;
  plus the `caas_perms.rds` genome-wide **permulation-corrected `p.perm`** for
  the CAAS module, mirroring RER's pathway permulation).
- `gate_fdr` / `gate_pvalue` — deterministic scheme-priority pick (US > GS4 > …)
  for the significance gate (memory `bugfix_scoring_scheme_priority_tie`).

### 6.4 Consequence for Tier 0

The power question — "did the pipeline rank the planted genes / positions
highly?" — is a question about `position_scores.tsv` / `gene_scores.tsv` /
`gene_lists/`, i.e. it **needs SCORING to run**. The null question ("is `p.perm`
uniform?") is answerable from the CAAS-permulation / RER outputs alone, but the
*FCS* `p.perm` null also needs SCORING. So SCORING almost certainly belongs in
the Tier 0 module set — see Q-T0-D.

---

## 7. Implications for Tier 0 (to resolve WITH Miguel, not unilaterally)

1. **The planted foreground must be expressible as a continuous (or binary)
   `--my_traits` column**, because contrast selection / RER / FADE all start
   there. The simulator currently plants a discrete tip set; it needs to also
   emit the trait table. Two archetypes to build (Miguel's framing):
   `echo`-like binary 1/0, and `bodysize`-like continuous.
2. **contrast selection does not drop clades** — it forms `N` Dunn-independent
   pairs, one representative fg tip + one bg partner per independent origin. If
   the planted origins are non-nested and well separated (as `sample_foreground`
   makes them) and the convergent substitutions are planted on clade **stems**
   (so any representative carries them), the operative fg/bg should match the
   planted design closely. Residual mismatch (an origin lost to `min_contrasts`
   or Dunn) should be **measured**, not assumed away — that is the open Q2.
3. **Three different foregrounds per replicate**: CAAS sees the contrast-selection
   pairs; RER sees value==1 / the continuous vector; FADE sees the quantile
   extremes. Tier 0 scoring must score each module against *its own* operative
   foreground, and the harness must record all three.
4. **`min_contrasts`** silently skips the whole CT chain — Tier 0 must count and
   report skip rate, not just crash-or-pass.
5. The permulation null is **already phylogenetically structured** (BM/lambda +
   Dunn re-selection). A flat label-shuffle expectation is the wrong reference;
   the KS-uniformity gate is against the pool-derived `p.perm`, which is the
   right test.
