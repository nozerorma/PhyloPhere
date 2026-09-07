# `derived_agreement` harvest-wide: modal vs fractional — verdict

**Verification run:** `2.Primates .../cancer_no_prune_multi05_toy/malignant_prevalence_toy_complete`
(K=3 Voronoi domains, 4782 `filtered_discovery` rows → 792 `(Gene,Position,caap_group)`
rows → 267 scored positions, 32 genes; `tau=0.8`, `diversity_floor=0.75`).

Script: [`da_frac_verify.R`](da_frac_verify.R). Machine output: [`out/REPORT.md`](out/REPORT.md),
`out/rows_caap_group.csv`, `out/positions.csv`, `out/gene_caas.csv`,
`out/convergence_schemes_cmp.csv`, `out/top10_domain_breakdown.txt`.

## Reconstruction is exact

The modal pipeline (`fop_pool.R::apply_fop_pooling` → §2c–2g of `scoring_compute.R`)
was reproduced **bit-for-bit** from the flat TSVs: max |Δ| vs `position_scores.tsv`
is `0` for `derived_agreement`, `asr_score` and `CAAS_score`; the row-level asr
algebra matches the stored value to `3e-17`. So every number below is a like-for-like
swap of one factor.

### Caveat found along the way — FIXED

`position_scores.tsv` only reproduces with **equal-weight** pooling
(`hyp_pairs_path = NULL`). Root cause: `CHECK_MIN_CONTRASTS`
(`subworkflows/CT/ct_check_min_contrasts.nf`) rebuilt `traitfiles_ok_dir` with
only the `traitfile_H*.tab` files, dropping `contrast_hypotheses_pairs.tsv` —
and that dir is the source `main.nf` reads FOP PSS weights from. So **every FOP
run silently scored equal-weight**. Fixed: CHECK_MIN_CONTRASTS now carries the
PSS file through; `scoring.nf` gained a fallback chain
(`--scoring_hypotheses_pairs` → outdir auto-discovery); `conf/scoring.config`
declares the param. Verified: PSS weighting now moves `da` on 19/792 rows
(max |Δ| 0.083). The modal-vs-fractional numbers below use equal weights as the
common baseline (that is what the pre-fix run scored); the PSS-weighted
`da_frac` is reported alongside.

## What the fractional rule changes

| level | modal → fractional |
|---|---|
| `(Gene,Position,caap_group)` rows with \|Δda\|>1e-9 | **24 / 792 (3.0%)** |
| scored positions with da moved | **10 / 267 (3.7%)** |
| positions with `CAAS_score` moved | 10 / 267; max \|Δ\| **0.082**, mean over moved **0.030** |
| top-1% / top-5% position `CAAS_score` (Jaccard) | **1.00 / 1.00** (unchanged) |
| top-10% position `CAAS_score` (Jaccard) | 0.93 — one swap (`PRIMPOL 502` out, `CARD14 651` in) |
| `gene_caas_score` Spearman(modal, frac) | **0.9985** (n=32) |
| top-1/5/10% genes (Jaccard) | **1.00** — no gene enters or leaves any tier |
| `convergence_schemes` positions changed | **2 / 267** (`ZNF518B 705`: `GS1,GS3`→`""`; `ZNF518B 1018`: `GS4`→`""`) |

Back-compatibility holds exactly: of the **754** rows with no split domain, **0**
move (`|Δda| ≤ 1e-9`). Every one of the 24 differing rows has ≥1 domain whose
contrast pairs reconstructed to different residues across hypotheses.

## Why they differ — the two archetypes (from `top10_domain_breakdown.txt`)

1. **The modal fabricates agreement (22 of 24 rows, da falls, usually 1.0→0.75).**
   `PRIMPOL/502`: domain 1 is unanimously `D` across all 6 hypotheses; domain 2
   is `S,D,S,S,D,D` — a **3/3 tie**. The modal's `names(sort(table))[1]` picks `D`
   (alphabetical), so both domains "agree" → da = 1.0. The fractional rule sees
   domain 2 as `{D:0.5, S:0.5}`, giving `(1 + 0.5)/2 = 0.75`. This is spec-3.4's
   *genuine ambiguity* case: the arbitrary alphabetical tie-break invented a
   unanimous convergence signal out of a coin-flip domain.

2. **The modal's tie-break lands on the *wrong* residue (2 rows, da rises).**
   `ANK3/3881` (GS2): domain 1 → `x`; domain 3 is `D,Y,D,Y,Y,D` — a tie the modal
   breaks to `a` (`D`), which *disagrees* with domain 1 → da = 0.5. The fractional
   rule credits domain 3's half-weight on `x` → `(1 + 0.5)/2 = 0.75`. Spec-3.4's
   *partial alignment* case.

`convergence_schemes` inherits the same defect: `ZNF518B/705` domain 3
(`T,M,M,T`) is a tie the modal resolves to `M`, which co-encodes with domain 2's
`L` under GS1 → a spurious "chemical convergence" call that the fractional rule
(0.75 < 0.8) correctly withholds.

**18 of 792 rows carry an exact alphabetical tie-break that is load-bearing** —
it moves `da`, and in 2 cases a published `convergence_schemes` label.

PSS weighting (spec 3.1), had it been wired, would further move `da` on 19 rows
(max extra shift ~0.09, both directions) — it is *not* a no-op, though it changes
no ranking here.

## Recommendation

**Adopt the fractional rule. Do not fall back to `wmean(df$derived_agreement, row_w)`.**

- The `wmean` fallback is the *worst* of the three: it averages the
  already-within-hypothesis `da` values (each ≈1.0 here, since a domain has one
  pair per hypothesis) and is structurally blind to the between-hypothesis
  residue split — exactly the blindness POINT 3 was added to remove. Reverting to
  it would push all 24 of these rows back to ≈1.0, *more* inflated than the modal.

- The modal is a real defect, not just a stylistic one: its alphabetical
  tie-break is arbitrary, load-bearing in ~2.3% of scored rows, produces `da`
  values and `convergence_schemes` labels that the run currently publishes
  and cannot defend, and it is the only aggregate in `fop_pool.R` that ignores
  PSS. It also contradicts `residue_descriptors.py`, which keeps the full
  residue set for the same pair.

- The fractional rule is the correct primitive (the modal is its degenerate
  case), is exactly back-compatible on unanimous domains (0/754 violations here),
  is cheap, and is PSS-weightable like Job A/B.

- **Expectation management:** on data like this run the downstream effect is
  small — no change to top-5% position or gene rankings, `CAAS_score` moves by
  ≤0.08 on 3.7% of positions, 2 `convergence_schemes` flips. The value is
  *correctness and defensibility* of `da` / `convergence_schemes` where the modal
  is currently indefensible, not a re-ordering of hits. The split-domain case
  that motivates it appears in ~5% of scored rows, enough to be worth fixing
  properly rather than papering over.

- **Separately:** fix the PSS-file wiring into SCORING (background task
  `task_a5ce83ab`). It changed nothing decision-relevant in this toy run, but
  Job A/B and `da_frac` both depend on it and that is not guaranteed at genome
  scale.

### Implementation — SHIPPED 2026-09-07

- `fop_pool.R`: `.collect_changed_pair_dists(df, node_cols, top_cols, bot_cols, pair_pss)`
  → named list `"i|side" → residue→weight` (weight `pair_pss(hyp_id, i)`, else 1),
  each normalised to 1. `rebuild_derived_agreement_frac(dists, scheme)` sums
  `p_i^enc(g)` across a side's domains, `max_g / |D_s|`, mean over qualifying
  sides. `.make_pair_pss(hyp_pairs)` is the shared Job-A instrument.
  `pool_group`'s POINT 3 block and `.position_descriptors` (now taking
  `pair_pss`, passed through from `apply_fop_pooling`) both use them. The MODAL
  `.collect_changed_pairs` + `rebuild_derived_agreement` are retained for the
  `convergence_schemes` `>= 2 changed pairs` gate and the identity ("US") branch.
  `.wmean` fallback unchanged for `*_aa`-less inputs.
- `fop_pool.py` (null mirror): `_domain_side_dists` + `_da_frac_from_dists`,
  same math, PSS from `pss_by_hyp_domain`. Replaces the old
  `(domain,side,raw_aa)` dedup — R and py now agree bit-for-bit on this axis.
- VEP: `map_to_{cosmic,primateai}.py` take an optional `position_scores.tsv`
  arg; positions with `convergence_schemes == ""` are skipped;
  `anc_der_from_descriptor` narrows `der_aas -= anc_aas` (keep-all safety).
  Wired via `params.vep_position_scores` (standalone) / `SCORING.out.position_scores`
  (`--scoring`, VEP moved after SCORING; `VEP_STANDALONE` alias for the no-scoring
  path).
- Tests: `test_fop_pool.R` (+4 fixtures), `test_fop_pool.py` (+3),
  new `subworkflows/VEP/local/src/test_map_gate.py` (8). All green.
- Verified: `docs/da_frac/da_frac_verify.R` re-run — retained modal path
  byte-identical to the run's stored `derived_agreement` (max |Δ| 0);
  new `apply_fop_pooling` output == the script's independent reconstruction
  (max |Δ| 0); 24/792 rows move, 2 `convergence_schemes` flips, as predicted.
  DAG compiles with `--scoring --vep` and `--vep` alone.
