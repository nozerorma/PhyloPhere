# Validation suite — decision log

Running log of design decisions, so they can be reviewed and reversed. Newest
first. "ad-hoc" = made by Claude under general direction; "MIGUEL" = explicitly
chosen by Miguel.

---

## Tier 0

### D-T0-Q — evolutionary-model planting (2026-08-31, MIGUEL, supersedes the D-T0-P foreground)
The D-T0-P "farthest-point anchors + engineered trait" made contrast recovery a
tautology (we cherry-picked independent lineages and built the trait around them).
Replaced with:

- **Latent trait = BM on the tree rescaled by Pagel's lambda** (`pheno._lambda_rescale`
  + `_bm_latent`) — the same family the permulation's `simpermvec` draws its nulls
  from. lambda=0 star -> scattered extremes; lambda=1 full BM -> clade-clumped;
  lambda=0.5 ~ the real cancer trait (fitted ~0.74). **lambda in {0, 0.5, 1}** is a
  gate axis (`--lambdas`).
- **True foreground** = the independent top/bottom tail origins of that latent
  (`_independent_pairs`, Dunn >= 1). MIGUEL: the number **used is FIXED**
  (`n_pairs = 4`, matching the pipeline's `contrast_max_iter = 3`) so replicates
  are comparable across lambda — the latent is resampled until it supports >= 4.
  `notes['n_possible_pairs']` (how many it *could* support: ~20 at lambda 0, ~5
  at lambda 1) is REPORTED as the structure diagnostic.
- **Observed trait** = latent + sampling noise: binary threshold + flipped bits;
  rate `c ~ Binomial(n, rate(latent_pctile))`, n ~ U(25,70).
- **Plant on the TRUE lineages, not the operative pairs.** Contrast recovery
  (operative `traitfile.tab` vs true pairs) is now an honest lambda-dependent
  measurement: anchor recall / partner recall / exact-pair / lineage Jaccard,
  reported as a curve, NOT gated.
- **Top vs bottom**: each planted gene planted in one direction; the pipeline
  scores `change_top`/`change_bottom`, `gene_caas_score_top/bottom` symmetrically.
  Adapter reports `site_directional_recall` (detected AND on the right side).
- **Score-level separation** reported (not gated): position `CAAS_score` AUC,
  position `-log10 pvalue` AUC (the hypergeometric p separates strongly),
  `gene_caas_score` AUC (does not — size-decorrelated by design).
- VERDICT: occurrence separation at every cell + precision@k / site_recall at the
  lowest lambda.

### D-T0-P — FINAL design (2026-08-30, supersedes D-T0-A / D-T0-C / D-T0-M / D-T0-H)
After inspecting the real cancer run and the permulation code, the Tier 0 design
was rebuilt:

- **The pipeline is a prioritisation engine.** Its permulation p-values
  (`pvalue_boot`, `null_pvalue_boot`) are conditional foreground-specificity
  scores — 0% exact zero on the real run, but ~78% ≤ 0.05 because discovery
  pre-selects positions that align with the foreground. They are *never*
  Uniform(0,1), on a real trait or a null. `phen_score = 1 - percent_rank(...)`
  is a well-behaved uniform *ranking* input (verified against the real
  `position_scores.tsv` to machine precision). The FCS `p.perm` (Wilcoxon-AUC vs
  the 1000 Dunn-resample columns) *is* empirically uniform — the one calibrated
  null in the system. See [[project_tier0_scoping_reframe]].
- **Drop the KS `calibrate` command.** No "null permulation p ~ U(0,1)" test.
- **Foreground = `n_pairs` (default 4) explicit contrast pairs.** Real run =
  3 pairs from `traitfile.tab`, NOT 10v10 (`top_species.txt`/`bottom_species.txt`
  are FADE's percentilised proxy). `pheno.make_paired_foreground`: farthest-point
  anchor tips + nearest bg neighbour each; planted signal on the anchor terminal
  edges; contrast selection recovers those pairs.
- **Two archetypes** (D-T0-04 kept, reshaped): `binary` (0/1 code,
  `--trait_type ordinal`, CLASS 2) and `rate` (`c/n` + `n_pop`/`n_cases` count
  columns, CLASS 1 Jeffreys-CI). `echo`/`bodysize` (BM-threshold, clade-sample)
  are gone. `binary` needs the binary-trait contrast-selection fix (landed on
  `main` `1fbdfbf`).
- **Scoring** (`harness score`): (1) prioritisation — planted in
  `slice_global1`/`slice_global5`; (2) **null-vs-power separation** — AUC of
  planted-gene `gene_caas_score` (power) vs every gene's score (matched null),
  `separated` = AUC ≥ 0.9 AND planted p50 > null p95; (3) site recall by
  mechanism × scheme (`identical_aa→US`, `grouped_caap→GS`, US leakage);
  (4) contrast recovery (`traitfile.tab` fg vs planted pairs).
- Scale: ~120-150 genes, ~350-400 sites, `concentration ≈ 2`, production
  permulation settings. Gate = `primate` × {binary, rate} × {null, power}.

**GATE RESULT (2026-08-31, `runs/tier0_gate2`, 10 reps/cell, 120 genes): PASS.**
- null-vs-power separation AUC(detected-CAAS count) = **1.00** both archetypes
  (planted genes carry ~24 CAAS positions, null-replicate genes ~1).
- precision@k 0.83 (binary) / 0.95 (rate); site recall 0.99; `identical_aa→US`
  1.00; contrast recovery Jaccard 1.00, all pairs.
- **`gene_caas_score` does NOT separate across replicates** (`caas_score AUC
  ≈ 0.45`) — it's a within-run percentile rank, size-decorrelated by design;
  detected-CAAS count is the separating quantity. So `score`'s separation metric
  is on `n_positions`, not `gene_caas_score`.
- **`grouped_caap` also trips US ~86%**: with `patterns=1,2,3` + conserved
  disjoint bg, class-shared/residue-divergent fg satisfies US pattern 3. Clean
  GS-only positive needs a class-heterogeneous bg — deferred. GS recall (0.99)
  still > US (0.86).
- `slice_global5` membership is ~0.06 (mechanically capped: 20% planted vs 5%
  slice); `slice_global25` ≈ 0.26 ≈ its theoretical max. Not a gate criterion.

### D-T0-06 — pipeline studied before further integration design  (2026-08-30)
After D-T0-05 it became clear the integration was being designed on a wrong model
of contrast selection / the permulation null. Stopped, read the CAAS / contrast /
permulation / RER / FADE / disambiguation code + a real result, wrote
`docs/PIPELINE_MODEL.md`. Integration design resumes from there.

### D-T0-05 — MIGUEL: module set = contrast_selection + CAAS + FADE + RER + **disambiguation**
Disambiguation (and bundled CT_POSTPROC + ASR robustness) is IN, to bring the
`asr_path_score` into Tier 0 coverage (the star/ladder pathological trees).
Needs `ct_disambig_asr_mode = "compute"` (no precomputed ASR for simulated data),
`asr_model = "lg"`.
**SUPERSEDED IN PART by Q-T0-D**: originally said "Scoring / accumulation /
enrichment / VEP stay OFF", but PIPELINE_MODEL.md §6 shows SCORING is what
produces the gene/position rankings — without it Tier 0 has nothing to score
power against. VEP stays off. Accumulation / enrichment / scoring: see Q-T0-D.

### D-T0-04 — MIGUEL: two phenotype archetypes, both through contrast selection
- `echo`-like: binary 1/0, contrast selection extracts one pair per non-nested `1`.
- `bodysize`-like: continuous, percentilized into contrasts.
Both are separate replicate sub-sets. Rationale: mirror the two archetypal
PhyloPhere analyses. contrast selection is central and must always be in the loop
(MIGUEL, emphasised twice).

### D-T0-03 — ad-hoc: generative model = WAG exchangeabilities + site-specific Dirichlet profiles + Gillespie
Deliberately not LG (pipeline ASR/FADE default). Site heterogeneity (Dirichlet
`concentration` knob) is the property single-profile inference lacks and is what
makes background convergence arise. `data/{wag,dayhoff,lg}.dat` vendored from the
cluster PAML 4.10.10. Alternative considered: `iqtree2 --alisim` (rejected — less
control over the profile-shift mechanism, and would share the inference tool).

### D-T0-02 — MIGUEL (REVISED): two planted mechanisms — `identical_aa` + `grouped_caap`
- `identical_aa` = same-residue convergence (caastools US), hard-set on each
  origin edge, held by a near-point-mass Q_fg. Target = the least-favoured
  residue under the site's background profile. A US **and** GS1-GS4 positive.
- `grouped_caap` = same-*class* convergence (D-T0-N, now built). On fg edges the
  equilibrium is ~uniform over one GS class (default GS1), the class chosen as
  the one most disfavoured in the site background; a fresh residue from the class
  is seeded per origin edge, so origins share the class but usually differ in
  residue. A GS-only positive — US should **not** call it. This exercises what
  distinguishes the grouped schemes; Miguel: "that's more important".
- `profile_shift` stays **stripped** (neither a US nor a clean CAAP positive —
  a random alternate profile has no shared class). `SimConfig` hook dead at 0.

Default load: 12 `identical_aa` + 12 `grouped_caap` per planted gene.
GS class tables vendored in `validation/tier0/groups.py` from the pipeline's
`biochem/grouping.py` (keep in sync).

### D-T0-O — MIGUEL: separate the gate from the characterisation
Miguel: "isn't it too much for tier0". The **gate** (pass/fail, what certifies
the pipeline) is one config: `primate` tree × {echo, bodysize} × {null, power}
at `min_divergent_fraction 0.5` / `planted_fraction 0.25`, 20 reps. That is the
`run_replicates.py` default. Everything else — `primate_x5` / `mammal` /
`star` / `ladder`, and the `0.5,1.0` × `0.1,0.25` axis sweep — is **optional
characterisation**, run only if the gate passes and we want the degradation
curves for the paper. Kept in `grid.json` / behind flags, off by default.

### D-T0-N — DONE (see D-T0-02): grouped-CAAP planted mechanism built
`score` now reports `identical_aa_us_recall` (planted same-residue sites called
by US), `grouped_caap_gs_recall` (planted same-class sites called by any GS,
read from `position_scores.tsv:scheme_set`), and `grouped_caap_us_leakage` (of
those, the fraction that also tripped US — want ~0).

### D-T0-01 — MIGUEL: trees  (REVISED — full primate tree, no prune)
`primates_233_subst.tree` used **whole** (237 tips = 233 primates + 4 outgroups,
total length 1.30) — exactly what production runs on. The earlier 50-tip prune
left production contrast selection unable to form ≥ `min_contrasts` (3)
Dunn-independent pairs (50 tips is too sparse for the independent-contrasts
algorithm on a shallow tree; the run skipped at CHECK_MIN_CONTRASTS). Miguel:
disambiguation on 233 species is ~minutes at Tier 0 gene counts, no reason to
prune. `primate_x5` = same tree, branch lengths ×5 (**synthetic**, non-degenerate
power curve). `mammal` = `speciesTree_speciesname_pruned.nh` pruned to 60 (deep
null stress). `star`/`ladder` (in code) = ASR path-score edge cases.
`sample_foreground` auto-excludes the 4 outgroups (Mus/Oryctolagus/Tupaia/
Galeopterus — terminal branch > 5× the 95th pct) from being foreground.

---

### D-T0-D — MIGUEL: Tier 0 module set, FINAL (revised — RER/FADE dropped)
`contrast_selection` + `CAAS` (discovery, resample, bootstrap, permulation) +
`CT_DISAMBIGUATION` (+ CT_POSTPROC + ASR robustness) + `SCORING`.
**OFF**: `RER` + `FADE` ("not much need to validate rer and fade" — published,
author-validated methods; Tier 0's target is the novel machinery: CAAS discovery
+ permulation + disambiguation / ASR-path-score + SCORING integration),
`CT_ACCUMULATION` (little meaning at low gene count), `ENRICHMENT` / STRING /
DOMINO (Tier 2), `VEP`, `reporting` (exploration Rmds assume real taxonomy;
contrast_selection runs its own dataset_exploration + emits trait_stats anyway).
Consequences: SCORING runs with just the CAAS axis (`gene_caas_score`); its
`rer_*` / `fade_*` / `accum_*` columns come out NA (full_join + file_exists
guards). Removes the RER gene-trees problem (no ML tree builder in the env) and
supersedes D-T0-F (rer_transform).

### D-T0-A — ad-hoc (reversible): scoring truth = operative foreground per module
Score each module against the foreground **it actually used** (CAAS: contrast
pairs from `traitfile.tab`; RER: `value==1` / the continuous vector; FADE:
quantile-extreme species), read from the run's own outputs. Additionally emit a
**contrast-selection recovery** metric per replicate: precision/recall of the
planted origins vs the pairs `traitfile.tab` contains, plus the `low_contrasts.skip`
rate. End-to-end recall = product of the two, recoverable but not the headline.
Rationale in PIPELINE_MODEL.md §6.4 / §7.2 — the gap is expected to be small
because planting is on clade stems and origins are non-nested.

*Implemented* (`validation/harness/tier0_adapter.py`, `harness.cli score`): per
power replicate it emits gene recovery (`score_metrics` on `gene_caas_score`:
ROC/PR-AUC, precision@k, planted ranks), planted-**site** detection
precision/recall split by mechanism (SCORING only lists *detected* positions, so
no full ROC there), and contrast recovery = Jaccard(operative fg from
`traitfile.tab` `label==1`, `truth.phenotype.foreground_tips`) + fraction of
pair ids whose fg member is a true origin. RER/FADE arms deferred with those
modules. Smoke run: gene precision@k 1.0, ranks 1-4 exact; site recall
`identical_aa` 1.0 / `profile_shift` 0.0 (CAAS is by construction the
identical-residue detector); contrast Jaccard 0.8 (one origin's sister tip
picked instead of the origin tip).

### D-T0-M — ad-hoc: what the null set actually tests + selectable p-value source
The pipeline's null **is** the BM permulation (D-T0-B): `CAAS_PERMS_DISAMBIGUATE`
replays N phenotype relabelings through detection + ASR, and `pvalue_boot` /
`null_pvalue_boot` are that permulation p-value. So the Tier 0 null set (no
planted signal, random foreground) is a *calibration check on the permulation*:
when nothing convergent is associated with the trait, are those permulation
p-values ~ U(0,1)? That is exactly `null_calibration` (KS uniformity + observed
type-I at nominal α). It certifies **trait-association** type-I control; it does
not by itself certify robustness to profile-heterogeneity background convergence
except insofar as those chance convergences sometimes align with the random
foreground across many replicates (which the KS test integrates over).

`harness.cli calibrate --pcol` chooses the column, now just two:
`meta_caas_boot` (default, **the gate** — `US_meta_caas.tsv:pvalue_boot`, the
permulation p that drives the significance call) and `perm_pos_boot`
(`perm_pos_pval.tsv:null_pvalue_boot`, the raw per-position permulation p before
meta-aggregation and FDR — diagnostic, shows the null before smoothing).
Dropped: `meta_caas_hyp` / `position_hyp` (the analytic hypergeometric p is
anticonservative by design — testing it for uniformity is uninformative) and
`position_boot` (just `meta_caas_boot` carried into SCORING). Pooled across every
null context = whole null-set replicates **plus** the unplanted genes of power
replicates (primary = the null-set replicates). `perm_pos_boot` is
degenerate-discrete at low `n_cycles`.

### D-T0-B — MIGUEL: `perm_strategy = BM` only
BM is the shipped default and matches RER's own BM-only null (`getPermsContinuous`),
so the two modules' nulls stay comparable. lambda/FGBG not run.

### D-T0-F — MIGUEL: `rer_transform = "auto"`
echo runs RER binary (no transform). bodysize is a real-valued BM trait (~normal,
centred at 0), so `auto` -> raw values. (log10 is the natural allometric transform
but N/A for a zero-centred BM trait unless we simulate exp(BM); noted, not done.)

### D-T0-G — MIGUEL: CAAS detection = production defaults **+ strict-divergent variant**
Production: `caap_mode` ON (US + GS1/GS2/GS3), `patterns` 1,2,3,
`ct_disambig_convergence_mode` = focal_clade, `fade_method` = Variational-Bayes
(FADE config advanced defaults). Run each replicate at **two**
`min_divergent_fraction` values as separate sub-sets:
- `0.5` (production; `max_conserved` = half the pairs)
- `1.0` (strict; `max_conserved` = 0, no pair may share the residue)
`run_replicates.py --divergent-fractions 0.5,1.0`. Nothing else changes between
the two.

### D-T0-E — MIGUEL: relax the `& trait </> median` extreme-filter guard (option A)
The strict median guard makes a minority-foreground binary/bimodal trait
un-labelable (median == lower threshold == 0, nothing can be `< median`).
Relaxed `<` -> `<=` and `>` -> `>=` at the median in **two** places, both no-ops
for continuous traits (upper_thresh >= median >= lower_thresh always), both
committed on `validation`:
1. `lean_contrast_selector.R:143-144` — else `low_sp` empty -> CAAS permulation
   pool cannot fill.
2. `subworkflows/TRAIT_ANALYSIS/local/src/stats.R:188-189` (`global_label`) and
   `206-207` (`taxa_label`) — else `trait_stats.csv` has no `low_extreme` and
   FADE's `EXTRACT_EXTREME_SPECIES` crashes ("No bottom species found").
Production `3.CI-composition.Rmd` categorises with no median gate, so this brings
the three into line. **Coordinate**: `perms_lambda` branch has its own rewrite of
`lean_contrast_selector.R` that keeps the strict form — reapply there before merge.

### D-T0-C — MIGUEL-approved (via D-T0-E): trait encoding per archetype
- `echo` archetype: **exact 0/1** `--my_traits` column. RER auto-detects binary
  and runs `foreground2Tree` / the binary path. Contrast selection's production
  discrete path categorises trivially; the CAAS permulation works because of the
  D-T0-E fix. `top_quantile` / `bottom_quantile` set so the cut lands on the 1s
  (top) / 0s (bottom), which is also FADE's `EXTRACT_EXTREME_SPECIES` input.
- `bodysize` archetype: a **simulated Brownian-motion continuous trait** on the
  tree; foreground = its top-quantile tips (emergent, not pre-sampled); the
  molecular signal is planted on exactly those tips' lineages. RER continuous.
  DESIGN doc's "BM phenotype with known λ" cell.
- Per replicate, `top_quantile` / `bottom_quantile` are set so FADE's
  `EXTRACT_EXTREME_SPECIES` top set and contrast selection's high set both land
  on the operative foreground.

### D-T0-H — MIGUEL: planted-gene fraction is a sweep axis (0.1 / 0.25)
`run_replicates.py --planted-fractions 0.1,0.25`. Applies only to the `power` set
(null has no planted genes). Checks whether null calibration is sensitive to how
much real signal is mixed into a replicate. Sub-set dir gets a `pfNN` tag.

### D-T0-I — MIGUEL: disambiguation cost is not a concern
`ct_disambig_asr_mode = "compute"` on every replicate. Miguel: disambiguation on
the full 16k-primate / 233-species production set runs in ~6h on the cluster, so
the Tier 0 volume (tens of genes × tens of species × ~hundreds of replicates) is
comfortably within budget. Earlier "days" estimate was wrong. No gene/site
shrinking, no disambiguation subsetting.

### D-T0-J — MIGUEL: first run = local smoke
2 replicates per sub-set, `--trees primate`, small `n_genes` / `n_sites`, on the
laptop (native `phylophere` micromamba env, `-profile local`, no containers).
Goal: shake out pipeline wiring / crashes, not calibration. Then size up on the
cluster.

---

### D-T0-K — ad-hoc: synthetic support files per replicate + `gene_filter_mode = dubious`
Real primate runs always ship inputs the synthetic data lacks; `build_replicate`
now writes them: `ali_sp_names.txt`, `gene_ensembl.tsv` (fake coords, 3 fake chrs,
length = n_sites·3), `taxid.tsv` (fake taxids 9000001+; disambiguation renames
tips to numeric ids before PAML, dodging `_write_phylip`'s single-space
name/sequence delimiter — a latent PAML-format bug, flagged), `family` column in
`my_traits.tsv`. `gene_filter_mode = "dubious"` (production default, NOT "none":
the `filtered_discovery.tsv` emitter is `when: gene_filter_mode != 'none'` and
SCORING starves without it; "dubious" only drops IQR-outlier genes with
spatially-clustered CAAS, which planted genes are not).

### D-T0-L — MIGUEL: fix the `asr_mode=compute` permulation race in the PIPELINE
Miguel: "ideally on compute, caas_perms_disambiguate should run after, reusing
the computed asr." Done (commit `8d8b7e7`):
- `caas_permulation.nf` `CAAS_PERMULATION` gains an `asr_ready` gate; the replay's
  tree input is held via `.combine(asr_ready)` until it emits.
- `main.nf` passes `disambiguation_results.master_csv` as the gate when
  `ct_disambiguation && asr_mode==compute`; `NO_GATE` sentinel otherwise.
- `gene_wrapper.py`: the perms replay now WARNs per uncached gene (was silent),
  WARNs on partial coverage, ERRORs on a zero-row pass A.
Follow-up (commit `fe5279b`, Miguel-approved): `_load_gene_asr_context` now uses
`run_asr_pipeline(skip_if_exists=True)` — loads the cache on a hit, **computes
ASR into the cache on a miss**. So the perms replay is robust to null-only genes
(appear only under a shuffled labeling, never observed-significant, never cached
by the main disambiguation). Gate + compute-on-miss together: the gate ensures
the observed genes are cached first (no duplicate work / no race), compute-on-miss
handles the null-only tail. `asr_mode=compute` now yields a **complete**
permulation null. Verified on the smoke: null went from 4 → 9 contributing genes
(4 planted + 5 null-only), `asr_cache/` ends with all 15 genes, `perm_pos_pval`
30 → 208 rows, `null_pvalue_boot` no longer degenerate-flat.

## Pipeline bugs surfaced by Tier 0 (all on `validation`)
- `lean_contrast_selector.R:143` + `stats.R:188,206` — `& trait </> median` guard
  breaks minority-fg binary traits (D-T0-E, FIXED, Miguel-approved).
- `disambiguation_main.py:200` — gene name = `GenePos.split("_")[0]`, breaks on
  any id containing `_`; the `Gene` column exists and should be used. Worked
  around in Tier 0 (gene ids `gNNNN`); real fix pending Miguel.
- `reconstruct.py:388` (`_write_phylip`) — single space between PAML name and
  sequence; PAML needs 2+. Latent (production renames tips to numeric taxids).
  Worked around in Tier 0 via `taxid.tsv`; real fix pending.

## Still genuinely open

*(none blocking the build — reopen any D-T0-A/B/C if the first calibration run
shows they matter)*
