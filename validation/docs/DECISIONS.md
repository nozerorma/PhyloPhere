# Validation suite — decision log

Running log of design decisions, so they can be reviewed and reversed. Newest
first. "ad-hoc" = made by Claude under general direction; "MIGUEL" = explicitly
chosen by Miguel.

---

## Tier 0

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

### D-T0-02 — MIGUEL (REVISED): one planted mechanism — `identical_aa` only
`identical_aa` = same-residue convergence (Zou & Zhang / caastools US), hard-set
on each origin edge, held by a near-point-mass Q_fg. Target = the least-favoured
residue under the site's background Dirichlet profile, so the signal is
unambiguous. It is also a grouped-CAAP (GS1-GS4) hit for free (shared residue ⊂
shared group). `profile_shift` (the earlier PCOC-style preference shift) is
**stripped** — Miguel: it is neither a strict-CAAS positive nor a clean CAAP
positive because CAAP schemes are *grouped* on `amino_encoded`, and a random
alternate profile does not converge to a shared physicochemical group. It only
inflated the recall denominator with sites CAAS was never meant to call. The
`SimConfig.n_planted_profile_shift` hook stays (dead-defaulted to 0).

### D-T0-N — deferred: a real grouped-CAAP planted mechanism
To exercise GS1-GS4 as more than a trivial by-product of `identical_aa`: on fg
edges force the equilibrium onto a shared *amino_encoded group* (e.g. hydrophobic)
but leave the residue free within it, so US sees divergent residues while the
grouped schemes see convergence. Would let `score` separate US recall from GS
recall. Needs the pipeline's group definitions (the `amino_encoded` recoding
tables) mirrored in `model.py`. Not built — reopen if CAAP-scheme calibration
becomes a Tier 0 goal.

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

`harness.cli calibrate --pcol` chooses the column: `meta_caas_boot` (default,
gate — `US_meta_caas.tsv:pvalue_boot`), `meta_caas_hyp`, `perm_pos_boot`
(`perm_pos_pval.tsv:null_pvalue_boot`), `position_boot`, `position_hyp`. Pooled
across every null context = whole null-set replicates **plus** the unplanted
genes of power replicates (primary source = the null-set replicates; the
unplanted-gene tail is a secondary top-up). `perm_pos_boot` on the smoke run is
degenerate-discrete (`n_cycles` 50) — needs a real null-set run to read.

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
