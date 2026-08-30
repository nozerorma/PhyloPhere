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

### D-T0-02 — ad-hoc: two planted-positive mechanisms, scored separately
`profile_shift` (PCOC-style preference shift on fg edges) and `identical_aa`
(Zou & Zhang same-residue convergence, seeded at each origin edge). Recorded
separately in `truth.json` because CAAS strict mode and CAAP/FADE target
different definitions.

### D-T0-01 — MIGUEL: trees
`primates_233_subst.tree` (233 tips, total length 1.30) pruned to 50 by
farthest-point = the regime actually run + `primate_x5` (branch lengths ×5,
**synthetic**, for a non-degenerate power curve) + `speciesTree_speciesname_pruned.nh`
(719 tips, length 29.3) pruned to 50 = deep-tree null stress + star/ladder
(built in code) = ASR path-score edge cases. Fetched from cluster
`external_phylogenies/`, gitignored.

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
