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

### D-T0-D — MIGUEL: Tier 0 module set, final
`contrast_selection` + `CAAS` (discovery, resample, bootstrap, permulation) +
`CT_DISAMBIGUATION` (+ CT_POSTPROC + ASR robustness) + `RER` + `FADE` + `SCORING`.
**OFF**: `CT_ACCUMULATION` (little meaning at low gene count — revisit only if we
decide to validate accumulation itself), `ENRICHMENT` / STRING / DOMINO (pathway
level is Tier 2's job, "does not belong to Tier 0"), `VEP`.
Consequence: `scoring_compute.R` handles absent `accumulation_all_*` gracefully
(memory: SCORING Issue 1 fix), so SCORING runs fine with accumulation off — the
`accum_*` gene columns are just NA.

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

### D-T0-E — MIGUEL: fix `lean_contrast_selector.R` (option A)
The `& trait_vec > th$med` / `& trait_vec < th$med` extreme filters relaxed to
`>= th$med` / `<= th$med`. A minority-foreground 0/1 vector has `median == lower
== 0`, so the strict form left `low_sp` empty and the CAAS permulation pool could
not fill. Production `3.CI-composition.Rmd` has no median gate; this brings the
lean selector into line and is a no-op for any continuous trait. Committed on
`validation`. **Coordinate**: `perms_lambda` branch has its own rewrite of this
file that keeps the strict form — the same relax must be applied there before
that branch merges, or it will regress.

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

---

## Still genuinely open

*(none blocking the build — reopen any D-T0-A/B/C if the first calibration run
shows they matter)*
