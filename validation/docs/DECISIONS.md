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

### D-T0-B — ad-hoc (reversible): `perm_strategy = BM` for the primary Tier 0
BM is the shipped default (`CaasConfig.perm_strategy = "BM"`). `lambda` (the
`perms_lambda` branch, uncommitted; λ̂≈0.74, higher acceptance) and `FGBG` are
added as variant sub-sets once BM is calibrated. RER's own null is BM-only
(`getPermsContinuous`), so BM keeps the two modules' nulls comparable.

### D-T0-C — ad-hoc (reversible): trait encoding per archetype
- `echo` archetype: **bimodal continuous**, NOT exact 0/1. Foreground tips get
  U(0.85, 1.0), background U(0.0, 0.15), clear gap. Reason: the lean permulation
  selector (`lean_contrast_selector.R`) filters extremes with `trait > median`
  and `trait < median`; a minority-foreground 0/1 vector has median 0 so its low
  group is empty and `permulations.R` cannot fill the pool → the run fails. The
  bimodal encoding is morally binary for contrast selection / FADE but keeps the
  median split well-defined. **Cost**: RER runs continuous (not binary) for this
  archetype, so the RER binary-mode path is not exercised by Tier 0.
  **Flag for Miguel**: the `& trait_vec < th$med` / `& trait_vec > th$med`
  clauses in `lean_contrast_selector.R::evaluate_lean_contrast_selection` look
  like they make a genuinely binary trait un-permulable in production too — the
  production `3.CI-composition.Rmd` discrete path has no such median clause. Not
  fixed here; noted.
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
