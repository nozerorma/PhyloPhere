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
`asr_model = "lg"`. Scoring / accumulation / enrichment / VEP stay OFF.

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

## Open questions (do NOT resolve unilaterally)

- **Q-T0-A — scoring truth when operative fg/bg ≠ planted design.** After
  PIPELINE_MODEL.md §6.2 the gap looks small (contrast selection forms N
  independent pairs, does not drop clades; planting on clade stems makes the
  representative choice irrelevant). Proposal: score each module against *its
  own* operative foreground, and additionally report contrast selection's
  recovery of the planted origins as a first-class Tier 0 metric. Needs Miguel's
  sign-off.
- **Q-T0-B — `perm_strategy` for Tier 0.** BM vs lambda vs FGBG. lambda accepts
  far more designs; BM is the stricter default. Which is the calibration target?
- **Q-T0-C — continuous trait encoding for the `echo`-like binary archetype.**
  Feed raw 0/1 (RER auto-detects binary; contrast selection discretises trivially;
  FADE quantile-extremes on a 0/1 vector is degenerate) vs feed a bimodal
  near-0/near-1 continuous vector. Affects FADE's foreground most.
