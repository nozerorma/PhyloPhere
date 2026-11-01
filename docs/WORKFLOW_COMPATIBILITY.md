# Workflow compatibility repair

This repair covers SCORING, RERconverge, FADE, SELECTION_PREP, CT_POSTPROC,
ASR_DIAGNOSTICS, CT_META_CAAS and CONTRAST_SELECTION. The VEP/MAP repair is
documented separately. Missing provenance must not be replaced with invented
scientific settings or a claim that an old run completed successfully.

## Scientific requirements retained

Observed and null inference must still share their score rule, model, filters,
strict/exploratory status and tested universe. Cycle rosters, gene and position
axes, duplicate detection, matrix values and the original multiple-testing
families remain checked. An unavailable p-value remains unavailable.

RER continues to use weighted residual estimation. A numerical failure never
triggers an unweighted retry. FADE still requires a root with two descendants:
its HyPhy implementation tests `terms.trees.rooted`, whose installed tree
reader defines this from the root's two children. The repair does not reroot
input trees automatically.

## Legacy input handling

* **RER:** inference from an unstamped matrix reconstructs weighted residuals
  from the supplied trees and current phenotype species/minimum-species setting.
  Old matrix values are discarded; an absent stamp is not proof of its estimator.
  Native repeated/empty path display labels are accepted when they match the
  native tree axes. Reconstructed matrices get unique node-pair path identifiers
  and record `recomputed_from_legacy_inputs`. Contradictory recorded identities
  and mismatched supplied axes still fail. This reconstruction costs one residual
  fit in each inference task using an unstamped matrix; building and reusing a
  current matrix avoids that cost.
* **FADE / SELECTION_PREP:** an empty foreground/background pool is a valid
  biological outcome. Empty foregrounds produce audited exclusions; malformed
  species records, overlapping pools and technical failures still fail.
* **CT_POSTPROC:** a missing observed completion record permits table validation
  and descriptive filtering, with scientific provenance unknown. Report settings
  and universe ledgers are optional when regenerating an old report from its
  saved tables. Actual removals come from the input/retained difference; their
  causes stay unknown when the original settings are absent. New inference still
  requires sufficient observed/null scientific provenance.
* **ASR_DIAGNOSTICS:** the shared VEP/MAP repair already adds descriptive reading
  of an old master table without a completion record. Model/design provenance
  remains unknown; coverage is checked against a record when it exists.
* **SCORING / CT_META_CAAS reports:** absence of original evidence metadata no
  longer prevents descriptive rendering. Original scores and p/q are retained;
  thresholds and gene statistic stay unknown. Missing thresholds do not produce
  new significance calls. Literal identifiers, including `NA`, `001` and quote
  characters, survive the significance report reader.
* **CONTRAST_SELECTION:** if `observed_selection.rds` is absent, reconstruct it
  from current staged traits, tree and settings. Composition and reconstruction
  share the same interval/category constructor, including the original Jeffreys
  interval implementation. A present bundle with contradictory provenance still
  fails. This is a new computation, not certification of the old composition.

The report regeneration GUI treats these descriptive metadata as optional while
retaining the underlying tables as required inputs.

## File integrity and cache policy

Normal execution does not repeatedly hash entire RER null blocks, residual
matrices, FADE alignments, CAAS block matrices or calibration null stores just to
prove byte identity. Scientific structure and settings are still validated.
FADE's baseline is shared only for the same resolved alignment and tree paths;
its numerical fitting and foreground assignment are unchanged.

Set `PHYLOPHERE_VERIFY_CHECKSUMS=1` for optional file integrity checks and checksum
production in the repaired producers. Existing checksum records remain readable;
missing historical digests are not represented as verified. The relevant RER,
FADE, SCORING and FCS preparation task fingerprints include the integrity setting.
Small source-code hashes remain useful for resume invalidation. The small
calibration-universe digest and RER observed/manifest identity digests are retained:
they bind scientific families and observed/null settings, rather than rescanning
large databases. Integrity hashes do not establish scientific compatibility.

## Regression checks

Run serially in the local `phylophere` environment:

```sh
conda run --no-capture-output -n phylophere python -B -m unittest discover -s checks -v
```

The checks run native weighted RER estimation and its inference/table CLIs,
native HyPhy FADE and shared-baseline reuse, complete SCORING, CT_POSTPROC and
CAAS-significance HTML renders, contrast reconstruction and interval equivalence,
FCS historical-manifest compatibility, opt-in checksum rejection and scientific
axis/statistic mismatch rejection. They use small synthetic inputs with real
installed tools. They are regression evidence, not a full cluster production run.
