# POSENRICH contract

POSENRICH consumes the same declared tested positions, cleaned gene background,
observed postproc and canonical null as SCORING. It does not change SCORING scores,
p/q values, FCS inputs or gene filters. Main/alternative entrypoints,
configuration/help, WorkflowMap, GUI/generated runners and the full ENRICHMENT
review remain separate checkpoints.

## Fixed opportunity and identities

U is the intersection of `background.output` tested zero-based trimmed alignment
columns and the cleaned gene background. It includes tested positions with no
observed or null detection. Finite non-detection does not establish structural
ineligibility. An exclusion for structural ineligibility must be declared
independently of observed/null hits. Neither observed positives nor null-positive
positions select U, and no overlap heuristic can override the tested input.
Observed rows outside U fail. Null rows outside U are excluded from its statistics.
All cleaned genes are recorded, including those with zero tested positions.
`background.output` is the positional universe of the observed design (see SCORING_CONTRACT): the null it is
checked against must carry the same per-gene universe digests, so SCORING and POSENRICH share one positional family.

Gene IDs are literal, including `001`, `NA`, dots and quotes. Postproc uses
standard CSV quoting with a tab separator; SCORING uses unquoted literal TSV.
POSENRICH publishes quoted TSV and accumulation CSV. Readers select the appropriate
format explicitly. Duplicate positional keys, fractional/negative coordinates,
unknown sides, nonfinite scores and scores outside [0,1] fail. Postproc `none`
means score zero; positive SCORING rows are top/bottom only. SCORING and postproc
positive keys/scores must agree; positional p/q/nd are checked against an audit
recomputed with SCORING’s counter functions on U and the validated cycles; the registry never reconstructs missing scores.

`positions.tsv` is a diagnostic coordinate registry, which may include annotation
or MAP positions outside U. `in_tested_universe` distinguishes inferential
opportunity from those diagnostic rows. `position_evidence.tsv` preserves sides.
Characterization and proximity explicitly pool by position using the best side;
shared position p/q are one max-over-sides test, not directional tests.

## Null identity and silent cycles

Inference requires canonical blocks, complete `cycles.txt` and `cycle_roster.txt`,
aggregate counts/positions, `scoring_artifacts.json`, and postproc
`<observed>.scoring.json`. POSENRICH reuses SCORING's checks of artifact checksums,
scientific/filter provenance, unique cycles/position-cycle keys and aggregate
consistency. Workflow validation compares global matrix maxima to position rows.
Directional null scores must be positive finite fractions when present; missing
sides are absent detections, and `all` must equal their maximum.

Valid silent cycles remain in N. An entirely silent valid null is supported;
a missing/empty roster is not treated as a negative experiment. Inferential
consumers verify that their block contents equal those validated for U. A changed
background, tested input, annotation or reference must be revalidated. Old nulls
without the current contracts must be regenerated.

## Two accumulation models

Each side has its own fixed family: all cleaned genes with at least one position
in U. Zero observed positives yield p=1 and remain in BH. Genes without positions
are listed as `no_tested_positions`, with missing p/q.

`ActualCount` is the number of distinct observed positive positions of a gene.
`all` counts a position positive on both sides once.

* `PValue`, `QValue`, `MeanNull`, `SD`: hypergeometric upper tail under uniform
  allocation without replacement of the observed total K over U, with parameters
  M=|U| and n=gene opportunity. This conditions on K and assumes exchangeable
  positional opportunities. Its q-values belong to that model.
* `EmpiricalP`, `EmpiricalQ`, `EmpiricalMeanNull`, `EmpiricalSD`: upper-tail
  comparison of the same observed count with the count in every canonical cycle.
  Total K and upstream filtering may vary by cycle. Calibration is
  `(1 + cycles with count >= observed)/(N + 1)`.

The empirical count test differs from SCORING's gene test, whose statistic is
the raw maximum score. These models answer different questions. Their q-values do not control a union
of discoveries selected across models; such a union needs its own predeclared
multiple-testing policy. Neither
establishes positive selection, independent origins or causality.

Adding positions that are zero in the observed and every cycle changes the
hypergeometric opportunity. It does not by itself change count or maximum-score
exceedances if genes, scores, cycles and filters are held fixed. It can change
BH families. Within-gene centered statistics can change because their opportunity
ratios change. These effects are checked jointly with SCORING rather than
assuming that every p must change or that every result must remain constant.

## Annotation coordinates, references and coverage

MAP selected columns must be unique and contiguous; raw columns are explicit.
Raw/trimmed alignment projection is separate from biological reference transfer.
Protein/genomic coordinates and strand require the same verified `.vep.json`
reference contract as VEP. A missing sidecar is `unverified_reference`, with no
inferred strand or transferred reference coordinates. Malformed MAPs/references
fail. Protein-space source annotations additionally require matching
`reference_identities` per literal gene in their source manifest. Removed raw
columns are audited exclusions, not projected to neighboring columns.

Pfam uses declared direct CAAS columns or verified protein coordinates. UCR uses
explicit raw/trimmed/CAAS coordinates and the sliding method. FUBAR requires
explicit coordinates, records all tested sites including neutral calls and
rejects contradictory duplicate calls. FADE consumes its upstream selected-site
CSV without inventing a new BF cutoff. Custom TSV uses literal zero-based columns.
Overlapping memberships are sets, not repeated statistical units.

A source may be descriptive without supplying enough provenance for inference.
For within-gene use, `<source>.manifest.json` must checksum the complete table,
identify `source_uri` and `version` (or a verified projection version), bind the
productive alignment inventory or both `tested_sha256`/`background_sha256`, and
declare coverage independently of hits. Supported coverage declarations:

* `all_productive_alignment_columns` or `tested_positions_in_declared_genes`,
  with `coverage_genes` explicitly including covered genes with no hits;
* `listed_tested_positions` for complete FUBAR tested-site tables.

`coverage_positions` can restrict each gene to explicit zero-based scanned
columns. Pfam generation publishes its exact query columns, excluding reference
gaps/X. Verified protein-space Pfam coverage is restricted to reference residues.
Absence of hits is not inferred to establish coverage. Sources lacking a verified
coverage contract are marked `unverified_source` and excluded from inference;
their descriptive annotations remain visible. `annotation_coverage.tsv` stores
compact covered intervals separately from positive memberships.

Pfam autogeneration and first-use cache preparation are preserved. Generated
snapshots record HMM/clan content hashes, source URI, cutoffs and exact query
coverage. Cache preparation is serialized and downloads are atomic. Prepared
resource contents and transitive helpers enter the task fingerprint. UCR
regeneration records its productive input inventory and algorithm parameters.

## Within-gene concentration and reports

For each term/gene pair, n is the number of source-covered positions in U, m its
members and k the positive positions in that covered opportunity. Its centered
residual is `positive_members - k*m/n`. The term sums residuals across eligible
genes, which require both members and nonmembers. `min_genes` counts genes
actually supplying such a contrast. Eligibility/families depend on declared
annotation and coverage, never observed positives. Empty or unidentifiable
terms have explicit states and no p/q.

Observed and null cycles use the same residual function, each with its own k.
This is centering, not exact conditioning on observed gene counts; genes with
more detections can contribute larger residuals. Empirical upper-tail calibration
uses every valid cycle and inclusive ties, with BH separately within each
collection. Holdout cycles are excluded from the calibration denominator.
Related Pfam/UCR/FUBAR/FADE evidence is not independent external validation.

Proximity remains descriptive, measures reference residues separately from
nucleotide anchors, and reports missing reference coverage rather than infinite
distances. Characterization denominators distinguish tested opportunity from
positive-score availability; missing annotation is not a biological negative.

Manifests checksum inputs and output tables. Rendering verifies snapshots and
uses producer p/q; it does not rerun inference or select a new BH family.
The workflow publishes `report_settings.json`; regeneration recovers the saved
BH display threshold unless the caller explicitly overrides it.
Nextflow uses deep content caches and explicit source fingerprints. No GUI or
runner regeneration redesign is included in this checkpoint.
