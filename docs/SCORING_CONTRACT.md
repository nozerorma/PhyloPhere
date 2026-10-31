# SCORING contract

SCORING preserves literal gene IDs (including `001`, `NA`, internal dots and quotes).
CT postproc tables use standard CSV quoting with a tab separator; FADE, RER and
SCORING tables use literal, unquoted TSV fields. Tabs and newlines are forbidden
inside identifiers. Gene tables have one row per gene; position tables have one
row per `(Gene, Position, side)`. Positions are zero based.

The position statistic is the upstream domain-pooled `caas_pos_species_pss_v4`
score. The pooled position test uses the best side. The gene statistic is its
raw maximum, globally or within a side. Descriptive residue groupings do not
enter either statistic. Neither score is a calibrated probability of independent
origins. RER correlation sign is independent of CAAS top/bottom; FADE maximum BF
is descriptive site evidence, not a calibrated gene-wide test.

Positional universe. `background.output` lists, per gene, the positions that b_0 examined and that pass the CAAS test
(patterns, `max_conserved`, `min_non_overlapping`, gap, missing-species and `miss_pair` limits of the run) under at least
one orientation of the pairs of at least one observed hypothesis, i.e. with the members of each pair swapped between
foreground and background. No pair is added or re-matched, so the Dunn index of each hypothesis is unchanged. The
universe is fixed by the observed design and is also the only set the null scores: each null cycle drops positions
outside it before the gene maximum and the position counts are taken (detections still count for gene removal and
cluster trains). Genes with no universe position keep a `NULL` line in `background.output` and are not in `background_genes.output`, so the
gene family (that list minus the genes of the gene filter, `cleaned_background_main.txt`) is the genes that have positions to
test and every consumer reads the same family. The
strict analysis has its own universe, computed with `max_conserved` 0. Null tiles record a digest of every gene's
universe, and SCORING rejects a null whose digests differ from `background.output` or are absent.

Position calibration requires the declared tested positions (`background.output`)
and cleaned gene background. Their intersection fixes the BH family before
looking at observed/null detections. All observed positives must belong to it.
Null-only positions outside that family are excluded. Observed-undetected members
contribute p=1 to BH without materializing a large vector. Never-redetected
observed positives retain `(1+k)/(1+N)` with k=0. Valid silent cycles remain in N;
a valid entirely silent null is distinct from an absent/failed null. The old
50% overlap and observation-count guards have been removed.

Gene calibration uses the complete cleaned background, zeros for valid
non-detections, three separate BH families (global/top/bottom), all canonical
cycles together and native leave-one-out null ranks. Missing evidence is exposed
by availability fields and missing p/q; no analytic substitute is used.
Disabling pathway enrichment does not disable SCORING calibration.

Canonical block artifacts carry `scoring_artifacts.json` content checksums. The
validated input identities are published under `scoring/scoring_inputs/`.

New inference from a saved null requires canonical blocks and
`filtered_discovery.tsv.scoring.json`, produced by CT_POSTPROC. Null tile/merge
records carry model, tree, weights, taxon mapping, cluster and gene-filter
settings, MAP/alignment hashes and length values. SCORING verifies the observed
checksum and compatibility of these records. It checks unique cycles and
position/cycle keys, reconstructs aggregate detection counts, compares position
maxima with global matrices, and enforces `0 <= k_emp <= nd <= N`.
Legacy nulls without scientific/filter provenance must be regenerated; matching
a score-rule string alone does not establish compatibility. The public workflow
signature and downstream table/list channels are retained.

Calculation, counters, calibration and rendering use content-sensitive caches
and source fingerprints. Report metadata records table hashes and effective
thresholds. Rendering verifies its input snapshots and preserves the producer's
RER empirical BH family. The compatibility gene-correlation table is header-only.
The GUI regeneration change only stages the shared report style; full GUI,
configuration/help, main-entry and WorkflowMap reviews remain separate.
