# VEP characterization contract

VEP accepts exactly one SCORING `position_scores.tsv`. Literal unquoted TSV IDs,
including `001`, `NA`, dots and quotes, are preserved. Required columns are
`Gene, Position, side, CAAS_score, caas, ancestral_aa, derived_aa`. Positions are
zero-based trimmed alignment columns, side is top/bottom, positive scores are
finite and at most one. Duplicate keys, malformed rows and incompatible schemas
fail before annotation. `Position + 1` indexes MAP `prot_ali_col`.

## Scientific scope

SCORING pools ancestral and derived states across domains. These sets do **not**
recover domain-paired substitutions. VEP therefore labels all transfers
`human_reference_candidate` / `pooled_ASR_state_candidate`. A human reference AA
in the ancestral set permits direct candidate changes to derived states; one in
the derived set permits reverse candidate changes to ancestral states. Both can
apply. The opposite phenotype group is never used to infer an ancestor. Missing
ASR states are unresolved; gaps/ambiguous residues are excluded explicitly.
Only nuclear chromosomes 1–22, X and Y with the standard genetic code are
supported; mitochondrial/alternate-contig references require a future explicit
coordinate and translation contract. Only single nucleotide changes that translate to a requested candidate AA are
queried. Multi-nucleotide substitutions, indels and stop variants are outside
this contract, not negative annotation results. A hit describes a human variant;
it is not evidence that a substitution occurred in another species, pathogenicity
in that species, positive selection, or cancer causality.

## Verified reference sidecar

For each uniquely resolved `<gene>[.<suffix>].map.tsv`, supply
`<map filename>.vep.json` alongside it. Missing MAPs and missing sidecars are
explicit per-target statuses, not guessed references. Ambiguous filenames,
duplicate coordinates, malformed metadata and mismatched checksums fail.
Legacy MAPs alone cannot verify the reference protein, codon boundaries or strand.
Sidecars must be derived from the same pinned annotation and reference sequence
used to construct the MAP; this workflow verifies their internal consistency and
checksums, not the biological correctness of an untrusted source declaration.
Example (checksums abbreviated here only):

```json
{
  "schema": "phylophere.vep_reference.v1",
  "gene": "001",
  "species": "homo_sapiens",
  "assembly": "GRCh38",
  "annotation_release": 115,
  "source_uri": "reference producer and pinned annotation/sequence location",
  "map_sha256": "SHA256 of the MAP bytes",
  "gene_id": "ENSG00000100001.1",
  "transcript_id": "ENST00000100001.2",
  "protein_id": "ENSP00000100001.2",
  "protein_sequence": "A",
  "protein_sha256": "SHA256 of the exact protein_sequence text",
  "strand": "+",
  "transcript_aliases": [],
  "codons": {
    "1": {"positions": ["chr1:100", "chr1:101", "chr1:900"], "coding_bases": "GCT"}
  }
}
```

Protein, transcript and gene IDs are versioned. Optional transcript aliases must
refer to the same verified transcript, not other isoforms. The optional existing
gene/protein TSV is read literally; duplicate genes and protein incompatibility
fail. Missing protein IDs exclude that target when the TSV is provided.
The MAP has named columns `status, prot_ali_col, hg38_aa_pos, hg38_nt_coord`.
Selected columns must be unique and contiguous 1..n; removed columns have no
trimmed coordinate. Reference residue numbers are explicit, unique and one-based.
Missing reference residues indicate human alignment gaps. `hg38_nt_coord`, when
present, is the first codon base **in coding order** and must agree with the
sidecar. The sidecar supplies all three genomic one-based forward coordinates,
in coding order (decreasing on minus strand); they may cross an exon boundary.
`coding_bases` are in coding orientation, translate under the standard nuclear
code to the declared protein, and are complemented to make minus-strand genomic
alleles. No strand or contiguous codon is inferred. Missing codons stay unmapped.
The current annotation contract is restricted to human GRCh38.

## Database snapshots

Each supplied PrimateAI or COSMIC file needs `<database path>.vep.json`:

```json
{
  "schema": "phylophere.vep_database.v1",
  "source": "primateai",
  "species": "homo_sapiens",
  "assembly": "GRCh38",
  "coordinates": "genomic_1based_forward",
  "version": "pinned snapshot identifier",
  "annotation_release": 115,
  "source_uri": "original snapshot location and preparation provenance",
  "sha256": "SHA256 of the complete database file",
  "fields": {}
}
```

Files are literal TSV (plain or gzip). `fields` maps canonical keys to physical
column names. PrimateAI defaults: `chrom=chr, pos=pos, ref=non_flipped_ref,
alt=non_flipped_alt, ref_aa=ref_aa, alt_aa=alt_aa,
aa_pos=change_position_1based, transcript_id=gene_name,
gene_id=ensembl_gene_id, score=score_PAI3D, percentile=percentile_PAI3D,
prediction=prediction`. Gene identity must be added during snapshot preparation
if absent in the distributed source, using the same pinned transcript annotation;
never use a gene symbol as an Ensembl ID. Scores and percentiles are fractions
in [0,1]; predictions use the source's gene-specific benign/pathogenic call.
No new pathogenicity cutoff is applied.

COSMIC defaults: `chrom=CHROMOSOME, pos=GENOME_START, ref=GENOMIC_WT_ALLELE,
alt=GENOMIC_MUT_ALLELE, aa_change=MUTATION_AA,
transcript_id=TRANSCRIPT_ACCESSION, gene_id=GENE_ID`.
An adapter may rename columns through `fields` and normalize documented source
allele orientation/coordinate units before producing the checksummed snapshot.
COSMIC additionally requires `scope_gene_ids`: the explicit Ensembl gene IDs of
the snapshot's Census coverage, including genes with no matching mutations.
Its absence cannot be inferred from observed hits. Noncovered genes are
`outside_source_scope`. Only exact genomic alleles and verified transcript,
gene, residue number and AA changes can match. Transcript aliases require
reference-side verification. Versions, when supplied by a source, must match;
unversioned source IDs are accepted only within the matched pinned annotation
release. Identical source rows are deduplicated; distinct source records retain
all columns under `db_` names and a record hash. No isoform is chosen by file order.

## Ensembl calculation and resources

Input is genomic VCF, not protein HGVS. All ready references share one release;
the VEP executable must match it. Query uses explicit transcript versions and
reference checking. Only consequences on the verified reference transcript are
transferred. MANE/canonical flags remain descriptive; they do not authorize
transferring a different isoform. Raw input, VEP output, version and logs are
published separately from the final table. Partial/absent matching transcript
results are not negative consequences. VEP nonzero exit or warnings are technical
failures and terminate; details remain in the work directory.

No resources are installed or downloaded during annotation. Offline mode requires
a prepared `vep_resource.json` in the cache directory with schema
`phylophere.vep_cache.v1`, `species`, `assembly`, numeric `release`, `source_uri`,
`fasta` (relative path), and `files` (every other file relative path -> SHA256,
including cache data, reference FASTA and prebuilt indexes). Expected Ensembl
layout is `homo_sapiens/<release>_GRCh38`. Missing/extra/changed files fail;
use an immutable prepared cache. The FASTA must have its prebuilt `.fai` index.
The resource preparation step is deliberately separate from calculation.

`auto` chooses a prepared offline manifest if present; otherwise database mode.
It never probes connectivity to decide whether to download. Database mode pins
the declared release and records host/port, but cannot checksum a remote MySQL
snapshot, so its Ensembl task is deliberately not reused with `-resume`.
Offline calculation fingerprints all resource content and sources. Public
parameter names remain; this document supersedes older auto-download help until
the separate configuration/help checkpoint is reviewed.

## Outputs, coverage and caching

Public channels and filenames remain `primateai_tsv` (`primateai_mapped.tsv`),
`cosmic_tsv` (`cosmic_scores.tsv`), `ensembl_vep_tsv`
(`ensembl_vep_mapped.tsv`). Every table carries `(Gene, Position, side)`,
`status`, `coverage`, candidate/annotated SNV counts, reference identities and
transfer semantics. Annotated targets have one row per matching variant/source
record; nonannotated targets have one status row. `no_hit` is valid only for a
covered COSMIC gene after searching the complete declared snapshot. It means no
matching record in that snapshot, not a biological negative. Missing PrimateAI
records are `missing_annotation`; a benign summary requires complete candidate
coverage. `unqueryable_alt_aas` records desired states that require multiple
nucleotide changes; mixed reachable/unreachable states retain partial coverage
even when all candidate SNVs have database records. Status counts in manifests count targets, not variant/transcript rows.
Empty SCORING tables still publish valid headers and complete audits.

Prepared targets store each protein/codon reference once in a registry to avoid
quadratic duplication across positions. Input identities and prepared targets live in `vep/vep_inputs/`. Each source
publishes `<output>.manifest.json` with output and input checksums, versions,
coverage counts and metadata. Content fingerprints include MAPs, sidecars,
databases, manifests, source code, and offline cache contents. Source scripts
are copied for execution, but their hashes are explicit cache inputs. Repeated
rendering consumes final tables and does not rerun annotation. Downstream readers
use literal IDs, side-specific keys and explicit status. Unverified legacy tables
are excluded with a warning; regenerate them. No annotation modifies SCORING
scores, positional/gene calibration or BH families. Main/alternative entries,
configuration/help, WorkflowMap and GUI/runners remain separate checkpoints.
