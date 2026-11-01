# VEP inputs and interpretation

VEP consumes exactly one SCORING `position_scores.tsv`. Gene IDs are literal TSV
fields (including `001`, `NA`, internal dots and quotes). The target key is
`(Gene, Position, side)`; positions are zero-based productive alignment columns.
Scores must be finite and positive, side must be top/bottom, and duplicate or
malformed keys fail before annotation. Public workflow inputs, output channel
names and output filenames remain unchanged.

## Existing inputs are supported

The original MAP, PrimateAI-3D and COSMIC files are accepted directly. Private
`*.vep.json` reference/database files and `vep_resource.json` cache manifests are
**not required**. Old optional reference sidecars can still be used when present;
contradictions in explicitly supplied metadata are errors. Supplied databases
are checked for existence and nonempty size before task dispatch, not after
hashing the MAP collection. Database schema checks read their headers first.

MAP has named columns `status, prot_ali_col, hg38_aa_pos, hg38_nt_coord` (POSENRICH
also reads `ori_codon_col`). `Position + 1` indexes `prot_ali_col`. Removed columns
have no productive index; gaps, missing coordinates, out-of-range positions and
duplicates have explicit handling. Literal MAP filenames are resolved once per
process, preserving dots; an exact name takes priority and ambiguous alternatives
fail. Gene/protein tables are optional for PrimateAI/COSMIC and needed for legacy
Ensembl protein resolution. Productive human alignment residues are checked when
`--alignment` is available.

MAP `hg38_nt_coord` denotes the first nucleotide **in coding order**, at one-based
forward-genome coordinates. It does not establish the other two bases. A legacy
codon is reconstructed only when a consistent coordinate order and the next
reference residue's anchor exactly three bases away establish contiguity. The
next residue can be a removed alignment column because MAP contains original
columns too. Otherwise only the declared anchor is searched, with partial
coverage. A gene's final codon or a codon crossing an exon boundary is never
completed by guessing `+1/+2` or `-1/-2`.

The paired Ortholog Characterizator producer now appends `hg38_nt_coord_2`,
`hg38_nt_coord_3`, `strand`, `human_protein_id` and `hg38_codon` to its MAP output.
The three coordinates come from the existing ordered CDS list, so split codons
are preserved; coding bases come from the reference CDS alignment. Older columns
remain readable. Regenerating this export requires its original TRACK/CDS/GFF
inputs, not rerunning BMGE or fitting a model. Existing cluster files are not
modified by local code changes.

## Original annotation sources

PrimateAI uses its actual fields: `chr, pos, non_flipped_ref, non_flipped_alt,
gene_name` (Ensembl transcript), `change_position_1based, ref_aa, alt_aa,
score_PAI3D`, with optional `percentile_PAI3D, refseq, prediction`. No added ENSG
column or exact release equality is demanded. Scores/predictions are preserved,
not recalibrated. Candidate joins require MAP genomic location and protein residue
number, compatible AA states and the productive reference AA when available.
Declared transcript identities and nucleotide alleles of a fully resolved codon
are also checked. If the MAP does not identify the isoform, source transcript IDs
remain source annotations, not a claim of independently verified isoform identity.

COSMIC uses `CHROMOSOME, GENOME_START, GENOMIC_WT_ALLELE,
GENOMIC_MUT_ALLELE, MUTATION_AA, TRANSCRIPT_ACCESSION`; original `GENE_SYMBOL`,
`COSMIC_GENE_ID`, `HGVSP` and other fields are retained. `COSMIC_GENE_ID` is never
interpreted as ENSG. The numeric residue in `p.X123Y` must match the MAP protein
position; `HGVSP` protein identity is checked when a protein ID was supplied.
Unlocated source rows and non-SNV/non-missense records are outside this join.
Exact duplicate source records are deduplicated; distinct records remain separate
annotation rows, retaining the complete target key.

A missing match is `missing_annotation`, not benign/not-cancer-relevant evidence.
A Census scope is not needed to report a positive coincidence. Missing source
coverage is never inferred from positive records. Both files are streamed; if no
target has usable coordinates the schema is checked without scanning the file.
No new snapshot preparation step is required. A sibling database manifest can be
provided to the script for optional source metadata/field aliases, but ordinary
workflow use does not require one.

## Ensembl and reference resolution

Ensembl VEP consumes genomic VCF, not protein HGVS. If legacy targets lack exact
codon bases, database mode resolves the supplied human protein through Ensembl
REST: protein/CDS sequence, parent transcript and translation-to-genome mapping.
It checks the productive reference AA and MAP anchor before generating SNVs and
keeps the responses in the task audit. Unresolvable proteins and conflicting
reference coordinates remain explicit unavailable targets. Service failures are
technical failures, not negative annotations. This lookup happens only when
Ensembl annotation is enabled, not for the PrimateAI/COSMIC workflows.

SNV enumeration supports the nuclear standard genetic code. Mitochondrial MAP
coordinates can be joined to source records, but are not translated with that
code; unresolved Ensembl mitochondrial codons remain explicitly unavailable.

Offline mode uses the standard `homo_sapiens/<installed-release>_GRCh38` cache and
a preindexed genomic FASTA. No private checksum manifest is needed. Offline targets
need already resolved codon bases/coordinates, for example extended MAP columns.
Missing codons are reported explicitly. `auto` prefers an existing standard cache,
otherwise database mode. Resources are prepared separately; annotation does not
install or download a multi-GB cache. VEP nonzero exit and errors fail the task;
variant warnings are published and unreturned compatible consequences remain
unresolved. Only compatible protein/transcript consequences are transferred.

## Meaning, coverage and caching

SCORING's ancestral/derived columns pool states across domains; they do not retain
paired ancestral substitutions. All matches are labelled human reference candidate
transfers from pooled ASR states. An ancestor is never inferred from the opposite
phenotype group. A human reference state can support direct or reverse candidate
changes; neither establishes an event observed in another species, pathogenicity
there, positive selection or cancer causality.

Tables preserve `Gene, Position, side`, status, mapping basis, reference status
and coverage. `candidate_count_known` distinguishes resolved SNV candidate sets
from a lower bound based on source matches. MAP-only results retain partial
coverage; a benign summary requires a complete known candidate set. Positions
without a reconstructed ancestor, with gaps, unknown coordinates or ambiguous
residues are explicit exclusions/unavailable annotations. A partial match can be
useful descriptive evidence without certifying all possible substitutions.

Normal runs use Nextflow's ordinary name/size/mtime cache policy. Directory inputs
have a metadata inventory; their data bytes are not hashed by the controller.
Small source scripts retain content fingerprints so code edits invalidate their
tasks. Database mode avoids reusing a mutable remote Ensembl query. Each source
writes a lightweight audit of paths, mapping limitations and statuses, separate
from raw Ensembl artifacts. No full database checksum is computed on every run.
The usual metadata-cache limitation applies to in-place edits preserving size and
mtime: use a fresh resource version/path or explicitly rerun the affected stage.
