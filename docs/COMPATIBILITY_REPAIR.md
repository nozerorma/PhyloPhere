# Compatibility repair: existing references and execution cost

This repair separates scientific compatibility from artifact integrity. The
former remains required for inference; ordinary execution no longer requires
private snapshots that existing producers never emitted, or repeated whole-data
hashes. This is not a full scientific re-audit of every workflow.

## Findings and decisions

| Area | Existing-input problem | Repair and retained check |
| --- | --- | --- |
| VEP | Mandatory MAP/database sidecars; native COSMIC IDs treated as missing Ensembl IDs | Read native PrimateAI/COSMIC schemas and declared MAP coordinates. Check location, protein residue and compatible AA states; check supplied isoform identities. Missing matches are unavailable annotation. |
| POSENRICH | VEP reference requirement also removed declared protein/genomic coordinates from positional annotations | Preserve explicit MAP mappings. Accept native FUBAR sites and existing projected Pfam/UCR tables. Keep scanned coverage separate from positive memberships. |
| CT_POSTPROC | Missing completion record prevented even descriptive processing | Validate legacy tables; mark scientific provenance unknown. Unknown settings cannot support calibrated inference. |
| ASR_DIAGNOSTICS | Missing completion record prevented legacy diagnostics | Validate table structure and preserve unknown model/threshold provenance. |
| SCORING | Additional artifact manifest and repeated MD5/SHA reads were prerequisites | Optional integrity checking; retain observed/null scientific settings, filters, universe, cycle and matrix consistency checks. |
| Task caching | 46 deep-cache directives across 10 workflow families, plus controller content scans | Ordinary Nextflow caching; data directory fingerprints use names, sizes and modification times. Small implementation sources remain content-sensitive. |

The scientific coordinate concern was justified: the legacy MAP producer exports
only a codon's first base. A split codon with CDS coordinates `[100,101,900]` must
not become `[100,101,102]`. Requiring nonexistent sidecars for all annotation was
not a workable repair. Legacy anchors and protein positions remain usable;
unresolved bases limit coverage explicitly. The local Ortholog Characterizator
exporter now appends all three ordered coordinates, strand, protein ID and coding
triplet from inputs it already reads. Its existing cluster exports are unchanged.

`A1BG.Homo_sapiens.prot` is a source record, not a special gene in pipeline logic.
The compatibility adapter explicitly recognizes `.prot` and `.Homo_sapiens`
producer suffixes for any gene. Exact IDs win; ambiguous alternatives fail;
arbitrary dotted IDs remain literal. This is an explicit naming convention, not
a universal identifier normalization rule.

## Validation and limits

Targeted checks cover native database headers/records, real MAP and projected
FUBAR/Pfam/UCR/PP2A samples, literal IDs, ambiguous aliases, split/minus-strand
codons, schema/residue mismatches, legacy diagnostics, scientific null rejection,
extended MAP export, and local Nextflow execution/resume. A small complete
POSENRICH run includes native FUBAR, custom PP2A, null calibration and HTML output.
The final serial batch passed all 33 checks in 50.24 seconds; `git diff --check`
also passed.

The serial local regression script is
`validation/test_compatibility_repair.py` (the repository ignores `validation/`).
Run with `conda run --no-capture-output -n phylophere python -B -m pytest
validation/test_compatibility_repair.py -q -p no:cacheprovider`. Real external
samples are optional via `PHYLOPHERE_REPAIR_SAMPLES`; they are not distributed in
the repository. Local verification logs are `/tmp/phylophere-repair.log` and
`/tmp/phylophere-repair.xml`.

Ensembl tests use an executable fixture and mocked reference responses: no real
installed VEP or full scientific dataset has been validated here. Exact SNV
enumeration currently supports nuclear standard-code references. Legacy MAP-only
annotation can have partial coverage. Old nulls with unknown scientific settings
remain unsuitable for calibration. Metadata caching does not detect edits that
preserve both size and mtime; version resource paths or explicitly rerun those
tasks. `PHYLOPHERE_VERIFY_CHECKSUMS=1` enables supported optional integrity checks.

The broader decision audit of SCORING, RER, FADE, SELECTION_PREP, CT_POSTPROC,
ASR_DIAGNOSTICS, CT_META_CAAS and CONTRAST_SELECTION continues separately.
