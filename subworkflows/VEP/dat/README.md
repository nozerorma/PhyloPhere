# VEP dat/ — Static Reference Data

This directory holds the large static reference files required by the VEP
characterization workflow.  Files are **not** bundled with the repository
(sizes range from hundreds of MB to ~2 GB).  Place or symlink each file here
before running `--vep`, or override the corresponding parameter on the command
line to point to files stored elsewhere.

## Required files

| File | Source | Default param |
|------|--------|---------------|
| `Homo_sapiens.cds.fa.gz` | Ensembl CDS FASTA for Homo sapiens | `--vep_hs_cds` |
| `Homo_sapiens.sorted.gff` | Ensembl GFF annotation (sorted by transcript ID) | `--vep_gff` |
| `proteiID_gene_equivalences.txt` | Two-column TSV: gene_name → Ensembl protein ID | `--vep_gene_equiv` |
| `transvar/` | TransVar index directory (hg38.ensembl.gtf.gz, hg38.ccds.txt, etc.) | `--vep_transvar_db` |
| `PrimateAI-3D.hg38.txt.gz` | PrimateAI-3D pathogenicity scores (hg38) | `--vep_primateai_db` |

## Origin of files

These files correspond directly to the contents of:

- `to_integrate/prot2pos/needed_files/`  → `Homo_sapiens.cds.fa.gz`, `Homo_sapiens.sorted.gff`,
  `proteiID_gene_equivalences.txt`
- `to_integrate/prot2pos/fa_transvar_ref/` → `transvar/`
- `to_integrate/prot2pos/PrimateAI/`       → `PrimateAI-3D.hg38.txt.gz`

## Per-project files (NOT stored here)

The per-gene CDS FASTA files (`GENE.Homo*.fasta`) and codon-filtering HTML
tracking files (`GENE*.html`) are project-specific outputs from the upstream
codon-filtering step.  Supply them via:

```
--vep_cds_dir   /path/to/CDS/
--vep_track_dir /path/to/TRACK/
```

## TransVar database structure

The `transvar/` directory must contain the index files built by TransVar for
hg38.  Run `transvar config --download_anno --refversion hg38` to populate it,
or copy the files from `fa_transvar_ref/`. The workflow writes a temporary
`transvar.cfg` and points `TRANSVAR_CFG` at it so TransVar resolves the staged
CCDS `.transvardb` files through its supported configuration mechanism.
