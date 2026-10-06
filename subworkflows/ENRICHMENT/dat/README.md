# Enrichment data

Files read by the enrichment module when the run does not supply its own.

## Gene sets (`*.gmt`)

Default content of `gmt_dir`: every `*.gmt` in this directory is loaded by the FCS enrichment. The source and download date of each file are not recorded here. With `auto_fetch_gmt`, `bin/resolve_gmts.py` adds the current GO (Enrichr) and WikiPathways sets to a copy of them and writes `gmt_source.json` with the origin and SHA-256 of every file; without it the copy holds only these files.

## eggNOG orthogroups (POSENRICH)

| File | Content |
|---|---|
| `9443_members_human.tsv.gz` | eggNOG 5.0 Primates-level (taxid 9443) orthogroup members, human (`9606.ENSP`) members only |
| `9443_annotations_human.tsv.gz` | annotations of the same orthogroups |

Source: `http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_{members,annotations}.tsv.gz`, downloaded 2026-09-20, license CC-BY 4.0 (eggNOG). `build_position_gmt.py` reads only the orthogroup id, its description and the members whose taxon prefix is `9606.ENSP`, so the files keep the 18,408 of the 23,677 Primates orthogroups with a human member, with only the human members in each row, and the annotations of those orthogroups. This reduces the pair from about 1.8 MB to about 470 KB compressed and does not change what the pipeline reads.

`bin/resolve_eggnog.py` copies this pair when `egg_members_file` and `egg_annotations_file` are blank, for tax level 9443 and reference species 9606. With `auto_fetch_eggnog` it downloads the pair of another clade and applies the same filter. Each run writes `eggnog_source.json` with the origin and the SHA-256 of the files used.

To refresh the versioned pair, run `bin/resolve_eggnog.py --fetch --output-dir <dir>` with the default tax level and reference species, replace the two files here with the outputs (renamed `9443_members_human.tsv.gz` and `9443_annotations_human.tsv.gz`) and update the download date above.
