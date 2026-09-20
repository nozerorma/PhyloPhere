# Vendored eggNOG Primates orthogroups (human subset)

Offline fallback for `--egg_members_file` / `--egg_annotations_file` when left
blank and a live fetch (`bin/resolve_eggnog.py`) fails or isn't attempted.
Fetching is tried first — these vendored copies exist so the pipeline still
runs on an offline cluster node, not to be the primary source of truth.

| File | Source | License | Fetched |
|---|---|---|---|
| `9443_members_human.tsv.gz` | [eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_members.tsv.gz](http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_members.tsv.gz) (Primates, taxid 9443) | CC-BY 4.0 (eggNOG) | 2026-09-20 |
| `9443_annotations_human.tsv.gz` | [eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_annotations.tsv.gz](http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_annotations.tsv.gz) | CC-BY 4.0 (eggNOG) | 2026-09-20 |

`subworkflows/ENRICHMENT/local/src/build_position_gmt.py` only ever reads the
orthogroup id, its description, and members whose taxon prefix is `9606.ENSP`
(human) — every other member and annotation column is discarded on load. So
these vendored copies are filtered down to that subset at the row level: of
the 23,677 Primates-level orthogroups in the full upstream file, the 18,408
that contain at least one human member are kept, non-human members are
dropped from each row's member list, and the annotations file is filtered to
the same 18,408 orthogroup ids. This shrinks the pair from ~1.8MB to ~470KB
compressed without changing what the pipeline actually uses.

To refresh: re-run the fetch URLs above and re-apply the same filter (see
`bin/resolve_eggnog.py`), then replace the files and update the table above.
