# Vendored GMT gene sets

Offline fallback for `--gmt_dir` when left blank and a live fetch
(`bin/resolve_gmts.py`) fails or isn't attempted. Fetching is tried first —
these vendored copies exist so the pipeline still runs on an offline
cluster node, not to be the primary source of truth.

| File | Source | License | Fetched |
|---|---|---|---|
| `go_biological_process.gmt` | [Enrichr `GO_Biological_Process_2023`](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023) (derived from the Gene Ontology) | CC-BY 4.0 (Gene Ontology) | 2026-09-19 |
| `go_molecular_function.gmt` | [Enrichr `GO_Molecular_Function_2023`](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2023) | CC-BY 4.0 (Gene Ontology) | 2026-09-19 |
| `reactome_pathways.gmt` | [reactome.org/download/current/ReactomePathways.gmt.zip](https://reactome.org/download/current/ReactomePathways.gmt.zip) (human only, `R-HSA-*`) | CC0 | 2026-09-19 |
| `wikipathways.gmt` | [data.wikipathways.org/current/gmt](https://data.wikipathways.org/current/gmt/) (`Homo_sapiens`) | CC0 | 2026-09-19 |

All four use HGNC gene symbols directly (no ID conversion needed) —
`load_gmts()` (subworkflows/ENRICHMENT/local/src/posenrich_enrich.py) and
the FCS engine match genes by symbol.

**Not vendored: KEGG** — commercial-use-restricted, excluded per project
policy (see the external-file-dependency reduction plan).

To refresh: re-run the fetch URLs above (dates/versions change upstream
periodically) and replace the files, updating the table above.
