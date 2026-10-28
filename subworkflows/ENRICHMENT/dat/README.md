# Enrichment data

Files read by the enrichment module when the run does not supply its own.

## Gene sets (`*.gmt`)

Default content of `gmt_dir`: every `*.gmt` in this directory is loaded by the FCS enrichment. The source and download date of each file are not recorded here. With `auto_fetch_gmt`, `bin/resolve_gmts.py` adds the current GO (Enrichr) and WikiPathways sets to a copy of them and writes `gmt_source.json` with the origin and SHA-256 of every file; without it the copy holds only these files.
