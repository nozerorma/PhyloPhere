#!/usr/bin/env Rscript
# =============================================================================
# Download GMT gene-set files for the ORA engine
# =============================================================================
# The ORA engine (ORA_general.Rmd) runs
# clusterProfiler against EVERY *.gmt file found in `gmt_dest_dir`. WebGestalt is
# no longer used as a runtime backend -- it is used here only as a *downloader*
# to materialise the curated WebGestalt databases as symbol-based GMT files.
#
# This script produces two families of GMTs in `gmt_dest_dir`:
#   1. WebGestalt curated databases  -> "<dbname>.symbols.gmt"
#      (genes mapped entrezgene -> HGNC symbol so they join the symbol-based
#       gene lists / backgrounds the pipeline uses).
#   2. MSigDB collections (latest release) -> "<coll>.all.v<rel>.symbols.gmt".
#
# All gene identifiers are HGNC symbols, matching the pipeline's gene lists.
# =============================================================================

suppressMessages({
  library(WebGestaltR)
})

# ── Destinations ─────────────────────────────────────────────────────────────
script_args <- commandArgs(trailingOnly = TRUE)
gmt_dest_dir <- if (length(script_args) >= 1 && nzchar(script_args[1])) {
  script_args[1]
} else {
  "/home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/ORA/dat"
}
cache_dir <- file.path(tempdir(), "wg_cache")

if (!dir.exists(gmt_dest_dir)) dir.create(gmt_dest_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(cache_dir))    dir.create(cache_dir,    recursive = TRUE, showWarnings = FALSE)

organism <- "hsapiens"

# ── 1. WebGestalt curated databases ──────────────────────────────────────────
# CURATED, non-redundant subset (current WebGestaltR 0.4.6 identifiers). The FCS
# engine BH-corrects per GMT, so we keep only complementary, non-overlapping
# ontologies to avoid inflating the number of per-database test families and to
# avoid double-counting GO/pathways that also live in MSigDB.
webgestalt_dbs <- c(
  # Gene Ontology (non-redundant; MSigDB C5 is intentionally dropped as duplicate)
  "geneontology_Biological_Process_noRedundant",
  "geneontology_Cellular_Component_noRedundant",
  "geneontology_Molecular_Function_noRedundant",
  # Pathways (MSigDB C2 intentionally dropped as duplicate pathway source)
  "pathway_KEGG",
  "pathway_Reactome",
  "pathway_Wikipathway",          # was pathway_WikiPathways (renamed upstream)
  "pathway_Wikipathway_cancer",   # cancer-specific (renamed from WikiPathways_cancer)
  # Disease / phenotype (cancer-relevant)
  "disease_Disgenet",
  "disease_OMIM",
  "phenotype_Human_Phenotype_Ontology",
  # Networks: protein complexes + regulatory targets
  "network_CORUMA",                       # CORUM all (NOT network_CORUM)
  "network_Transcription_Factor_target",
  "network_Kinase_target"
  # Dropped from the previous "extra" set as noisy/redundant for this study:
  #   pathway_Panther, network_PPI_BIOGRID, network_miRNA_target,
  #   chromosomalLocation_CytogeneticBand
)

# Write a data.frame (geneSet, description, gene...) to GMT format.
write_gmt <- function(df, path) {
  con <- file(path, open = "wt")
  on.exit(close(con))
  for (gs in unique(df$geneSet)) {
    sub <- df[df$geneSet == gs, ]
    desc <- sub$gmt_description[1]
    if (is.na(desc) || !nzchar(desc)) desc <- gs
    genes <- unique(sub$symbol[!is.na(sub$symbol) & nzchar(sub$symbol)])
    if (length(genes) == 0) next
    writeLines(paste(c(gs, desc, genes), collapse = "\t"), con)
  }
}

download_webgestalt_db <- function(db) {
  out_path <- file.path(gmt_dest_dir, sprintf("%s.symbols.gmt", db))
  message(sprintf("→ WebGestalt: %s", db))
  gs <- tryCatch(
    loadGeneSet(organism = organism, enrichDatabase = db, cache = cache_dir),
    error = function(e) { message(sprintf("  ✗ loadGeneSet failed: %s", e$message)); NULL }
  )
  if (is.null(gs) || is.null(gs$geneSet) || nrow(gs$geneSet) == 0) {
    message("  ✗ no gene sets returned, skipping"); return(FALSE)
  }
  mapping <- gs$geneSet
  # entrezgene -> HGNC symbol
  ent <- unique(as.character(mapping$gene))
  idm <- tryCatch(
    idMapping(organism = organism, dataType = "list", inputGene = ent,
              sourceIdType = "entrezgene", targetIdType = "genesymbol",
              cache = cache_dir),
    error = function(e) { message(sprintf("  ✗ idMapping failed: %s", e$message)); NULL }
  )
  if (is.null(idm) || is.null(idm$mapped) || nrow(idm$mapped) == 0) {
    message("  ✗ id mapping empty, skipping"); return(FALSE)
  }
  sym_lookup <- setNames(idm$mapped$geneSymbol, as.character(idm$mapped$userId))
  mapping$symbol <- sym_lookup[as.character(mapping$gene)]
  # Human-readable set descriptions (fall back to the URL/id otherwise)
  if (!is.null(gs$geneSetDes) && nrow(gs$geneSetDes) > 0) {
    des_lookup <- setNames(gs$geneSetDes$description, gs$geneSetDes$geneSet)
    mapping$gmt_description <- des_lookup[mapping$geneSet]
  } else {
    mapping$gmt_description <- mapping$geneSet
  }
  write_gmt(mapping, out_path)
  n_sets <- length(unique(mapping$geneSet))
  message(sprintf("  ✓ wrote %d gene sets -> %s", n_sets, basename(out_path)))
  TRUE
}

# ── 2. MSigDB collections (latest human release) ─────────────────────────────
# Curated MSigDB collections: H (hallmark, compact high-level), C6 (oncogenic
# signatures) and C8 (cell-type signatures) — cancer-relevant. C1 (positional),
# C2 (pathways, duplicate of WebGestalt), C3 (motif), C4 (computational),
# C5 (GO/HPO, duplicate of WebGestalt) and C7 (immunologic) are intentionally
# excluded to keep the per-GMT test families non-redundant.
msigdb_release <- "2026.1.Hs"
msigdb_collections <- c("h", "c6", "c8")
msigdb_base_url <- sprintf(
  "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/%s", msigdb_release)

download_with_retry <- function(url, destfile, tries = 3) {
  for (i in seq_len(tries)) {
    message(sprintf("  Downloading %s (attempt %d/%d)...", basename(destfile), i, tries))
    ok <- tryCatch({
      download.file(url, destfile = destfile, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) { message(sprintf("    error: %s", e$message)); FALSE })
    if (ok && file.exists(destfile) && file.size(destfile) > 1000) return(TRUE)
    Sys.sleep(3)
  }
  FALSE
}

download_msigdb <- function() {
  # Remove MSigDB GMTs that are either from an older release OR from a collection
  # no longer in the curated set, so the engine does not double-count gene sets.
  existing <- list.files(gmt_dest_dir, pattern = "\\.all\\.v.*\\.Hs\\.symbols\\.gmt$",
                         full.names = TRUE)
  wanted <- file.path(gmt_dest_dir,
                      sprintf("%s.all.v%s.symbols.gmt", msigdb_collections, msigdb_release))
  stale <- setdiff(existing, wanted)
  for (s in stale) { message(sprintf("→ removing stale/uncurated MSigDB file %s", basename(s))); file.remove(s) }

  for (coll in msigdb_collections) {
    fname <- sprintf("%s.all.v%s.symbols.gmt", coll, msigdb_release)
    dest  <- file.path(gmt_dest_dir, fname)
    if (file.exists(dest) && file.size(dest) > 1000) {
      message(sprintf("✓ %s already present", fname)); next
    }
    message(sprintf("→ MSigDB: %s", fname))
    ok <- download_with_retry(sprintf("%s/%s", msigdb_base_url, fname), dest)
    if (ok) message(sprintf("  ✓ %s", fname)) else warning(sprintf("✗ failed: %s", fname))
  }
}

# ── Run ──────────────────────────────────────────────────────────────────────
message("=== Downloading WebGestalt curated databases ===")
for (db in webgestalt_dbs) download_webgestalt_db(db)

message("\n=== Downloading MSigDB ", msigdb_release, " collections ===")
download_msigdb()

message("\nGMT download task completed. Destination: ", gmt_dest_dir)
