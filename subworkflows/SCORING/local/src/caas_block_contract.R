# Validate one CAAS block at a time; retain only IDs and metadata across blocks.
caas_ids <- function(x, what) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) stop(paste("Invalid", what, "IDs"))
  x
}
caas_block_index <- function(files, roster = NULL, universe = character(), gene_stat = NULL) {
  if (!length(files) || any(!file.exists(files))) stop("Missing CAAS matrix blocks")
  seen <- character(); genes <- caas_ids(universe, "universe"); entries <- list(); units <- 0
  for (f in files) {
    b <- readRDS(f)
    g <- caas_ids(b$genes, "gene"); cy <- caas_ids(b$cycles, "cycle")
    if (!length(cy) || any(!grepl("^b_[1-9][0-9]*(~H[0-9]+)?$", cy))) stop("Invalid null cycle IDs")
    if (any(cy %in% seen)) stop("cycles occur more than once within or across blocks")
    if (is.null(gene_stat)) gene_stat <- b$gene_stat
    if (length(gene_stat) != 1 || !identical(gene_stat, "max") ||
        !identical(b$gene_stat, gene_stat)) stop("CAAS gene_stat mismatch: only max is supported; rebuild legacy gene matrices")
    check_matrix <- function(m) {
      if (!is.matrix(m) || !is.numeric(m) || !identical(rownames(m), if(length(g)) g else NULL) ||
          !identical(colnames(m), cy) || any(!is.finite(m))) stop(paste("Invalid matrix dimensions, IDs, order or values in", f))
    }
    dirs <- c("global", "top", "bottom")
    if (!identical(names(b$stat), dirs)) stop("Missing CAAS rankings")
    for (d in dirs) {
      check_matrix(b$stat[[d]])
      if(any(b$stat[[d]]<0 | b$stat[[d]]>1)) stop('Invalid CAAS score domain')
    }
    if(!isTRUE(all.equal(as.numeric(b$stat$global),pmax(as.numeric(b$stat$top),as.numeric(b$stat$bottom)),tolerance=1e-12))) stop('CAAS global score differs from best direction')
    if (!is.null(b$parts)) {
      if (!all(c("max", "n") %in% names(b$parts))) stop("Missing CAAS parts")
      for (part in b$parts) for (d in dirs) check_matrix(part[[d]])
      for (d in dirs) {
        n <- b$parts$n[[d]]
        if (any(n < 0 | n != floor(n))) stop("Invalid CAAS parts")
        expected <- b$parts$max[[d]]
        if (!isTRUE(all.equal(as.numeric(b$stat[[d]]), as.numeric(expected), tolerance=1e-12))) stop("CAAS stat disagrees with parts/gene_stat")
      }
    }
    if (length(b$n_units) != 1 || !is.finite(b$n_units) || b$n_units < 0 || b$n_units != floor(b$n_units)) stop("Invalid n_units")
    local_roster <- file.path(dirname(f), "cycles.txt")
    if (file.exists(local_roster) && !identical(readLines(local_roster), cy)) stop("Block roster differs from matrix cycles")
    seen <- c(seen, cy); genes <- union(genes, g); units <- units + b$n_units
    entries[[length(entries)+1L]] <- list(md5=unname(tools::md5sum(f)), cycles=cy, genes=g)
    rm(b)
  }
  if (!is.null(roster)) {
    caas_ids(roster, "roster")
    if (!setequal(roster, seen)) stop("CAAS block cycles not in roster or roster cycles missing from blocks")
  }
  list(format="caas_blocks_v1", gene_stat=gene_stat, genes=sort(genes), cycles=sort(seen),
       n_units=units, blocks=entries)
}
caas_export_single <- function(files, index, output) {
  dirs <- c("global", "top", "bottom")
  mats <- setNames(lapply(dirs, function(d) matrix(0, length(index$genes), length(index$cycles),
                     dimnames=list(index$genes, index$cycles))), dirs)
  for (f in files) {
    b <- readRDS(f)
    for (d in dirs) mats[[d]][b$genes, b$cycles] <- b$stat[[d]]
  }
  # Preserve the historical empty-null export schema; the canonical roster retains N.
  if (index$n_units == 0) { message("[caas_perms] no scored rows; writing empty RDS"); mats <- list() }
  saveRDS(list(caas_corStat_byrank=mats, gene_stat=index$gene_stat), output)
}

# Whole-memory control scripts opt in to reconstruction; productive FCS uses blocks.
caas_control_null <- function(store) {
  files <- Sys.glob(file.path(store, "blocks/blk_*/matrices.rds"))
  if (dir.exists(file.path(store, "blocks"))) {
    roster <- file.path(store, "cycle_roster.txt")
    if (!file.exists(roster)) stop("Missing canonical CAAS cycle roster")
    index <- caas_block_index(files, readLines(roster))
    output <- tempfile("caas_export_", fileext=".rds")
    on.exit(unlink(output))
    caas_export_single(files, index, output)
    return(readRDS(output))
  }
  readRDS(file.path(store, "caas_perms.rds"))
}
