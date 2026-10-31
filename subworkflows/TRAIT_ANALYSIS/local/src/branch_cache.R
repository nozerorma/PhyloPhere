# branch_cache.R — Verify branch-color caches against their tree and source trait.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/

branch_cache_provenance <- function(tree, species, values, trait_name) {
  data <- data.frame(species = as.character(species),
                     value = numeric_trait_values(values, trait_name))
  data <- data[data$species %in% tree$tip.label & !is.na(data$value), , drop = FALSE]
  data <- data[order(data$species, data$value), , drop = FALSE]
  rownames(data) <- NULL
  list(version = 1L, trait = trait_name, data = data,
       tips = tree$tip.label, edge = unname(tree$edge),
       lengths = tree$edge.length, root_edge = tree$root.edge)
}

read_branch_cache <- function(path, provenance) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path) || !file.exists(paste0(path, ".branch_cache.rds"))) {
    message("Recomputing branch colors; cache or provenance sidecar is missing.")
    return(NULL)
  }
  tryCatch({
    saved <- readRDS(paste0(path, ".branch_cache.rds"))
    if (!identical(saved, provenance)) stop("tree or trait provenance differs")
    cached <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("node", "BR", "species") %in% names(cached))) stop("missing cache columns")
    n <- length(provenance$tips)
    nodes <- cached$node
    if (!is.numeric(nodes) || anyNA(nodes) || any(nodes != as.integer(nodes)) ||
        any(nodes < 1 | nodes > max(provenance$edge)) || anyDuplicated(nodes)) {
      stop("invalid cache node coordinates")
    }
    if (!identical(as.character(cached$species[match(seq_len(n), nodes)]), provenance$tips)) {
      stop("cache tip coordinates differ")
    }
    cached$BR <- numeric_trait_values(cached$BR, "cached BR")
    cached[, c("node", "BR")]
  }, error = function(e) {
    message("Recomputing branch colors; cache unavailable or unverifiable: ", conditionMessage(e))
    NULL
  })
}
