# Run from the repository root with the project's Rscript.
repo <- normalizePath(".")
src <- file.path(repo, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "percentile_flags.R"))
source(file.path(src, "fcs_enrich.R"))

same <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE)))

# The two rules are independent. The analytic rule needs raw p.perm below its
# threshold and the smallest eligible analytic level; sig_perm is BH on p.perm.
d <- data.frame(
  ranking = "global", database = "db", pathway = paste0("S", 1:8),
  stat = c(.2, .2, .2, .2, .2, .2, -.2, -.2),
  p.adj = c(.03, .08, .12, .20, .03, .04, .03, .03),
  p.perm = c(.01, .01, .02, .01, .30, .025, .01, NA_real_),
  sig_top = TRUE # Old input is accepted, but the removed column is not propagated.
)
r <- fcs_classify_significance(d, fdr_thr = .15, p_perm_thr = .025)
stopifnot(!"sig_top" %in% names(r),
          identical(r$redlich_level, c(.05, .10, .15, NA, NA, NA, .05, NA)),
          identical(r$sig_perm, c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)),
          identical(r$direction, c("positive", "positive", "positive", "positive", NA, "positive", "negative", NA)),
          is.na(r$q.perm[8]))

# User-supplied levels work independently and are sorted to report the lowest
# qualifying level even when their input order differs.
custom <- fcs_classify_significance(d, .15, .025, analytic_levels = c(.15, .05, .10))
same(custom$redlich_level, r$redlich_level)
custom_one <- fcs_classify_significance(d, .15, .025, analytic_levels = .10)
stopifnot(identical(custom_one$redlich_level, c(.10, .10, NA, NA, NA, NA, .10, NA)))

# Empty output has stable column names/types and never carries the legacy flag.
empty <- fcs_run_all(list(global = setNames(numeric(), character())),
                     list(db = list(genesets = list(), geneset.names = character())),
                     perms_file = "NO_FILE")
expected <- c("ranking", "database", "pathway", "stat", "pval", "p.adj", "p.perm", "q.perm",
              "num.genes", "gene.vals", "sig_perm", "direction", "redlich_level")
stopifnot(nrow(empty) == 0L, identical(names(empty), expected), !"sig_top" %in% names(empty))
types <- fcs_enrich_col_types()$cols
stopifnot("collector_logical" %in% class(types$sig_perm),
          "collector_double" %in% class(types$redlich_level))

cat("PASS: independent significance rules, custom levels, legacy sig_top removal, signed stats, missing null, empty schema\n")
