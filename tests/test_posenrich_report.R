# Render of the POSENRICH report (13.POSENRICH_report.Rmd) on a minimal pipeline fixture, with every component,
# without level 1 and level 2, and with the held-out-cycle validation table.
# Usage: Rscript tests/test_posenrich_report.R <repo> <scratch>
a <- commandArgs(TRUE)
root <- normalizePath(a[1]); scratch <- a[2]
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
src <- file.path(root, 'subworkflows/POSENRICH/local/src')
py <- Sys.getenv('PYTHON'); if (!nzchar(py)) py <- 'python3'
run_py <- function(script, args) {
  rc <- system2(py, c(file.path(src, script), args))
  if (rc != 0) stop(paste('failed:', script))
}

# 1. Build component fixtures under the explicit contracts.
rc <- system2(py, c(file.path(root,'tests/posenrich_report_fixture.py'), root, scratch))
if (rc != 0) stop('POSENRICH component fixture failed')
universe_table <- file.path(scratch,'universe/universe.tsv.gz')

# 2. Render. The Rmd is rendered in place, as the process does in its work directory.
suppressPackageStartupMessages(library(rmarkdown))
rmd <- file.path(scratch,'13.POSENRICH_report.Rmd')
file.copy(file.path(root,'subworkflows/POSENRICH/local/13.POSENRICH_report.Rmd'),rmd,overwrite=TRUE)
dir.create(file.path(scratch,'src'),showWarnings=FALSE)
file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src/report_presentation.R'),file.path(scratch,'src/report_presentation.R'),overwrite=TRUE)
render <- function(name, ...) {
  args <- modifyList(list(universe_dir = file.path(scratch, 'universe'), registry_dir = file.path(scratch, 'registry'),
                          characterization_dir = file.path(scratch, 'characterization'),
                          accumulation_dir = file.path(scratch, 'accumulation'),
                          within_gene_dir = file.path(scratch, 'within_gene'),
                          proximity_dir = file.path(scratch, 'proximity'), traitname = 'trait'), list(...))
  out <- file.path(scratch, name)
  env <- new.env(parent = globalenv())
  rmarkdown::render(rmd, params = args, output_file = out, envir = env, quiet = TRUE)
  stopifnot(file.exists(out))
  paste(readLines(out, warn = FALSE), collapse = '\n')
}
has <- function(html, ...) for (needle in c(...)) if (!grepl(needle, html, fixed = TRUE)) stop(sprintf('missing in report: %s', needle))
lacks <- function(html, ...) for (needle in c(...)) if (grepl(needle, html, fixed = TRUE)) stop(sprintf('unexpected in report: %s', needle))

# 2a. Every component
html <- render('full.html')
has(html, 'POSENRICH: position-level enrichment', 'Universe and coverage', 'Level 1: accumulation per gene',
    'Level 2: concentration within genes', 'Overlaps and proximity (descriptive)', 'Methods and limits', 'TERM_A',
    'Gene-list check', 'Genes above their share of U', 'Overlaps (exact identity)', 'Near positions')
lacks(html, 'mathjax.rstudio.com', 'not available')
# Numbers come from the tables: |U| = 12 positions, 12 observed positives
u <- read.delim(gzfile(universe_table))
stopifnot(nrow(u) == 12L, sum(u$observed_positive) == 12L, sum(u$observed_bottom) == 1L)
has(html, 'Positions in U')

# 2b. Without level 1 and level 2: each says so and the rest renders
html <- render('partial.html', accumulation_dir = file.path(scratch, 'none'), within_gene_dir = file.path(scratch, 'none'))
has(html, 'Level 1: not available', 'Level 2: not available', 'was not found: the module ran', 'Overlaps and proximity (descriptive)')
lacks(html, 'Genes above their share of U')

# 2c. Validation mode: held-out table instead of term results
html <- render('holdout.html', within_gene_dir = file.path(scratch, 'within_gene_holdout'))
has(html, 'False-positive rate with held-out cycles', 'validation mode', 'Rate minus nominal level')

# 2d. A term that is not identifiable shows no p or q
wg <- read.delim(file.path(scratch, 'within_gene/within_gene_results.tsv'))
stopifnot(nrow(wg) > 0, 'TERM_A' %in% wg$term_id,
          all(is.na(wg$p[wg$status != 'tested'])))
cat('PASS: POSENRICH report renders with all components, without levels 1 and 2, and in validation mode\n')
