# Render of the POSENRICH report (14.POSENRICH_report.Rmd) on a minimal pipeline fixture, with every component,
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

# 1. Fixture: two genes x six positions, a custom term localized inside G1 and G2, a whole-protein term.
maps <- file.path(scratch, 'maps'); dir.create(maps, showWarnings = FALSE)
for (gene in c('G1', 'G2')) {
  chrom <- if (gene == 'G1') 'chr1' else 'chr2'
  writeLines(c('ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord',
               sprintf('%d\tselected\t%d\t%d\t%s:%d', 1:6, 1:6, 1:6, chrom, 1000 * 1:6)),
             file.path(maps, paste0(gene, '.map.tsv')))
}
writeLines(c('G1', 'G2'), file.path(scratch, 'bg.txt'))
writeLines(c('G1\t0,1,2,3,4,5', 'G2\t0,1,2,3,4,5'), file.path(scratch, 'background.output'))
g1 <- c(0.9, 0.8, 0.7, 0.2, 0.1, 0.05); g2 <- c(0.5, 0.4, 0.3, 0.25, 0.2, 0.1)
rows <- c('Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\tnd\tnull_unseen\tcaas_nominal\tcaas_adjusted')
for (gene in c('G1', 'G2')) for (pos in 0:5) {
  sc <- if (gene == 'G1') g1[pos + 1] else g2[pos + 1]
  side <- if (gene == 'G2' && pos == 5) 'bottom' else 'top'
  p <- if (gene == 'G1' && pos < 2) '0.01' else '0.2'
  rows <- c(rows, sprintf('%s\t%d\t%s\t%s\t%s\t0.5\t2\tFALSE\tFALSE\tFALSE', gene, pos, sc, side, p))
}
writeLines(rows, file.path(scratch, 'position_scores.tsv'))
writeLines(c('G1\t0\tTERM_A\tlocalized high', 'G1\t1\tTERM_A\tlocalized high',
             'G2\t3\tTERM_A\tlocalized low', 'G2\t4\tTERM_A\tlocalized low',
             paste(sprintf('G1\t%d\tTERM_FULL\twhole protein', 0:5), collapse = '\n')),
           file.path(scratch, 'custom.tsv'))

# One block of 40 exchangeable null cycles, and the positions the null scores (null_positions.parquet).
blk <- file.path(scratch, 'blk_0001')
rc <- system2(py, c('-c', shQuote(paste(
  'import sys, os, numpy as np, pandas as pd',
  'rng = np.random.default_rng(1)',
  'rows = [(g, f"b_{c}", p, float(rng.uniform(0.05, 1))) for c in range(1, 41) for g in ("G1", "G2")',
  '        for p in range(6) if rng.random() < 0.4]',
  'os.makedirs(sys.argv[1], exist_ok=True)',
  'df = pd.DataFrame(rows, columns=["gene", "cycle", "position", "all"]).astype({"position": "int32"})',
  'df.to_parquet(sys.argv[1] + "/positions.parquet")',
  'open(sys.argv[1] + "/cycles.txt", "w").write("\\n".join(f"b_{c}" for c in range(1, 41)) + "\\n")',
  'nd = df.groupby(["gene", "position"]).size().rename("nd").reset_index()',
  'nd.to_parquet(sys.argv[2], index=False)', sep = '\n')), blk, file.path(scratch, 'null_positions.parquet')))
if (rc != 0) stop('null block construction failed')

run_py('posenrich_registry.py', c('--map_dir', maps, '--cleaned_background', file.path(scratch, 'bg.txt'),
                                 '--position_scores', file.path(scratch, 'position_scores.tsv'),
                                 '--custom_marker_file', file.path(scratch, 'custom.tsv'),
                                 '--output_dir', file.path(scratch, 'registry')))
run_py('posenrich_characterize.py', c('--registry_dir', file.path(scratch, 'registry'),
                                      '--output_dir', file.path(scratch, 'characterization')))
run_py('posenrich_universe.py', c('--background-genes', file.path(scratch, 'bg.txt'),
                                 '--background-positions', file.path(scratch, 'background.output'),
                                 '--null-positions', file.path(scratch, 'null_positions.parquet'),
                                 '--observed', file.path(scratch, 'position_scores.tsv'),
                                 '--output-dir', file.path(scratch, 'universe')))
universe_table <- file.path(scratch, 'universe/universe.tsv.gz')
run_py('accumulation_test.py', c('--universe', universe_table, '--blocks', blk, '--output-dir', file.path(scratch, 'accumulation')))
run_py('posenrich_within_gene.py', c('--universe', universe_table, '--observed', file.path(scratch, 'position_scores.tsv'),
                                     '--registry-dir', file.path(scratch, 'registry'), '--blocks', blk,
                                     '--min-genes', '1', '--output-dir', file.path(scratch, 'within_gene')))
run_py('posenrich_within_gene.py', c('--universe', universe_table, '--observed', file.path(scratch, 'position_scores.tsv'),
                                     '--registry-dir', file.path(scratch, 'registry'), '--blocks', blk,
                                     '--min-genes', '1', '--holdout', '10',
                                     '--output-dir', file.path(scratch, 'within_gene_holdout')))
run_py('posenrich_proximity.py', c('--registry_dir', file.path(scratch, 'registry'),
                                   '--output_dir', file.path(scratch, 'proximity')))

# 2. Render. The Rmd is rendered in place, as the process does in its work directory.
suppressPackageStartupMessages(library(rmarkdown))
rmd <- file.path(root, 'subworkflows/POSENRICH/local/14.POSENRICH_report.Rmd')
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
    'Gene-list check', 'Genes above the null share', 'Overlaps (exact identity)', 'Near positions')
lacks(html, 'mathjax.rstudio.com', 'not available')
# Numbers come from the tables: |U| = 12 positions, 12 observed positives
u <- read.delim(gzfile(universe_table))
stopifnot(nrow(u) == 12L, sum(u$observed_positive) == 12L, sum(u$observed_bottom) == 1L)
has(html, 'Positions in U')

# 2b. Without level 1 and level 2: each says so and the rest renders
html <- render('partial.html', accumulation_dir = file.path(scratch, 'none'), within_gene_dir = file.path(scratch, 'none'))
has(html, 'Level 1: not available', 'Level 2: not available', 'was not found: the module ran', 'Overlaps and proximity (descriptive)')
lacks(html, 'Genes above the null share')

# 2c. Validation mode: held-out table instead of term results
html <- render('holdout.html', within_gene_dir = file.path(scratch, 'within_gene_holdout'))
has(html, 'False-positive rate with held-out cycles', 'validation mode', 'Rate minus nominal level')

# 2d. A term that is not identifiable shows no p or q
wg <- read.delim(file.path(scratch, 'within_gene/within_gene_results.tsv'))
stopifnot(nrow(wg) > 0, 'TERM_A' %in% wg$term_id,
          all(is.na(wg$p[wg$status != 'tested'])))
cat('PASS: POSENRICH report renders with all components, without levels 1 and 2, and in validation mode\n')
