# Render del informe POSENRICH modular sobre un pipeline dirigido de fixture.
# Uso: Rscript tests/test_posenrich_modular_report.R <repo> <scratch>
a <- commandArgs(TRUE)
root <- normalizePath(a[1]); scratch <- a[2]
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
src <- file.path(root, 'subworkflows/ENRICHMENT/local/src')

# 1. Fixture: dos genes x seis posiciones, miembros localizados (reutiliza la lógica
#    del fixture de localización construyendo los archivos directamente).
maps <- file.path(scratch, 'maps'); dir.create(maps)
for (gene in c('G1', 'G2')) {
  chrom <- if (gene == 'G1') 'chr1' else 'chr2'
  lines <- c('ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord',
             sprintf('%d\tselected\t%d\t%d\t%s:%d', 1:6, 1:6, 1:6, chrom, 1000 * 1:6))
  writeLines(lines, file.path(maps, paste0(gene, '.map.tsv')))
}
writeLines(c('G1', 'G2'), file.path(scratch, 'bg.txt'))
g1 <- c(0.9, 0.8, 0.7, 0.2, 0.1, 0.05); g2 <- c(0.5, 0.4, 0.3, 0.25, 0.2, 0.1)
rows <- c('Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\tnd\tnull_unseen\tcaas_nominal\tcaas_adjusted')
for (gene in c('G1', 'G2')) for (pos in 0:5) {
  sc <- if (gene == 'G1') g1[pos + 1] else g2[pos + 1]
  p <- if (gene == 'G1' && pos < 2) '0.01' else '0.2'
  rows <- c(rows, sprintf('%s\t%d\t%s\ttop\t%s\t0.5\t2\tFALSE\tFALSE\tFALSE', gene, pos, sc, p))
}
writeLines(rows, file.path(scratch, 'position_scores.tsv'))
custom <- c('G1\t0\tTERM_A\tlocalized high', 'G1\t1\tTERM_A\tlocalized high',
            'G2\t3\tTERM_A\tlocalized low', 'G2\t4\tTERM_A\tlocalized low',
            paste(sprintf('G1\t%d\tTERM_FULL\twhole protein', 0:5), collapse = '\n'))
writeLines(custom, file.path(scratch, 'custom.tsv'))

py <- Sys.getenv('PYTHON'); if (!nzchar(py)) py <- 'python3'
run_py <- function(script, args) {
  rc <- system2(py, c(file.path(src, script), args))
  if (rc != 0) stop(paste('fallo', script))
}
run_py('posenrich_registry.py', c('--map_dir', maps, '--cleaned_background', file.path(scratch, 'bg.txt'),
                                 '--position_scores', file.path(scratch, 'position_scores.tsv'),
                                 '--custom_marker_file', file.path(scratch, 'custom.tsv'),
                                 '--output_dir', file.path(scratch, 'registry')))
run_py('posenrich_characterize.py', c('--registry_dir', file.path(scratch, 'registry'),
                                      '--output_dir', file.path(scratch, 'characterization')))
run_py('posenrich_localization.py', c('--registry_dir', file.path(scratch, 'registry'),
                                      '--output_dir', file.path(scratch, 'localization'),
                                      '--n_permutations', '199', '--seed', '7'))
run_py('posenrich_global.py', c('--registry_dir', file.path(scratch, 'registry'),
                                '--output_dir', file.path(scratch, 'global')))
run_py('posenrich_proximity.py', c('--registry_dir', file.path(scratch, 'registry'),
                                   '--output_dir', file.path(scratch, 'proximity')))

# 2. Render del informe modular.
suppressPackageStartupMessages(library(rmarkdown))
out_html <- file.path(scratch, 'POSENRICH_modular_report_trait.html')
rmarkdown::render(file.path(root, 'subworkflows/ENRICHMENT/local/posenrich_modular_report.Rmd'),
                  params = list(registry_dir = file.path(scratch, 'registry'),
                                characterization_dir = file.path(scratch, 'characterization'),
                                localization_dir = file.path(scratch, 'localization'),
                                global_dir = file.path(scratch, 'global'),
                                proximity_dir = file.path(scratch, 'proximity'),
                                traitname = 'trait', ranking = 'caas_score', aggregation = 'mean',
                                n_permutations = 199, seed = 7, radios = '0,5,10'),
                  output_file = out_html, quiet = TRUE)
stopifnot(file.exists(out_html))
html <- paste(readLines(out_html, warn = FALSE), collapse = '\n')
for (needle in c('POSENRICH modular', 'Localización dentro de proteínas',
                 'Asociación global', 'Proximidad', 'registro', 'TERM_A',
                 '0.625', 'not_identifiable', 'descriptive', 'not_validated')) {
  stopifnot(grepl(needle, html, fixed = TRUE))
}
# El informe distingue los estados de cada componente.
stopifnot(grepl('ejecutado', html, fixed = TRUE))
cat('PASS: render del informe POSENRICH modular con pipeline dirigido de fixture\n')
