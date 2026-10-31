# Saved report inputs make standalone regeneration use the original settings.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5)
dir.create('fade_inputs', showWarnings=FALSE)
file.copy(args[2], 'fade_inputs/gene_counts.tsv', overwrite=TRUE)
file.copy(args[3], 'fade_inputs/sites.csv', overwrite=TRUE)
if (!startsWith(basename(args[4]), 'NO_')) file.copy(args[4], 'fade_inputs/foreground.txt', overwrite=TRUE)
if (!startsWith(basename(args[5]), 'NO_')) file.copy(args[5], 'fade_inputs/settings.json', overwrite=TRUE)
rmarkdown::render('6.FADE_report.Rmd', params=list(
  gene_counts_file='fade_inputs/gene_counts.tsv', sites_file='fade_inputs/sites.csv',
  fg_list_file=if (file.exists('fade_inputs/foreground.txt')) 'fade_inputs/foreground.txt' else NULL,
  metadata_file=if (file.exists('fade_inputs/settings.json')) 'fade_inputs/settings.json' else NULL,
  direction=args[1]), output_file=paste0('6.FADE_report_',args[1],'.html'))
