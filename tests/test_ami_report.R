# Rscript tests/test_ami_report.R <repo root> <scratch directory>
# Edge cases of 13.AMI_analysis.Rmd against the minimal DOMINO fixture
# (tests/fixtures/ami): six nominal seeds, two BH hits, two network connectors;
# then an empty nominal list, which must skip cleanly instead of failing.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1])
work <- args[2]
fixture <- file.path(root,'tests','fixtures','ami')
render_ami <- function(wd, global_list) {
  dir.create(wd,recursive=TRUE,showWarnings=FALSE)
  writeLines(global_list,file.path(wd,'global.txt'))
  for (f in c('network.sif','background.txt','gene_scores.tsv')) file.copy(file.path(fixture,f),wd,overwrite=TRUE)
  file.copy(file.path(fixture,'modules'),wd,recursive=TRUE,overwrite=TRUE)
  file.copy(file.path(fixture,'cache'),wd,recursive=TRUE,overwrite=TRUE)
  file.copy(file.path(root,'subworkflows/ENRICHMENT/local/13.AMI_analysis.Rmd'),wd,overwrite=TRUE)
  owd <- setwd(wd); on.exit(setwd(owd))
  rmarkdown::render('13.AMI_analysis.Rmd',quiet=TRUE,
    params=list(background_file='background.txt',background_basename='background.txt',
      project_name='fixture',species=9606,domino_network_score_thr=700,
      gene_scores_file='gene_scores.tsv',string_db_dir=NULL,string_cache_dir='cache',
      domino_network_sif='network.sif',domino_modules_dir='modules',
      domino_edge_scores_file=NULL,seed='1998'),
    output_file='13.AMI_analysis_fixture.html',envir=new.env())
}

# Six nominal seeds (g1..g6), of which g1/g2 also pass BH; g7/g8 are connectors only.
render_ami(file.path(work,'nominal'),c('g1','g2','g3','g4','g5','g6'))
ev <- read.delim(file.path(work,'nominal','ami_networks','caas_global_gene_evidence.tsv'))
stopifnot(nrow(ev)==6L, all(ev$caas_nominal), sum(ev$caas_adjusted)==2L,
          !any(c('g7','g8') %in% ev$Gene), all(c('g1','g2') %in% ev$Gene[ev$caas_adjusted]))
desc <- read.delim(file.path(work,'nominal','ami_networks','module_descriptions_all_lists.tsv'))
stopifnot(nrow(desc)==1L, desc$n_genes==6L, !any(grepl('g7|g8',desc$genes)),
          file.exists(file.path(work,'nominal','13.AMI_analysis_fixture.html')))

# Empty nominal list: the report renders and skips the CAAS section instead of failing.
render_ami(file.path(work,'empty'),character())
stopifnot(file.exists(file.path(work,'empty','13.AMI_analysis_fixture.html')))
cat('PASS AMI report: nominal seeds vs BH hits vs connectors; empty list skips cleanly\n')
