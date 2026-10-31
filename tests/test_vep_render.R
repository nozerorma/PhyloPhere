args <- commandArgs(trailingOnly=TRUE)
root <- args[1];work <- normalizePath(args[2])
file.copy(file.path(root,'subworkflows/ENRICHMENT/local/14.Comparison_report.Rmd'),work,overwrite=TRUE)
dir.create(file.path(work,'src'),showWarnings=FALSE)
for (f in c('evidence_flags.R','report_presentation.R','ami_evidence.R','fcs_selection.R')) {
 file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src',f),file.path(work,'src'),overwrite=TRUE)
}
e <- new.env()
rmarkdown::render(file.path(work,'14.Comparison_report.Rmd'),output_dir=work,quiet=TRUE,envir=e,
 params=list(traitname='VEP_fixture',fcs_dir=file.path(work,'absent_fcs'),fcs_le_dir=file.path(work,'absent_le'),
 position_scores_file=file.path(work,'scores.tsv'),vep_primateai_file=file.path(work,'pai.tsv'),
 vep_cosmic_file=file.path(work,'cosmic_out.tsv'),vep_ensembl_file=file.path(work,'ens.tsv')))
old <- setwd(work)
on.exit(setwd(old))
res <- e$build_interesting_positions('Top')
stopifnot(identical(as.character(res$data$Gene),'G.internal.2'),res$data$COSMIC=='hit',
 res$data$VEP_consequence=='missense_variant',res$data$PAI3D_pathogenic=='pathogenic',
 res$data$COSMIC_status=='annotated',res$data$PAI3D_coverage=='complete')
# Re-render from published inputs without invoking annotation again.
rmarkdown::render(file.path(work,'14.Comparison_report.Rmd'),output_dir=work,quiet=TRUE,envir=new.env(),
 params=list(traitname='VEP_regenerated',fcs_dir=file.path(work,'absent_fcs'),fcs_le_dir=file.path(work,'absent_le'),
 position_scores_file=file.path(work,'scores.tsv'),vep_primateai_file=file.path(work,'pai.tsv'),
 vep_cosmic_file=file.path(work,'cosmic_out.tsv'),vep_ensembl_file=file.path(work,'ens.tsv')))
cat('VEP_RENDER_PASS\n')
