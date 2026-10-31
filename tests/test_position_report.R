# Rscript tests/test_position_report.R <repo root> <scratch directory>
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1])
work <- args[2]
dir.create(work,recursive=TRUE,showWarnings=FALSE)
file.copy(file.path(root,'subworkflows/ENRICHMENT/local/14.Comparison_report.Rmd'),work,overwrite=TRUE)
dir.create(file.path(work,'src'),showWarnings=FALSE)
file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src',c('evidence_flags.R','fcs_selection.R','report_presentation.R','ami_evidence.R')),file.path(work,'src'),overwrite=TRUE)
writeLines(c('Gene\tPosition\tside\tCAAS_score\tp.emp\tp.adj_bh\tcaas_nominal\tcaas_adjusted','G1\t0\ttop\t0.9\t0.01\t0.02\tTRUE\tTRUE','G1\t0\tbottom\t0.8\t0.01\t0.02\tTRUE\tTRUE','G2\t1\ttop\t0.7\t0.7\t0.7\tFALSE\tFALSE'),file.path(work,'scores.tsv'))
writeLines(c('Gene\tPosition\tside\tFeature\tConsequence\tIMPACT\tMANE_SELECT\tstatus\tcoverage\tsource\tvariant_id\tcandidate_count\tannotated_count',
'G1\t0\ttop\tENST1\tmissense_variant\tMODERATE\tNM1\tannotated\tcomplete\tensembl\tchr1:1:A:T\t1\t1',
'G1\t0\ttop\tENST2\tstop_gained\tHIGH\t-\tannotated\tcomplete\tensembl\tchr1:1:A:T\t1\t1',
'G1\t0\tbottom\tENST3\tsplice_region_variant\tLOW\t-\tannotated\tcomplete\tensembl\tchr1:1:A:T\t1\t1'),file.path(work,'vep.tsv'))
writeLines(c('gene,position,max_bf,target_aa','G1,0,100,S'),file.path(work,'fade.csv'))
writeLines('gene,position,max_bf,target_aa',file.path(work,'fade_bottom.csv'))
writeLines(c('Gene\tPosition\tpfam_domain\tpfam_clan\tucr_region\tposition_variability\tfubar_selection',
'G1\t0\tDomain1;Domain2\tClan1\tcore\t0.1\tnegative'),file.path(work,'char.tsv'))
env <- new.env()
rmarkdown::render(file.path(work,'14.Comparison_report.Rmd'),output_dir=work,quiet=TRUE,envir=env,
 params=list(traitname='annotation_fixture',fcs_dir=file.path(work,'absent_fcs'),fcs_le_dir=file.path(work,'absent_le'),
 position_scores_file=file.path(work,'scores.tsv'),vep_ensembl_file=file.path(work,'vep.tsv'),
 fade_sites_top_file=file.path(work,'fade.csv'),fade_sites_bottom_file=file.path(work,'fade_bottom.csv'),position_characterization_file=file.path(work,'char.tsv')))
env$params <- list(top_n=100)
old_eval_wd <- getwd(); setwd(normalizePath(work))
top <- env$build_interesting_positions('Top')
bottom <- env$build_interesting_positions('Bottom')
# Return value carries flags separately from its rendered table.
str(top,max.level=1)
saveRDS(list(top=top,bottom=bottom),file.path(work,'tables.rds'))

stopifnot(top$data$VEP_consequence[top$data$Gene=='G1']=='missense_variant',
 bottom$data$VEP_consequence[bottom$data$Gene=='G1']=='splice_region_variant',
 top$data$`Pfam domain`[top$data$Gene=='G1']=='Domain1;Domain2',
 top$data$FADE_overlap[top$data$Gene=='G1'],
 top$data$`FADE site BF`[top$data$Gene=='G1']==100,
 top$data$`FADE target AA`[top$data$Gene=='G1']=='S',
 top$data$`FUBAR selection`[top$data$Gene=='G1']=='purifying',
 !bottom$data$FADE_overlap[bottom$data$Gene=='G1'])
setwd(old_eval_wd)
writeLines('VEP_REPORT_PASS: MANE preference, side-specific joins, Pfam overlap and FADE sites',file.path(work,'validation.txt'))

# Sin position_characterization.tsv (POSENRICH no corrió) el informe renderiza y omite solo las columnas de la caracterización.
work2 <- file.path(work, 'without_characterization')
dir.create(file.path(work2, 'src'), recursive=TRUE, showWarnings=FALSE)
file.copy(file.path(work, '14.Comparison_report.Rmd'), work2, overwrite=TRUE)
file.copy(list.files(file.path(work, 'src'), full.names=TRUE), file.path(work2, 'src'), overwrite=TRUE)
for (f in c('scores.tsv', 'vep.tsv', 'fade.csv', 'fade_bottom.csv')) file.copy(file.path(work, f), work2, overwrite=TRUE)
env2 <- new.env()
rmarkdown::render(file.path(work2, '14.Comparison_report.Rmd'), output_dir=work2, quiet=TRUE, envir=env2,
 params=list(traitname='annotation_fixture', fcs_dir=file.path(work2, 'absent_fcs'), fcs_le_dir=file.path(work2, 'absent_le'),
 position_scores_file=file.path(work2, 'scores.tsv'), vep_ensembl_file=file.path(work2, 'vep.tsv'),
 fade_sites_top_file=file.path(work2, 'fade.csv'), fade_sites_bottom_file=file.path(work2, 'fade_bottom.csv'),
 position_characterization_file=file.path(work2, 'absent_characterization.tsv')))
env2$params <- list(top_n=100)
setwd(normalizePath(work2))
top2 <- env2$build_interesting_positions('Top')
stopifnot(top2$data$VEP_consequence[top2$data$Gene=='G1']=='missense_variant',
 !'Pfam domain' %in% names(top2$data), !'FUBAR selection' %in% names(top2$data))
setwd(old_eval_wd)
writeLines('NO_CHARACTERIZATION_REPORT_PASS', file.path(work2, 'validation.txt'))
