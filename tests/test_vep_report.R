args <- commandArgs(trailingOnly=TRUE)
root <- args[1];work <- args[2]
library(dplyr)
s <- readLines(file.path(root,'subworkflows/ENRICHMENT/local/14.Comparison_report.Rmd'))
.file_ok <- function(p) !is.null(p) && file.exists(p) && file.info(p)$size > 0
start <- match('# <vep_contract>',s);end <- match('# </vep_contract>',s)
eval(parse(text=s[(start+1):(end-1)]))
write_literal <- function(path,d) {
 writeLines(paste(names(d),collapse='\t'),path)
 for(i in seq_len(nrow(d)))cat(paste(as.character(d[i,]),collapse='\t'),'\n',file=path,append=TRUE,sep='')
}
mk <- function(source,status='annotated',coverage='complete',side='top',gene='"Q"') {
 data.frame(Gene=gene,Position='0',side=side,status=status,coverage=coverage,source=source,
 variant_id=if(status=='annotated')'chr1:100:G:A' else '',candidate_count='1',annotated_count=if(status=='annotated')'1' else '0',
 prediction='benign',stringsAsFactors=FALSE)
}
p <- file.path(work,'annotations.tsv')
write_literal(p,rbind(mk('primateai'),mk('primateai',side='bottom',gene='001')))
d <- read_vep_contract(p,'primateai')
stopifnot(identical(d$Gene,c('"Q"','001')))
base <- data.frame(gene=c('"Q"','"Q"','001','NA'),Position=0L,side=c('top','bottom','bottom','top'))
stopifnot(identical(vep_status_for(d,base),c('annotated','target_absent','annotated','target_absent')))
# Old unverified tables cannot participate silently.
writeLines(c('Gene\tPosition\tside\tConsequence','001\t0\ttop\tmissense_variant'),p)
warn <- FALSE
legacy <- withCallingHandlers(read_vep_contract(p,'ensembl'),warning=function(w){warn<<-TRUE;invokeRestart('muffleWarning')})
stopifnot(is.null(legacy),warn)
write_literal(p,rbind(mk('cosmic','no_hit','snapshot_searched'),mk('cosmic','unmapped_no_map','not_evaluated',side='bottom')))
vep_cosmic_cmp_df <- read_vep_contract(p,'cosmic');has_vep_cosmic_cmp <- TRUE
write_literal(p,mk('primateai',coverage='partial'))
vep_primateai_cmp_df <- read_vep_contract(p,'primateai');has_vep_primateai_cmp <- TRUE
vep_ensembl_cmp_df <- NULL;has_vep_ensembl_cmp <- FALSE
start <- grep('^  # Human-reference candidate annotations',s)
end <- grep('^  # Position characterization GMTs:',s)
block <- s[start:(end-1)]
eval(parse(text=block))
stopifnot(identical(base$COSMIC,c('no_hit','missing','missing','missing')),
          all(base$PAI3D_pathogenic=='missing'),all(base$VEP_consequence=='missing'))
# Complete benign is side-specific; partial pathogenic remains descriptive positive evidence.
write_literal(p,rbind(mk('primateai'),transform(mk('primateai',coverage='partial',side='bottom'),prediction='pathogenic')))
vep_primateai_cmp_df <- read_vep_contract(p,'primateai')
eval(parse(text=block))
stopifnot(identical(base$PAI3D_pathogenic,c('benign','pathogenic','missing','missing')))
# Duplicate status rows are invalid, whereas different annotated SNVs are allowed.
write_literal(p,rbind(mk('cosmic','no_hit','snapshot_searched'),mk('cosmic','no_hit','snapshot_searched')))
stopifnot(inherits(try(read_vep_contract(p,'cosmic'),silent=TRUE),'try-error'))
cat('VEP_REPORT_PASS\n')
