# Run from repository root; optional argument: saved marine_complete directory.
suppressPackageStartupMessages({ library(dplyr); library(Matrix) })
source('subworkflows/ENRICHMENT/local/src/evidence_flags.R')
source('subworkflows/ENRICHMENT/local/src/fcs_enrich.R')
source('subworkflows/ENRICHMENT/local/src/gene_p.R')
same <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr), error=identity), 'error'))
# The observed writers must use the same unbounded transform as the native null.
r <- data.frame(Rho=c(1,-1,0,1,0,NA),P=c(0,0,0,1e-310,1,0), row.names=paste0('g',1:6))
expected <- RERconverge::getStat(r)
nf <- readLines('subworkflows/RERCONVERGE/rer_gene_lists.nf')
line <- gsub(intToUtf8(92), '', trimws(grep('log_p    <-',nf,value=TRUE)), fixed=TRUE)
df <- r; eval(parse(text=line)); got <- fcs_build_vals(setNames(sign(df$Rho)*log_p,rownames(df)),rownames(df),score_source='rer')
same(got,expected); stopifnot(identical(names(got),c('g1','g2','g4','g5')),got['g4']>300,is.infinite(got['g1']),got['g5']==0)
s <- readLines('subworkflows/SCORING/local/src/scoring_compute.R')
gene_scores <- data.frame(rer_P=r$P); .col <- function(x,n) x[[n]]
eval(parse(text=trimws(grep('^  \\.rp <-',s,value=TRUE))))
same(.rp,r$P)
report <- readLines('subworkflows/RERCONVERGE/local/5.RERconverge_report.Rmd')
stopifnot(any(grepl('rer_df$rer_stat <- sign(rer_df$Rho) * -log10(rer_df$P)',report,fixed=TRUE)))
# Exact duplicates, including a missing score, must fail; suffixes remain distinct IDs.
reject(fcs_build_vals(setNames(c(NA,2),c('g','g')), 'g',score_source='rer'))
reject(fcs_build_vals(setNames(c(1,2),c('g','g')), 'g'))
reject(fcs_build_vals(setNames(c(1,2),c('g.1','g.2')), character(),score_source='rer'))
reject(fcs_check_rer_gene_ids(c('g','g.1'),2,'RER null matrix'))
same(names(fcs_build_vals(setNames(c(1,2),c('g.1','h.2')), character(),score_source='rer')),c('g.1','h.2'))
reject(fcs_align_null(matrix(1,2,2,dimnames=list(c('g','g'),NULL)),'g'))
reject(fcs_null_enrichstat_vectorized(matrix(1,2,2,dimnames=list(c('g','g'),NULL)),list(),list()))
cat('PASS: native extreme-P transform in observed writers; exact duplicates rejected; full IDs preserved\n')
# Native AUC with infinities and measured zeros remains well-defined.
vals <- setNames(c(Inf,-Inf,310,0,2,1,0,0),paste0('g',1:8))
gs <- list(A=names(vals)[1:4],B=names(vals)[5:8]); gmts <- list(db=list(genesets=gs,geneset.names=names(gs)))
obs <- fcs_run_ranking(vals,gmts,num_g=2,alternative='two.sided')
real <- list(db=data.frame(stat=obs$stat,row.names=obs$pathway))
x <- matrix(vals,ncol=1,dimnames=list(names(vals),NULL))
v <- fcs_null_enrichstat_vectorized(x,gmts,real,num_g=2)$db
same(v[,1],obs$stat); stopifnot(all(is.finite(v)))
cat('PASS: observed/native and vectorized null AUC agree with infinite gene scores\n')
a <- commandArgs(TRUE)
if(length(a)) {
 root <- a[1]; null <- readRDS(file.path(root,'caas_permulation/caas_perms.rds'))
 stats <- read.delim(file.path(root,'scoring/fcs_stats.tsv'))
 invisible(capture.output(gmts <- fcs_load_gmts(file.path(root,'core_inputs/gmt'))))
 js <- unique(round(seq(1,100,length.out=12)))
 cells <- 0; worst <- 0; p_worst <- 0
 for(rk in c('global','top','bottom')) for(mode in c('raw','gene_p')) {
  full <- null$caas_corStat_byrank[[rk]]
  obs <- setNames(stats[[paste0('score_',rk)]],stats$gene)[rownames(full)]
  if(mode=='gene_p') { gp <- gene_p(full,obs); obs <- setNames(gp$nlog10_obs,rownames(full)); full <- gp$nlog10_null }
  vals <- fcs_build_vals(obs,rownames(full)); m <- fcs_align_null(full[,js,drop=FALSE],names(vals))
  observed <- fcs_run_ranking(vals,gmts,num_g=10,alternative='greater')
  if (!nrow(observed)) {cat(sprintf('PASS: saved %s/%s has no evaluable pathways\n',rk,mode));next}
  real <- lapply(split(as.data.frame(observed),observed$database),function(d){rownames(d)<-d$pathway;d})
  vec <- fcs_null_enrichstat_vectorized(m,gmts,real,num_g=10)
  refs <- lapply(vec,function(x) x*NA_real_)
  for(j in seq_len(ncol(m))) {
   invisible(capture.output(native <- RERconverge::fastwilcoxGMTall(setNames(m[,j],rownames(m)),gmts,outputGeneVals=FALSE,num.g=10)))
   for(db in names(vec)) refs[[db]][,j] <- native[[db]]$stat[match(rownames(vec[[db]]),rownames(native[[db]]))]
  }
  for(db in names(vec)) {same(vec[[db]],refs[[db]]); cells<-cells+sum(!is.na(vec[[db]])); worst<-max(worst,abs(vec[[db]]-refs[[db]]),na.rm=TRUE)}
  counts <- fcs_null_counts_by_block(rk,mutate(observed,ranking=rk),m,gmts,'greater',num_g=10)
  blocked <- fcs_null_counts_by_block(rk,mutate(observed,ranking=rk),m,gmts,'greater',num_g=10,block_cols=5)
  counts<-arrange(counts,database,pathway); blocked<-arrange(blocked,database,pathway);same(counts,blocked)
  for(db in names(real)) {
   obsstat<-real[[db]]$stat
   expected <- (1+rowSums(refs[[db]] >= obsstat,na.rm=TRUE))/(1+rowSums(!is.na(refs[[db]])))
   got <- filter(counts,database==db); expected<-expected[match(got$pathway,rownames(refs[[db]]))]
   actual<-(got$count+1)/(got$valid+1);same(actual,expected);p_worst<-max(p_worst,abs(actual-expected),na.rm=TRUE)
  }
  cat(sprintf('PASS: marine %s/%s, native AUC and independent p counts, block size 5 vs whole\n',rk,mode))
 }
 cat(sprintf('TOTAL: %d valid native AUC comparisons, max delta %.3g; max p delta %.3g\n',cells,worst,p_worst))
}
