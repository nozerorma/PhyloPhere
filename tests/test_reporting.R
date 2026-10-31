# Reporting contracts: run from the repository root with Rscript tests/test_reporting.R.
options(phylo_debug = FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggtree)})
root <- normalizePath('.')
src <- file.path(root, 'subworkflows/TRAIT_ANALYSIS/local/src')
source(file.path(src, 'io_utils.R'))
source(file.path(src, 'stats.R'))
source(file.path(src, 'plotting_fun.R'))
source(file.path(src, 'palettes.R'))
source(file.path(src, 'report_style.R'))
source(file.path(src, 'report_plots.R'))
source(file.path(src, 'branch_cache.R'))
debug_log <- function(...) NULL
checks <- 0L
check <- function(name, expr) {
  force(expr)
  checks <<- checks + 1L
  cat('PASS:', name, '\n')
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({force(expr); NULL}, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message))
}
trait <- 'score'; taxon_of_interest <- 'family'; n_trait <- ''; c_trait <- ''
secondary_trait <- ''; branch_trait <- ''
fixture <- tempfile('reporting-'); dir.create(fixture)

check('numeric validation retains NA and rejects malformed or infinite values', {
  stopifnot(identical(numeric_trait_values(c('1', NA, '')), c(1, NA_real_, NA_real_)))
  for (bad in list(c('1', 'oops'), c(1, Inf), c(1, -Inf))) {
    expect_error(numeric_trait_values(bad, 'score'), 'score.*invalid/non-finite')
  }
})
check('descriptive summaries exist for one to three species', {
  for (n in 1:3) {
    d <- data.frame(species=paste0('S',1:n),family='F',score=seq_len(n))
    tab <- stats.f(d)
    stopifnot(nrow(tab)==n, all(tab$g_median==median(seq_len(n))),
              all(tab$g_q25==unname(quantile(seq_len(n),.25))))
  }
  expect_error(stats.f(transform(d, score=NA_real_)), 'No non-missing')
})
check('numeric validation precedes filtering in shared statistics', {
  d <- data.frame(species=paste0('S',1:4),family='F',score=c('1','2','bad','4'))
  expect_error(stats.f(d), 'invalid/non-finite')
  stopifnot(ci_overlap(0,1,1,2), !ci_overlap(0,1,2,3),
            is_ordinal_trait(c(0,1,0,1)), !is_ordinal_trait(c(.1,.2,.3,.4)))
})
check('configured headers and archive delimiters survive round trips', {
  d <- data.frame('sample id'=c('S1','S2'),family='F',score=c(1,2),check.names=FALSE)
  csv <- file.path(fixture,'input.csv'); txt <- file.path(fixture,'input.txt')
  out <- file.path(fixture,'original_trait_file.tsv')
  write.csv(d,csv,row.names=FALSE)
  write.table(d,txt,sep='\t',row.names=FALSE,quote=TRUE)
  stopifnot(identical(names(read_trait_table(txt)),names(d)))
  write_trait_archive(csv,out)
  expected <- as.data.frame(lapply(d,as.character),check.names=FALSE)
  stopifnot(identical(read_trait_table(out),expected), grepl('\t',readLines(out,n=1)))
})
check('unknown taxon colors use the shared fallback', {
  color_palette <- 'primates_palette'; clade_name <- 'primates'
  capitalized_trait <- 'Score'; capitalized_taxon <- 'Family'; capitalized_clade <- 'Primates'
  tableau_10 <- c(vermillion='#E15759',blue='#4E79A7',green='#59A14F')
  d <- stats.f(data.frame(species=paste0('S',1:4),family=c('Hominidae','Hominidae','NewFamily','NewFamily'),score=1:4)) %>%
    mutate(taxa=family)
  ordered_taxa <- unique(d$taxa)
  colors <- suppressWarnings(ggplot_build(contrast_plot.f(d))$data[[1]]$colour)
  stopifnot(colors[1]==primates_palette['Hominidae'],
            colors[3]==resolve_taxa_palette('NewFamily'), colors[3]!='grey50')
})
# Tree/trait provenance is separate from the public CSV schema.
tree <- ape::read.tree(text='((S1:1,S2:1):1,(S3:1,S4:1):1);')
values <- c(1,2,3,4)
check('branch cache verifies topology, lengths, trait values and schema', {
  p <- branch_cache_provenance(tree,tree$tip.label,values,'score')
  path <- file.path(fixture,'tree_with_traits.csv')
  cached <- data.frame(node=1:7,BR=1:7,species=c(tree$tip.label,rep(NA,3)))
  write.csv(cached,path,row.names=FALSE); saveRDS(p,paste0(path,'.branch_cache.rds'))
  stopifnot(nrow(read_branch_cache(path,p))==7)
  different <- ape::read.tree(text='(((S1:1,S2:1):1,S3:2):1,S4:3);')
  stopifnot(identical(tree$tip.label,different$tip.label))
  stopifnot(is.null(read_branch_cache(path,branch_cache_provenance(different,different$tip.label,values,'score'))))
  different <- tree; different$edge.length[1] <- 2
  stopifnot(is.null(read_branch_cache(path,branch_cache_provenance(different,different$tip.label,values,'score'))),
            is.null(read_branch_cache(path,branch_cache_provenance(tree,tree$tip.label,values+1,'score'))),
            is.null(read_branch_cache(path,branch_cache_provenance(tree,tree$tip.label,values,'other'))))
  unlink(paste0(path,'.branch_cache.rds')); stopifnot(is.null(read_branch_cache(path,p)))
})
check('ACE failure leaves uncertainty and derivedness unavailable', {
  tree8 <- ape::read.tree(text='(((S1:1,S2:1):1,(S3:1,S4:1):1):1,((S5:1,S6:1):1,(S7:1,S8:1):1):1);')
  d <- data.frame(species=paste0('S',1:8),family='F',score=c(1:7,100))
  fit <- function(tree,trait,model,...) list(opt=list(lnL=0,sigsq=1,alpha=.1,aicc=if(model=='BM')0 else 10))
  ace <- function(...) stop('injected ACE failure')
  result <- asr_tree.f(d,tree8,'score','family',fit_fun=fit,ace_fun=ace)
  stopifnot(all(is.na(result$node_table$ci_lower)),all(is.na(result$node_table$ci_upper)),
            all(is.na(result$tip_table$immediate_parent)),
            all(is.na(result$tip_table$cumulative_derived)))
})
# Exercise the real shared setup and pruning chunk in an isolated working directory.
old <- getwd(); setwd(fixture)
dir.create('src'); invisible(file.copy(list.files(src,full.names=TRUE),'src',recursive=TRUE))
tree8 <- ape::read.tree(text='(((S1:1,S2:1):1,(S3:1,S4:1):1):1,((S5:1,S6:1):1,(S7:1,S8:1):1):1);')
ape::write.tree(tree8,'tree.nwk')
write.csv(data.frame(species=paste0('S',1:8),family='F',score=c(1:7,100),N=10),'traits.csv',row.names=FALSE)
writeLines('S8','prune.txt')
params <- list(trait_file='traits.csv',tree_file='tree.nwk',output_dir='data_exploration',seed='1998',
               clade_name='primates',taxon_of_interest='family',sp_colname='species',traitname='score',
               n_trait='N',c_trait='',secondary_trait='',branch_trait='',trait_type='continuous',
               pss_top_pct='0.05',perm_strategy='best_model',max_contrasts='0',prune_list='prune.txt',prune_list_secondary='')
source('src/commons.R')
chunk <- function(file,label) {
  lines<-readLines(file); start<-grep(paste0('^```\\{r ',label,'[,}]'),lines)[1]+1
  end<-start+which(lines[start:length(lines)]=='```')[1]-2
  paste(lines[start:end],collapse='\n')
}
rmd <- file.path(root,'subworkflows/TRAIT_ANALYSIS/local/0.Data_pruning.Rmd')
check('pruned statistics describe the retained cohort', {
  dir.create(data_pruning_dir,recursive=TRUE)
  eval(parse(text=chunk(rmd,'set_ourdirs')))
  stopifnot(nrow(pruned_stats)==7,all(pruned_stats$g_median==4))
})
check('sample-size-only fan aggregation preserves N', {
  lines<-readLines(rmd); start<-grep('^if \\(isTRUE\\(has.n\\) && isTRUE\\(has.secondary\\)\\) \\{ # N and secondary trait',lines)[1]
  end<-start+grep('^tree_trait_data <- tree_trait_data %>%',lines[(start+1):length(lines)])[1]-1
  tree_trait_data<-stats_df %>% mutate(taxa=family,value=score) %>% group_by(species,taxa)
  eval(parse(text=paste(lines[start:end],collapse='\n')))
  stopifnot(all(tree_trait_data$n_sum==10))
})
setwd(old)
cat(checks,'reporting regression groups passed\n')
