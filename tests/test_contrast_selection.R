# Run from the repository root: conda run -n phylophere Rscript tests/test_contrast_selection.R
options(phylo_debug = FALSE)
repo <- normalizePath('.')
source('subworkflows/TRAIT_ANALYSIS/local/src/selection_algorithm.R')
source('subworkflows/TRAIT_ANALYSIS/local/src/contrast_inputs.R')
checks <- 0L
check <- function(name, expr) { force(expr); checks <<- checks + 1L; cat('PASS:', name, '\n') }
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message))
}
tr <- ape::read.tree(text='(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):2);')
v <- setNames(c(.1,.9,.2,.8,.3,.7),tr$tip.label)
cv <- ape::vcv.phylo(tr); D <- ape::cophenetic.phylo(tr)
params <- list(seed='1998',trait_type='continuous',perm_strategy='bm',pss_top_pct='0',max_contrasts='0',max_fop='3')
d <- data.frame(species=names(v),score=unname(v),n=100,c=round(v*100))
check('helpers load from a relocated checkout and an unrelated working directory', {
  clone <- tempfile('contrast-portable-');dir.create(clone)
  src_copy <- file.path(clone,'subworkflows/TRAIT_ANALYSIS/local/src')
  ct_copy <- file.path(clone,'subworkflows/CT/local/scripts')
  dir.create(src_copy,recursive=TRUE);dir.create(ct_copy,recursive=TRUE)
  file.copy(file.path(repo,'subworkflows/TRAIT_ANALYSIS/local/src/selection_algorithm.R'),src_copy)
  file.copy(file.path(repo,'subworkflows/CT/local/scripts',c('lean_contrast_selector.R','pss_core.R')),ct_copy)
  # A fresh R process avoids accidentally using the engine already loaded above.
  test_script <- file.path(clone,'load.R')
  writeLines(c('setwd(tempdir())',paste0('source(',deparse(file.path(src_copy,'selection_algorithm.R')),')'),
    'stopifnot(exists("calculate_pairwise_scores"),exists("fop_pair_sel.f"))'),test_script)
  status <- system2(file.path(R.home('bin'),'Rscript'),shQuote(test_script))
  stopifnot(status==0)
})
check('invalid count columns, domains, modes, fractions and caps fail explicitly', {
  validate_contrast_inputs(d,'score','n','c',params)
  validate_contrast_inputs(d,'score','n','c',modifyList(params,list(perm_strategy='auto')))
  expect_error(validate_contrast_inputs(d,'score','typo','c',params),'not found')
  expect_error(validate_contrast_inputs(d,'score','','c',params),'requires n_trait')
  for (bad in list(transform(d,n=0),transform(d,c=-1),transform(d,c=101),transform(d,c=.5)))
    expect_error(validate_contrast_inputs(bad,'score','n','c',params),'Invalid integer counts')
  for (field in c('seed','max_contrasts','max_fop','trait_type','perm_strategy','pss_top_pct')) {
    bad <- params;bad[[field]] <- 'invalid'
    expect_error(validate_contrast_inputs(d,'score','n','c',bad),'.')
  }
  expect_error(validate_contrast_inputs(d,'score','n','c',modifyList(params,list(max_fop='0'))),'max_fop')
})
check('zero percentile disables filtering and invalid fractions are rejected', {
  cc <- function(p) lean_candidate_df(v,D,1,tr,cv,cv,'BM',top_pct=p,ordinal=FALSE)
  stopifnot(nrow(cc(0)$cand_df)==15,nrow(cc(1)$cand_df)==15,nrow(cc(.2)$cand_df)==3)
  for (bad in c(-1,2,NA,Inf)) expect_error(cc(bad),'top_pct')
})
check('percentile ties use species names independently of tree order', {
  star <- ape::read.tree(text='(A:1,B:1,C:1,D:1);')
  vv <- setNames(c(0,1,0,1),star$tip.label); vc <- ape::vcv.phylo(star)
  choose <- function(tree,values,cov) lean_candidate_df(values,ape::cophenetic.phylo(tree),1,
    tree,cov,cov,'BM',top_pct=.25,ordinal=FALSE)$cand_df[,c('species1','species2')]
  a <- choose(star,vv,vc)
  star$tip.label <- rev(star$tip.label);vc <- ape::vcv.phylo(star)
  b <- choose(star,vv[star$tip.label],vc)
  stopifnot(identical(unname(as.matrix(a)),unname(as.matrix(b))))
})
check('final Dunn indices match the entire selected set and domain IDs survive export', {
  cc <- lean_candidate_df(v,D,1,tr,cv,cv,'BM',top_pct=0,ordinal=FALSE)
  sel <- greedy_dunn_select(cc$cand_df,D)
  stopifnot(nrow(sel$selected)>1, all.equal(sel$selected$Dunn_index,
    vapply(seq_along(sel$members),function(k)mod_dunn_lean(D,sel$members,k),numeric(1))))
  export <- gen_traitfiles(sel$selected)
  stopifnot(nrow(export)==2*nrow(sel$selected),all(table(export$pair)==2))
  stopifnot(nrow(gen_traitfiles(data.frame(species1=character(),species2=character())))==0)
})
ci <- data.frame(species=names(v),trait=unname(v),n_population=NA_real_,n_cases=NA_real_,
                 trait_ci_lb=NA_real_,trait_ci_ub=NA_real_)
check('small and invariant cohorts return zero candidates without model fitting', {
  for (n in 1:2) {
    tiny <- ape::keep.tip(tr,names(v)[seq_len(n)])
    result <- prepare_observed_contrasts(ci[seq_len(n),],tiny,0,FALSE,'none')
    stopifnot(is.null(result$candidates$cand_df),is.null(result$ctx$fits))
  }
  flat <- ci;flat$trait <- .5
  stopifnot(is.null(prepare_observed_contrasts(flat,tr,0,FALSE,'none')$candidates$cand_df))
})
check('model-fitting failures are explicit and never broaden the candidate pool', {
  old_fit <- fit_models
  fit_models <- function(...)stop('injected optimizer failure')
  expect_error(prepare_observed_contrasts(ci,tr,0,FALSE,'none','BM'),
               'PSS model fitting failed.*injected optimizer failure')
  fit_models <- old_fit
})
check('prepared candidates are reused by both canonical and FOP selection', {
  result <- prepare_observed_contrasts(ci,tr,0,FALSE,'none','BM')
  old <- lean_candidate_df;lean_candidate_df <- function(...) stop('unexpected rescore')
  fop <- fop_pair_sel.f(result$ctx,ordinal=FALSE,top_pct=0,max_fop=3,seed=1998,
                        candidate_result=result$candidates)
  lean_candidate_df <- old
  stopifnot(nrow(fop$canon_pairs)>0, identical(fop$cand_df,result$candidates$cand_df))
})
check('multi-domain FOP alternatives preserve the generated design and PSS metadata', {
  tree12 <- ape::read.tree(text='(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1,((I:1,J:1):1,(K:1,L:1):1):1);')
  values <- setNames(rep(c(.1,.9,.2,.8),3),tree12$tip.label)
  ci12 <- data.frame(species=names(values),trait=unname(values),n_population=NA_real_,
    n_cases=NA_real_,trait_ci_lb=NA_real_,trait_ci_ub=NA_real_)
  prepared <- prepare_observed_contrasts(ci12,tree12,0,FALSE,'none','BM')
  selected <- fop_pair_sel.f(prepared$ctx,ordinal=FALSE,top_pct=0,max_contrasts=3,
    max_fop=10,seed=1998,candidate_result=prepared$candidates)
  stopifnot(length(selected$hypotheses)>1)
  for(h in selected$hypotheses) {
    exported <- gen_traitfiles(h)
    stopifnot(nrow(h)==3,all(table(exported$pair)==2),!anyDuplicated(exported$species))
    for(i in seq_len(nrow(h)))
      stopifnot(exported$species[exported$pair==h$cluster[i] & exported$contrast_group==1]==h$species1[i])
  }
})
check('provenance rejects changed trait, tree and selection settings', {
  provenance <- contrast_provenance(d,tr,'score','n','c',0,FALSE,'none','BM',1998)
  path <- tempfile(fileext='.rds');saveRDS(list(provenance=provenance),path)
  read_observed_contrasts(path,provenance)
  for (field in c('top_pct','seed','force_model','version')) {
    changed <- provenance;changed[[field]] <- 'different'
    expect_error(read_observed_contrasts(path,changed),'provenance mismatch')
  }
  changed <- provenance;changed$data$score[1] <- 123
  expect_error(read_observed_contrasts(path,changed),'provenance mismatch')
  changed <- provenance;changed$tree$edge.length[1] <- 123
  expect_error(read_observed_contrasts(path,changed),'provenance mismatch')
})
check('separated count intervals cannot reverse phenotype orientation', {
  count_ci <- ci;count_ci$n_population <- 100
  count_ci$trait_ci_lb <- 1-count_ci$trait-.01
  count_ci$trait_ci_ub <- 1-count_ci$trait+.01
  expect_error(prepare_observed_contrasts(count_ci,tr,0,FALSE,'none','BM'),'conflict with phenotype orientation')
})
cat('Passed',checks,'contrast regression groups\n')
