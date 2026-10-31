# Weighted RER estimation: a technical error must never change the estimator.
a <- commandArgs(TRUE)
suppressPackageStartupMessages(library(RERconverge))
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
trait <- rer_trait(a[1]); trees <- readRDS(a[2])
min_sp <- rer_integer(if(length(a)>=4L) a[4] else 10,'min.sp',3)
if (length(intersect(names(trait$values),trees$masterTree$tip.label)) < min_sp)
  stop('Insufficient trait/master-tree species overlap')
identity <- rer_tree_contract(trees)
mat <- RERconverge::getAllResiduals(trees,useSpecies=names(trait$values),
  transform='sqrt',use.weights=TRUE,norm='scale',min.sp=min_sp)
rer_ids(rownames(mat),'Matrix genes')
if (!identical(colnames(mat),colnames(trees$paths))) stop('Native RER column order differs from tree paths')
colnames(mat) <- identity$paths
attr(mat,'rer_contract') <- list(version='rer_v1',tree=identity$tree,paths=colnames(mat),
  genes=rownames(mat),species=names(trait$values),use.weights=TRUE,transform='sqrt',norm='scale',min.sp=min_sp,
  package_version=as.character(packageVersion('RERconverge')))
saveRDS(mat,a[3])
