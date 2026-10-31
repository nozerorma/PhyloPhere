# Validate and prune gene trees using the same species contract as the trait.
a <- commandArgs(TRUE)
suppressPackageStartupMessages({library(ape); library(RERconverge)})
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
trait_col <- if(length(a)>=7L) a[7] else NULL
x <- rer_read_traits(a[2], a[3], trait_col)
species <- x[[a[3]]]
if (!length(species)) stop('No usable trait species')
records <- readLines(a[1],warn=FALSE)
records <- records[nzchar(trimws(records))]
named <- grepl('^[^(]*\t',records)
if (any(named)) {
  if (!all(named)) stop('Mixed named and unnamed gene tree records')
  delimiter <- regexpr('\t',records,fixed=TRUE)
  ids <- substring(records,1L,delimiter-1L)
  rer_ids(ids,'Tree genes')
  newicks <- substring(records,delimiter+1L)
  trees <- ape::read.tree(text=newicks,keep.multi=TRUE)
  if (length(trees)!=length(ids)) stop('Expected one Newick tree per named record')
  names(trees) <- ids
} else trees <- ape::read.tree(a[1],keep.multi=TRUE)
if (!inherits(trees,'multiPhylo') || !length(trees)) stop('No gene trees')
ids <- names(trees)
if (is.null(ids)) ids <- rep('',length(trees))
if (all(!nzchar(ids))) ids <- paste0('gene',seq_along(trees))
rer_ids(ids,'Tree genes')
names(trees) <- ids
mapping <- NULL
tax_path <- if(length(a)>=6L) a[6] else ''
if (nzchar(tax_path) && basename(tax_path) != 'NO_FILE') {
  m <- rer_table(tax_path)
  if (!all(c('tax_id','species') %in% names(m))) stop('TaxID table requires tax_id and species')
  m$species <- rer_species(m$species)
  m <- unique(m[,c('tax_id','species')])
  rer_ids(m$species,'TaxID species')
  if (any(!nzchar(m$tax_id))) stop('Missing tax_id')
  target <- m[m$species %in% species,,drop=FALSE]
  rer_ids(target$tax_id,'Trait tax_id (ambiguous canonical species)')
  canonical <- setNames(target$species,target$tax_id)
  mapping <- setNames(unname(canonical[m$tax_id]),m$species)
}
audit <- data.frame(gene=ids,status='complete',reason='',input_tips=integer(length(trees)),retained_tips=integer(length(trees)))
kept <- list()
for (i in seq_along(trees)) {
  t <- trees[[i]]
  t$tip.label <- rer_species(t$tip.label)
  rer_ids(t$tip.label,paste('Tree',ids[i],'tips'))
  if (is.null(t$edge.length) || any(!is.finite(t$edge.length) | t$edge.length < 0))
    stop('Tree ',ids[i],': finite nonnegative branch lengths required')
  audit$input_tips[i] <- length(t$tip.label)
  if (!is.null(mapping)) {
    mapped <- unname(mapping[t$tip.label]); hit <- !is.na(mapped)
    t$tip.label[hit] <- mapped[hit]
    rer_ids(t$tip.label,paste('Tree',ids[i],'tips after TaxID translation'))
  }
  tips <- intersect(t$tip.label,species)
  audit$retained_tips[i] <- length(tips)
  if (length(tips)<3L) {
    audit$status[i] <- 'excluded'; audit$reason[i] <- 'fewer_than_three_trait_taxa'
  } else kept[[ids[i]]] <- ape::keep.tip(t,tips)
}
rer_write_tsv(audit,paste0(a[4],'.audit.tsv'))
if (!length(kept)) stop('No usable gene trees after pruning; see tree audit')
# The native reader needs at least one tree spanning the union of taxa.
all_sp <- unique(unlist(lapply(kept,`[[`,'tip.label')))
if (!any(vapply(kept,function(t) setequal(t$tip.label,all_sp),logical(1))))
  stop('No gene tree spans the retained taxon union; cannot infer a master tree')
lines <- vapply(names(kept),function(id) paste(id,ape::write.tree(kept[[id]]),sep='\t'),'')
writeLines(lines,a[4])
# Native readTrees uses scan() with quote handling. Opaque IDs prevent literal
# apostrophes/quotes in gene IDs from changing its parser; restore every gene axis.
native_ids <- sprintf('rer_gene_%09d',seq_along(kept))
native_path <- tempfile('rer_native_',tmpdir=getwd())
writeLines(vapply(seq_along(kept),function(i) paste(native_ids[i],ape::write.tree(kept[[i]]),sep='\t'),''),native_path)
obj <- RERconverge::readTrees(native_path)
unlink(native_path)
if (!identical(names(obj$trees),native_ids)) stop('Native reader changed gene coverage')
names(obj$trees) <- names(kept)
rownames(obj$report) <- names(kept)
if (!is.null(rownames(obj$paths))) rownames(obj$paths) <- names(kept)
attr(obj,'rer_provenance') <- list(version='rer_v1',input_sha256=rer_file_hash(a[1]),
  trait_sha256=rer_file_hash(a[2]),taxid_sha256=if(nzchar(tax_path) && basename(tax_path)!='NO_FILE') rer_file_hash(tax_path) else NULL)
saveRDS(obj,a[5])
