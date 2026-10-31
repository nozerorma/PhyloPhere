# Compare each canonical global matrix to independently aggregated position rows.
a <- commandArgs(TRUE)
for(i in seq_along(a[-1])) {
  b <- readRDS(file.path(a[i+1],'matrices.rds'))
  x <- read.delim(file.path(a[1],sprintf('block_%06d.tsv',i)),colClasses='character',quote='',na.strings=NULL,check.names=FALSE)
  expected <- matrix(0,nrow(b$stat$global),ncol(b$stat$global),dimnames=dimnames(b$stat$global))
  if(nrow(x)) {
    g <- match(x$gene,rownames(expected));c <- match(x$cycle,colnames(expected))
    if(anyNA(g)||anyNA(c))stop('Position gene/cycle missing from matrix')
    expected[cbind(g,c)] <- as.numeric(x$max)
  }
  if(!isTRUE(all.equal(as.numeric(expected),as.numeric(b$stat$global),tolerance=1e-12)))stop('Global null matrix disagrees with positions')
}
