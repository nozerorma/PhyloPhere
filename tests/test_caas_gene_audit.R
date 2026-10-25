# Independent explicit pools for the experimental gene-score decomposition.
source('bin/audit_caas_gene_evidence.R')
S<-rbind(c(0,.2,.5,0),c(.1,.5,.5,.9),c(0,0,0,0),c(.2,.2,.7,0))
f<-fit_fact_gene(S)
stopifnot(identical(as.numeric(f$breaks),c(2,3)),identical(as.numeric(lengths(f$pools)),c(2,3,4)))
p<-fact_gene_p(c(.5,.9,0,.7),f$nd,f)
stopifnot(max(abs(p-c(.3,.4,1,.32)))<1e-14)
stopifnot(is.na(fact_gene_p(NA,2,f)),abs(fact_gene_p(.6,0,f)-1/15)<1e-14)
# Empty conditional pool has survival one; non-detection observed has p one.
e<-fit_fact_gene(matrix(0,2,4));stopifnot(fact_gene_p(.3,0,e)==.2,fact_gene_p(0,0,e)==1)
cat('PASS experimental gene audit: weighted classes, nd+1, inclusive tails, NA, empty pool, measured zero\n')
