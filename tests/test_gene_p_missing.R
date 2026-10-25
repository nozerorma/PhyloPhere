# Missing evidence must stay unavailable; measured absence (zero) remains a value.
source('subworkflows/ENRICHMENT/local/src/gene_p.R')
S <- rbind(empty=c(NA,NA,NA),missing_obs=c(.2,.4,.8),measured_zero=c(0,0,0),partial=c(.1,NA,.5))
o <- c(.7,NA,0,.5)
g <- gene_p(S,o)
stopifnot(is.na(g$p_obs[1]),is.na(g$p_obs[2]),g$p_obs[3]==1,abs(g$p_obs[4]-2/3)<1e-15,
          all(is.na(g$p_null[1,])),all(is.finite(g$p_null[2,])),all(g$p_null[3,]==1),is.na(g$p_null[4,2]))
# Independent count oracle, including observed/null ties within tolerance.
S2<-rbind(c(0,.4,.4+2e-16,.9),c(.1,NA,.5,.8));o2<-c(.4,.7);z<-gene_p(S2,o2)
for(i in seq_len(nrow(S2))) {
 x<-S2[i,];ok<-!is.na(x)
 stopifnot(abs(z$p_obs[i]-(1+sum(x>=o2[i]-1e-12,na.rm=TRUE))/(1+sum(ok)))<1e-15)
 for(j in which(ok))stopifnot(abs(z$p_null[i,j]-(1+sum(x[-j]>=x[j]-1e-12,na.rm=TRUE))/sum(ok))<1e-15)
}
cat('PASS: empty null, missing observed, measured zero, partial NA, ties, leave-one-out oracle\n')
