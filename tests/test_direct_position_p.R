source('subworkflows/SCORING/local/src/null_position_p.R')
obs<-data.frame(Gene=c('g','g','h'),Position=c(0,1,2),.obs=c(.7,.9,0))
ctr<-data.frame(Gene=c('h','g','g'),Position=c(2,1,0),k_emp=c(0,0,2),nd=c(1,0,3),c1=c(999,999,999))
x<-null_position_p(obs,ctr,NULL,5,9)
p<-c(.3,.1,1);expected<-p.adjust(c(p,rep(1,3)),'BH')[1:3]
stopifnot(isTRUE(all.equal(x$pos$p.emp,p)),isTRUE(all.equal(x$pos$p.adj_bh,expected)),x$n_family==6,!any(grepl('fact',names(x$pos))))
# Legacy fit/conditional counts cannot affect direct evidence.
y<-null_position_p(obs,transform(ctr,c1=-100),data.frame(garbage=NA),5,9)
stopifnot(identical(x,y))
empty<-null_position_p(obs[FALSE,],ctr[FALSE,],NULL,5,9)
stopifnot(nrow(empty$pos)==0,empty$n_family==5)
# Guarded coordinate mismatch remains unavailable, not a discovery.
z<-null_position_p(data.frame(Gene=letters[1:10],Position=0,.obs=1),data.frame(Gene=letters[1:10],Position=0,k_emp=0,nd=0),NULL,0,9)
stopifnot(all(is.na(z$pos$p.emp)),all(is.na(z$pos$p.adj_bh)))
cat('PASS direct position evidence: BH family, zero, empty, unseen, legacy counts, coordinate guard\n')
