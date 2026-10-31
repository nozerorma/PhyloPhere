# Optional saved-data CLI smoke (requires a verified rer_v1 matrix from build_matrix): Rscript tests/test_wp6_rer_real.R <marine_complete>
suppressPackageStartupMessages(library(RERconverge))
source('subworkflows/RERCONVERGE/local/rer_permulations.R')
root<-normalizePath(commandArgs(TRUE)[1]); wd<-tempfile('rer_cli_');dir.create(wd)
rs<-file.path(R.home('bin'),'Rscript')
trees_path<-file.path(root,'rerconverge/rer_objects/gene_trees.nwk.masterTree.output')
matrix_path<-file.path(root,'rerconverge/rer_objects/marine.RERmatrix.output')
trees<-readRDS(trees_path); mat<-readRDS(matrix_path)
run_cli<-function(kind,trait,extra) {
 out<-file.path(wd,paste0(kind,'.output'))
 cmd<-c(normalizePath(paste0('subworkflows/RERCONVERGE/local/',kind,'_rer.R')),trait,trees_path,
         file.path(wd,paste0(kind,'.paths')),matrix_path,out,extra)
 log<-system2(rs,cmd,stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(log,'status')))stop(paste(log,collapse='\n'))
 stopifnot(dir.exists(sub('[.]output$','.perms.rds',out)),length(list.files(sub('[.]output$','.perms.rds',out),pattern='^blk_.*[.]rds$'))==2)
 res<-readRDS(out);stopifnot(attr(res,'n_perms')==4,attr(res,'rer_batch_seed')==79)
 res
}
bin<-run_cli('binary',file.path(root,'rerconverge/rer_traits/marine.polished.output'),c(15,2,3,'terminal',2,2,79))
trait_env<-new.env();load(file.path(root,'rerconverge/rer_traits/marine.polished.output'),envir=trait_env)
fg<-names(trait_env$trait_vector)[!is.na(trait_env$trait_vector)&trait_env$trait_vector==1]
paths<-foreground2Paths(fg,trees,clade='terminal')
ref<-rer_score_paths(mat,paths,'binary',15,min.pos=2,winR=3)
stopifnot(isTRUE(all.equal(bin$P,ref$P,tolerance=1e-12)),isTRUE(all.equal(bin$Rho,ref$Rho,tolerance=1e-12)))
trait_vector<-setNames(seq_along(trees$masterTree$tip.label),trees$masterTree$tip.label)
n_vector<-c_vector<-NULL; trait<-file.path(wd,'continuous_trait.RData');save(trait_vector,n_vector,c_vector,file=trait)
cont<-run_cli('continuous',trait,c(15,1,2,2,2,'cc','none',79))
ref<-rer_score_paths(mat,char2Paths(trait_vector,trees),'continuous',15,winR=1,winT=2)
stopifnot(isTRUE(all.equal(cont$P,ref$P,tolerance=1e-12)),isTRUE(all.equal(cont$Rho,ref$Rho,tolerance=1e-12)))
cat('PASS real RER CLI binary/continuous, non-default harmonized options, reduced null blocks and gene p\n')
unlink(wd,recursive=TRUE)
