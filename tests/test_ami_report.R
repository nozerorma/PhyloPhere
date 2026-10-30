# Rscript tests/test_ami_report.R <repo root> <scratch directory>
# Edge cases of 12.AMI_analysis.Rmd against the minimal DOMINO fixture
# (tests/fixtures/ami): six nominal seeds, two BH hits, two network connectors;
# then an empty nominal list, which must skip cleanly instead of failing.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1])
work <- args[2]
fixture <- file.path(root,'tests','fixtures','ami')
render_ami <- function(wd, global_list, rich=FALSE, max_blocks=NULL) {
  dir.create(wd,recursive=TRUE,showWarnings=FALSE)
  writeLines(global_list,file.path(wd,'global.txt'))
  for (f in c('network.sif','background.txt','gene_scores.tsv')) file.copy(file.path(fixture,f),wd,overwrite=TRUE)
  file.copy(file.path(fixture,'modules'),wd,recursive=TRUE,overwrite=TRUE)
  file.copy(file.path(fixture,'cache'),wd,recursive=TRUE,overwrite=TRUE)
  if(rich) {
    # 40 background genes and four overlapping terms, so that the module (g1..g8) has enriched terms.
    ids <- paste0('9606.p',1:40); sy <- paste0('g',1:40)
    wgz <- function(lines,f){con<-gzfile(file.path(wd,'cache',f),'w');writeLines(lines,con);close(con)}
    wgz(c('#string_protein_id\tpreferred_name\tprotein_size\tannotation',paste(ids,sy,100,'fixture',sep='\t')),'9606.protein.info.v12.0.txt.gz')
    wgz(c('#string_protein_id\talias\tsource',paste(ids,sy,'fixture',sep='\t')),'9606.protein.aliases.v12.0.txt.gz')
    term <- function(gs,id,name) paste(paste0('9606.p',gs),'Process',id,name,sep='\t')
    wgz(c(term(1:6,'GO:1','Alpha term'),term(1:5,'GO:2','Beta term'),term(1:3,'GO:3','Gamma term'),
          term(7:8,'GO:4','Delta term'),term(20:30,'GO:5','Noise term')),'9606.protein.enrichment.terms.v12.0.txt.gz')
    writeLines(sy,file.path(wd,'background.txt'))
  }
  file.copy(file.path(root,'subworkflows/ENRICHMENT/local/12.AMI_analysis.Rmd'),wd,overwrite=TRUE)
  dir.create(file.path(wd,'src'),showWarnings=FALSE)
  file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src',c('ami_evidence.R','report_presentation.R','term_blocks.R','run_domino_modules.py')),file.path(wd,'src'),overwrite=TRUE)
  g<-read.delim(file.path(wd,'gene_scores.tsv'))
  for(sf in c('','_top','_bottom')) {
    g[[paste0('caas_p.gen',sf)]]<-c(rep(.01,6),rep(1,nrow(g)-6))
    g[[paste0('caas_p.gen.adj',sf)]]<-c(.02,.03,rep(.2,nrow(g)-2))
  }
  for(sf in c('','_top','_bottom')) {g[[paste0('caas_nominal',sf)]]<-g[[paste0('caas_p.gen',sf)]]<.05;g[[paste0('caas_adjusted',sf)]]<-g[[paste0('caas_p.gen.adj',sf)]]<.05}
  write.table(g,file.path(wd,'gene_scores.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  jsonlite::write_json(list(gene_stat='max',nominal_threshold=.05,fdr_threshold=.05),file.path(wd,'metadata.json'),auto_unbox=TRUE)
  stamp<-paste0(paste(sort(unique(global_list)),collapse='\n'),if(length(global_list))'\n' else '')
  dir.create(file.path(wd,'modules/inputs'),showWarnings=FALSE)
  writeLines('g1 g2 g3 g4 g5 g6 g7 g8',file.path(wd,'modules/inputs/slices.txt'))
  sha<-function(f)digest::digest(file=f,algo='sha256')
  jsonlite::write_json(list(contract_version=2L,p_value_scale='raw_upper_tail',p_adjustment='Bonferroni_n_putative',module_threshold=.05,slice_threshold=.3,n_putative=3L,final_modules=1L,
    active_genes_sha256=digest::digest(stamp,algo='sha256',serialize=FALSE),network_sha256=sha(file.path(wd,'network.sif')),
    producer_sha256=sha(file.path(wd,'src/run_domino_modules.py')),slices_sha256=sha(file.path(wd,'modules/inputs/slices.txt')),
    modules_sha256=sha(file.path(wd,'modules/global_domino_modules.tsv')),stats_sha256=sha(file.path(wd,'modules/global_domino_module_stats.tsv')),
    background_sha256=digest::digest(paste0(paste(sort(unique(trimws(readLines(file.path(wd,'background.txt'))))),collapse='\n'),'\n'),algo='sha256',serialize=FALSE)),
    file.path(wd,'modules/global_domino_provenance.json'),auto_unbox=TRUE)
  owd <- setwd(wd); on.exit(setwd(owd))
  pl <- list(background_file='background.txt',background_basename='background.txt',
      project_name='fixture',species=9606,domino_network_score_thr=700,
      gene_scores_file='gene_scores.tsv',gene_metadata_file='metadata.json',string_db_dir=NULL,string_cache_dir='cache',
      domino_network_sif='network.sif',domino_modules_dir='modules',
      domino_edge_scores_file=NULL,seed='1998')
  if(!is.null(max_blocks)) pl$ami_max_blocks <- max_blocks
  rmarkdown::render('12.AMI_analysis.Rmd',quiet=TRUE,params=pl,
    output_file='12.AMI_analysis_fixture.html',envir=new.env())
}

# Six nominal seeds (g1..g6), of which g1/g2 also pass BH; g7/g8 are connectors only.
render_ami(file.path(work,'nominal'),c('g1','g2','g3','g4','g5','g6'))
ev <- read.delim(file.path(work,'nominal','ami_networks','caas_global_gene_evidence.tsv'))
stopifnot(nrow(ev)==6L, all(ev$caas_nominal), sum(ev$caas_adjusted)==2L,
          !any(c('g7','g8') %in% ev$Gene), all(c('g1','g2') %in% ev$Gene[ev$caas_adjusted]))
# The certificate carries the hash of the background of the network: the matching background is accepted by the render above,
# and another one is refused with its reason.
source(file.path(root,'subworkflows/ENRICHMENT/local/src/ami_evidence.R'))
nom <- file.path(work,'nominal'); other <- file.path(nom,'other_background.txt'); writeLines(c('g1','g2'),other)
chk <- ami_module_provenance(file.path(nom,'modules/global_domino_provenance.json'),file.path(nom,'global.txt'),file.path(nom,'network.sif'),other)
stopifnot(!chk$available, grepl('background differs',chk$reason))
# g1/g2 use custom corroboration shapes. Their BH borders must be visible in
# the rendered image, rather than merely present in the exported evidence table.
img <- png::readPNG(file.path(work,'nominal','ami_networks','caas_global_domino_network.png'))
gold <- abs(img[,,1] - 217/255) < .04 & abs(img[,,2] - 154/255) < .04 & img[,,3] < .08
stopifnot(sum(gold) > 500L)
desc <- read.delim(file.path(work,'nominal','ami_networks','module_descriptions_all_lists.tsv'))
stopifnot(nrow(desc)==1L, desc$n_genes==6L, desc$n_nodes==8L, !any(grepl('g7|g8',desc$genes)),
          file.exists(file.path(work,'nominal','12.AMI_analysis_fixture.html')))

# Enriched module: term blocks. Alpha, Beta and Gamma annotate nested sets of the module genes
# (pairwise Jaccard 0.5 to 0.83) and form one block represented by the lowest-FDR term; Delta
# annotates other genes. ami_max_blocks only limits the blocks shown in the label and the summary.
for(mb in list(NULL,1L)) {
  wd <- file.path(work,paste0('rich',if(is.null(mb)) '' else mb))
  render_ami(wd,c('g1','g2','g3','g4','g5','g6'),rich=TRUE,max_blocks=mb)
  tt <- read.delim(file.path(wd,'ami_networks','ami_term_threshold_membership_all_tools.tsv'))
  stopifnot(nrow(tt)>=3L, !'term_raw' %in% names(tt),
            all(c('category','term_id','fdr','block','is_representative','block_size','block_min_jaccard') %in% names(tt)),
            !'Noise term' %in% tt$term,
            sum(tt$is_representative)==length(unique(tt$block)),
            tt$is_representative[which.min(tt$fdr)])
  if(is.null(mb)) {
    stopifnot(length(unique(tt$block))==2L,
              identical(tt$block[tt$term %in% c('Alpha term','Beta term','Gamma term')],rep(1L,3)),
              tt$term[tt$is_representative & tt$block==1L]=='Alpha term',
              all(tt$block_min_jaccard[!is.na(tt$block_min_jaccard)]>=0.2-1e-9))
    d <- read.delim(file.path(wd,'ami_networks','module_descriptions_all_lists.tsv'))
    stopifnot(identical(d$description,'Alpha term / Delta term'), all(tt$shown))
  } else {
    d <- read.delim(file.path(wd,'ami_networks','module_descriptions_all_lists.tsv'))
    stopifnot(length(unique(tt$block))==2L, identical(d$description,'Alpha term'),
              identical(tt$shown,tt$block<=1L))
    bs <- read.delim(file.path(wd,'ami_networks','term_block_summary.tsv'))
    stopifnot(bs$n_blocks==2L, bs$blocks_shown==1L, bs$module_genes_covered_by_shown<1)
  }
  stopifnot(file.exists(file.path(wd,'ami_networks','term_block_summary.tsv')))
}

# Empty nominal list: the report renders and skips the CAAS section instead of failing.
render_ami(file.path(work,'empty'),character())
stopifnot(file.exists(file.path(work,'empty','12.AMI_analysis_fixture.html')))
html<-paste(readLines(file.path(work,'nominal','12.AMI_analysis_fixture.html'),warn=FALSE),collapse='\n')
stopifnot(!grepl('<h4>Modules</h4>',html,fixed=TRUE),!grepl('<h5>global</h5>',html,fixed=TRUE))
cat('PASS AMI report: nominal seeds vs BH hits vs connectors; empty list skips cleanly\n')
