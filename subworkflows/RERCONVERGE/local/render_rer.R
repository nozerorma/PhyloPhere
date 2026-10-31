a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
top <- rer_integer(a[3],'top_n_labels')
rmarkdown::render('5.RERconverge_report.Rmd',params=list(continuous_rds=a[1],traitname=a[2],
  top_n_labels=top,output_dir=''),output_file='5.RERconverge_report.html',quiet=TRUE)
