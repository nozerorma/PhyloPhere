# Count-interval eligibility remains mandatory when PSS matching is enabled.
suppressPackageStartupMessages(library(ape))
source('subworkflows/CT/local/scripts/lean_contrast_selector.R')
tree <- read.tree(text='((A:1,B:1):1,C:2);')
D <- cophenetic(tree)
trait <- c(A=.8, B=.6, C=.1)
lb <- c(A=.70, B=.45, C=.02)
ub <- c(A=.90, B=.75, C=.20)
calculate_pairwise_scores <- function(...) data.frame(
  Species1=c('A','A','B'), Species2=c('B','C','C'),
  TraitValue1=c(.8,.8,.6), TraitValue2=c(.6,.1,.1),
  PatristicDistance=c(2,4,4), FinalScore=c(1,1.1,3))
result <- evaluate_lean_contrast_selection(trait, D, 1, tree, D, D, 'OU',
    ci_lb=lb, ci_ub=ub, pss_profile=1, pss_tol=.25)
stopifnot(result$tier == 1L, result$mode == 'ci', result$fg == 'A', result$bg == 'C')
# The overlapping A/B pair matches perfectly, but cannot rescue this strict tolerance.
rejected <- evaluate_lean_contrast_selection(trait, D, 1, tree, D, D, 'OU',
    ci_lb=lb, ci_ub=ub, pss_profile=1, pss_tol=.05)
stopifnot(rejected$tier == 0L, grepl('not matched', rejected$reason))
cat('CI eligibility and PSS tolerance verified\n')
