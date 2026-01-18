
# VALIDATION
## Phylowilcoxtest
#Phylogenetic ANOVA, often associated with the R package phylolm, is a method used to test the association between a continuous trait and one or more predictors, taking into account the phylogenetic relationships among the observations. It's an adaptation of the classic ANOVA that incorporates a phylogenetic tree to control for non-independence among species due to their shared evolutionary history.

#It can be used for hypothesis testing in a phylogenetic context and is a way to perform phylogenetic ANOVA when you need to validate a discovered attribute (like amino acids - AAs in our case) within a phylogenetic tree.# Define functions 
# Implement a phylogenetic t-test-like procedure with empirical p-value calculation
#----------------------------------------------------------------------------------
# Implement a phylogenetic U-Mann-Whitney-like procedure with empirical p-value calculation
#Null hypothesis: two samples come from the same population (i.e. have the same median) or, alternatively, observations in one sample tend to be larger than observations in the other
#-------------------------------------------------------------------------------

phylwilcoxtest <- function(tree, x, y, nsim=10000, data, sp_col) {
  # Ensure that tree is a phylo object
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".")
  }
  
  # Ensure that x and y are column names in data and not separate vectors
  if (!(x %in% names(data)) || !(y %in% names(data))) {
    stop("x and y must be column names in data.")
  }
  
  # Extract x and y values from data, ordered by tree$tip.label
  x_values <- data[[x]]
  y_values <- data[[y]]
  
  # Convert x to factor
  x_values <- as.factor(x_values)
  # Check if factor has exactly 2 levels
  if (length(levels(x_values)) != 2) {
    stop("Factor does not have 2 levels.")
  }
  
  # Compute trait rate for y
  sig2 <- mean(pic(y_values, multi2di(tree, random = FALSE))^2)
  
  # Calculate the observed difference
  result <- wilcox.test(y_values ~ x_values, alternative = "greater", conf.int = TRUE, conf.level = 0.95, digits.rank = 7, exact = TRUE, correct = TRUE)
  
  W.obs <- result$statistic
  
  # Simulate
  sims <- fastBM(tree, sig2 = sig2, nsim = (nsim - 1), bounds = c(0,1))
  W.null <- vector()
  W.null[1] <- W.obs
  
  for (i in 2:nsim) {
    # Recalculate group medians for the simulated data
    simresult <- wilcox.test(sims[,i-1] ~ x_values, alternative = "greater", conf.int = TRUE, conf.level = 0.95, digits.rank = 7, exact = TRUE, correct = TRUE)
    W.sim <- simresult$statistic
    
    W.null[i] <- W.sim
  }
  
  # Calculate p-value
  P.W <- sum(W.null >= W.obs) / nsim
  
  # Calculate phylogenetic signal (Pagel's Lambda)
  names(y_values) <- data[[sp_col]]
  lambda <- phylosig(tree, y_values, method="lambda", test = TRUE)
  
  # Return result
  obj <- list(Wobs = W.obs, Pval = P.W, Lambda = lambda$lambda, LambdaPval = lambda$P, LambdaSig2 = lambda$sig2, LambdaLogL = lambda$logL, LambdaLogL0 = lambda$logL0)
  class(obj) <- "phylt-test"
  # obj
  obj
}

## Define function to find outliers
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}