# Phylogenetic modelling from Vincze et al Nature paper from 2022 (Cancer Risk Across Mammals)
# Extrzct this into a different scirpt, create function for creatging quick tables for each of the different analyses
## Similarly as before, perform everything in duplicate 

# Zero inflated phylogenetic model from VinczeEtal2021Nature
ZIlogis <- function(formbin,   # model formula for binomial model
                    formlogis, # model formula for logistic model
                    databin,   # dataset for binomial GLM repsonse variable
                    datalogis, # dataset with logistic response variable
                    phyloTree, # phylogeneitc tree for logistic regression
                    phylobin,  # phylogenetic tree for binomial regression
                    l0) {      # starting value of Pagel's lambda in PGLS regression
  out <- list()
  out$bin <- binaryPGLMM(formbin, data = databin, phy = phylobin)
  out$wlogis <- gls(formlogis,
                    correlation = corPagel(l0, phy = phyloTree, fixed = FALSE, form=~species),
                    data = datalogis, weights = ~1/log(necropsy_count) )
  if(out$wlogis$modelStruct[1] < 0){ # refit model with lambda=0 if lambda converged to negative
    out$wlogis <- gls(formlogis,
                      correlation = corPagel(0, phy = phyloTree, fixed = TRUE, form=~species),
                      data = datalogis, weights = ~1/log(necropsy_count))}
  class(out) <- "ZILogis"
  return(out)}

# Summary function for the zero-inflated phylogenetic model
sumZILogis <- function(x,...){
  x1 <- cbind(as.data.frame(x$bin[[2]]), as.data.frame(x$bin[3]), as.data.frame(x$bin[5]), as.data.frame(x$bin[6]))
  x1[,1:3] <- round(x1[,1:3],2)
  x1[,4] <- round(x1[,4],4)
  x1[,1] <- paste(x1[,1]," (", x1[,2],")", sep='')
  x1 <- x1[,-2]
  x1b <- cbind(row.names(x1),x1[,1:3]);
  row.names(x1b) <- NULL
  x1b[nrow(x1)+1, 4] <- paste("n=", nrow(as.data.frame(x$bin[9])))
  x1b[nrow(x1)+1, 3] <- paste("s2=", round(as.numeric(x$bin[7]),2))
  x1b[nrow(x1)+1, 2] <- paste("P_s2=", round(as.numeric(x$bin[8]),4))
  x1b[,1] <- as.character(x1b[,1])
  x1b[nrow(x1)+1, 1] <- "ModelStats"
  x2w <- as.data.frame(summary(x$wlogis)$tTable); names(x2w) <- names(x1)
  x2w[,1:3] <- round(x2w[,1:3],2)
  x2w[,4] <- round(x2w[,4],4)
  x2w[,1] <- paste(x2w[,1]," (", x2w[,2],")", sep='')
  x2w <- x2w[,-2]
  x2wb <- cbind(row.names(x2w),x2w[,1:3]);
  row.names(x2wb) <- NULL
  x2wb[nrow(x2w)+1, 2] <- paste("AIC=", round(summary(x$wlogis)$AIC,2))
  x2wb[nrow(x2w)+1, 4] <- paste("n=", length(x$wlogis$residuals))
  x2wb[nrow(x2w)+1, 3] <- paste("L=", round(unlist(summary(x$wlogis))$modelStruct.corStruct,2))
  x2wb[,1] <- as.character(x2wb[,1])
  x2wb[nrow(x2w)+1, 1] <- "ModelStats"
  names(x2wb) <- names(x1b)
  names(x1b)[1] <- ""
  x3 <- x1b[1,]; x3[,1:4] <- ""; x3[1,1] <- "Probability of zeros"
  names(x2wb)[1] <- ""
  x5 <- x3; x5[1,1] <- "Weighted logistic reg."
  X <- rbind(x3,x1b,x5,x2wb)
  names(X)[1:4] <- c('_',"b (SE)","z/t-value", "p-value")
  return(X)        }

sumGLS <- function(model){ # summary function GLS
  su <- as.data.frame(summary(model)$tTable)
  su[,1:3] <- round(su[,1:3],2); su[,4] <- round(su[,4],4)
  su[,1] <- paste(su[,1],'_(', su[,2],')', sep='')
  su <- su[,c(1,3,4)]
  x <- nrow(su)
  #su[x+1,1] <- paste("ShaP=",round(as.numeric(unlist(shapiro.test(resid(model)))[[2]][1]),4), sep="")
  su[x+1,2] <- paste("AIC=", round(AIC(model),2), sep="")
  su[x+1,3] <- paste("n=", length(resid(model)), sep='')
  su[x+1,1] <-paste('lambda=',round(as.numeric(model[[1]])[1], 2), sep='') # put lambda value in bottom right corner
  row.names(su)[nrow(su)] <- "ModelStats"
  return(su)
}