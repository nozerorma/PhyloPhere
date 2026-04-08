########################
###function cluster#####
#####permulations#######
########################
library(RERconverge)
library(ape)
library(RRPP)
library(dplyr)
library("tidyr")
library("phylolm")
library(caper)

#READ ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

TABLE_AA <- read.csv(args[1], header=F, sep = "\t")
print(colnames(TABLE_AA))
colnames(TABLE_AA) <- c("ID", "spc", "Pos", "AA", "Top", "Trait")
TRAIT <- read.csv(args[2], header=F, sep = "\t")
print(colnames(TRAIT))
colnames(TRAIT) <- c("spc", "Value", "Trait")
TREE <- read.newick(args[3])
ALL_DF <- left_join(TABLE_AA, TRAIT, by=c("spc" = "spc", "Trait" = "Trait")) %>%  drop_na()
OUTPUT_TABLE <- args[4]
print(ALL_DF)

#SUBSET TRAITS

IDs <- c()
traits <- c()
pvalues <- c()
slopes <- c()
positions <- c()
#print(ALL_DF)
list_df <- split(ALL_DF, list(ALL_DF$Trait, ALL_DF$Pos), drop=TRUE)
#print(list_df)
#ESTIMATE RRPP IF CONTINUOUS
for (i in 1:length(list_df)) {
 DF <- list_df[[i]]
#print(DF)
#print(DF)
 # assign(new_names[i], iris_split[[i]])
#}
rownames(DF) <- DF$spc
#print(length(rownames(ALL_DF)))
#print(length(TREE$tip.label))
labels <- match(rownames(DF), TREE$tip.label)
labels_good <- labels[!is.na(labels)]
#print(labels_good)
#print(TREE$tip.label[-labels_good])
pruned_tree<-drop.tip(TREE,TREE$tip.label[-labels_good])
#print(pruned_tree)
corr <- vcv.phylo(pruned_tree)
#print(length(rownames(corr)))
#print(length(rownames(DF)))
#print(corr)

#print(rownames(corr))
#print(setdiff(rownames(corr), rownames(DF)))
#print(setdiff(rownames(DF), rownames(corr)))

SORTED_DF <- DF[match(rownames(corr), rownames(DF)), ]
#print(row.names(SORTED_DF))
#print(row.names(corr))

#print(unique(SORTED_DF$Value))

if (length(unique(SORTED_DF$Value)) > 2){
#ESTIMATE RRPPm CONTINUOUS
#print(ALL_DF)

primate_data <- comparative.data(phy = pruned_tree, data = SORTED_DF, names.col = "spc", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
 #plot(trait~gene_pred)
 #abline(a = coef(fit)[1], b = coef(fit)[2])
  if (class(try(pgls(as.numeric(as.character(Value))~as.numeric(as.character(Top)),
                     data = primate_data, lambda="ML"))) == "try-error"){
    lambda <- NA
    pvalue <- NA
    r_sq <- NA
    slope <- NA
  }
  else{
     fit<-pgls(as.numeric(as.character(Value))~as.numeric(as.character(Top)),
             data= primate_data, lambda="ML")
    #plot(as.numeric(as.character(X$trait))~as.numeric(as.character(X$gene)))
    #abline(a = coef(fit)[1], b = coef(fit)[2])
    res<- residuals(fit, phylo = TRUE)
    res<- res/sqrt(var(res))[1]
    #check standardized outliers
    rownames(SORTED_DF) <- SORTED_DF$spc 
    X_nooutliers<-SORTED_DF[which(abs(res)<3),]
    primate_data <- comparative.data(phy = pruned_tree, data = X_nooutliers, names.col = "spc", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
    if (class(try(pgls(as.numeric(as.character(Value))~as.numeric(as.character(Top)),
                       data = primate_data, lambda="ML"))) == "try-error"){
      lambda <- NA
      pvalue <- NA
      r_sq <- NA
      slope <- NA
    }
 else{

      fit_nooutliers <- pgls(as.numeric(as.character(Value))~as.numeric(as.character(Top)),
                  data=primate_data, lambda="ML")
      slope <- summary(fit_nooutliers)[5][[1]][2]
      r_sq <- summary(fit_nooutliers)[10]
      lambda <- summary(fit_nooutliers)[6][[1]][2]
      pvalue <- summary(fit_nooutliers)[5][[1]][2,4]
      print(summary(fit_nooutliers)[5][[1]])
      print(slope)
}}
p_val <- pvalue
slope_curr <- slope
print(p_val)
print(slope_curr)
#print(SORTED_DF)
} 
else {

 if (class(try(phyloglm(as.numeric(as.character(Value)) ~ as.numeric(as.character(Top)), data = SORTED_DF,
pruned_tree, btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot = 0, full.matrix = TRUE))) == "try-error")
{
    p_val <- NA
    slope_curr <- NA
  }
  else{

binaryGLS <- phyloglm(as.numeric(as.character(Value)) ~ as.numeric(as.character(Top)), data = SORTED_DF, 
pruned_tree, btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot = 0, full.matrix = TRUE)

#print(summary(binaryGLS))
#print(summary(binaryGLS)[2][[1]][1])
#print(summary(binaryGLS)[2][[1]][2])
#print(summary(binaryGLS)[2][[1]][3])
#print(summary(binaryGLS)[2][[1]][4])
print(summary(binaryGLS)[2][[1]])
#print(summary(binaryGLS)[2][[1]][2])

p_val <- summary(binaryGLS)[2][[1]][8]
slope_curr <-  summary(binaryGLS)[2][[1]][2]
#print(p_val)
print(summary(binaryGLS)[2][[1]])
print(slope_curr)
}}
pvalues <- c(pvalues, p_val)
slopes <- c(slopes, slope_curr)
traits <- c(traits, unique(as.character(SORTED_DF$Trait))[1]) 
positions <- c(positions, unique(as.character(SORTED_DF$Pos))[1])
IDs <- c(IDs, unique(as.character(SORTED_DF$ID))[1])
}

output_df <- as.data.frame(cbind(traits, positions, pvalues, slopes, IDs))
write.table(output_df, OUTPUT_TABLE, row.names=F, col.names=F, sep="\t", quote=F)
