# WGCNA
library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(parallel)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Example: NeuN 
# Load tables 
load("WGCNA_EXONS_NEUN.RData")

tab <- expAdj
datExpr <- as.data.frame(t(tab));
names(datExpr) <- rownames(tab);
rownames(datExpr) <- names(tab);

powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,
  powerVector=powers,
  verbose = 1, 
  blockSize= 14000, 
  networkType = "signed",
  RsquaredCut = 0.85) 
PWR=sft$powerEstimate


TOM = TOMsimilarityFromExpr(datExpr, 
      power= PWR,
      corType = "bicor",
      networkType="signed",
      TOMType="signed",
      TOMDenom = "mean",
      nThreads = 15,
      verbose = 5, 
      indent = 0)

hub = rowSums(TOM)
permTab <- data.frame(Gene = names(datExpr),Original = rowSums(TOM))

## Permutation
P=100  ## select number of permutation resamples
index.b <- list()
Y.b <- list()
TOM.b <- list()
hub.b <- list()
datPerm <- list()
sft <- list()
result=matrix(nrow=ncol(datExpr), ncol=P)
for (i in 1:P)
  {
  set.seed(i*100+1)
  print(i)
  Y.b <- mclapply(1:P,mc.cores =12, function(i){apply(tab, 2, function(col){sample(col)})}) #Randomize 100 times the within_data
  datPerm[[i]] <- as.data.frame(t(Y.b[[i]]))
  TOM.b[[i]] = TOMsimilarityFromExpr(datPerm[[i]],
          power=PWR,
          corType = "bicor",
          networkType="signed",
          TOMType="signed",
          TOMDenom = "mean",
          nThreads = 15,
          verbose = 5, 
          indent = 0)
  hub.b[[i]] = rowSums(TOM.b[[i]])
  result[,i]<-hub.b[[i]]
}

for (i in 1:nrow(permTab)){
permTab$p <- sum(abs(result[i,] >= abs(permTab$Original[i])))/P
}

save(permTab,result,file="PERM_WGCNA.RData")


