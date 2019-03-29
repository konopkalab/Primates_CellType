rm(list=ls())
# Load libraries
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DMwR))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))

source("Utility_Functions.R")

####################################
# Load Inputs and filter for Neun+ #
####################################
load("OUTPUTS_EXONS_OLIG2/WGCNA_EXONS_OLIG2.RData")

exp <- as.data.frame(expQuant)
pd <- pd

B=100  ## select number of times for leave-multiple-out method
Y.b <- mclapply(1:B,mc.cores =12, function(i){apply(exp, 2, function(col){sample(col)})}) #Randomize 100 times the within_data

for (i in 1:B)
  {
    colnames(Y.b[[i]]) <- colnames(exp)
    rownames(Y.b[[i]]) <- rownames(exp)
    Y.b[[i]] <- as.data.frame(t(Y.b[[i]]))
}

# ANOVA
Model <- vector("list", length = B)
permutation_OLIG2 <- vector("list",length = B)
for (i in 1:B)
    {
        permutation_OLIG2[[i]] <- data.frame(matrix(nrow=ncol(Y.b[[i]]), ncol=1, dimnames = list(colnames(Y.b[[i]]), c("Anova"))))
        {
        for (j in 1:ncol(Y.b[[i]])) 
            {   
            Model[[i]]=tryCatch(anova(lm(as.formula(paste("Y.b[[i]][,j] ~ ", paste(colnames(pd),collapse = " + "))), data = pd)),warning =  function(w) w)
                if (j %% 7560 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    permutation_OLIG2[[i]][j,"Anova"] = Model[[i]]$"Pr(>F)"[1]
                }
            }
        }
    }
save(permutation_OLIG2,file="permutation_OLIG2.RData")


myfiles <- vector("list", length = B)
for (i in 1:B)
    {
      permutation_OLIG2[[i]] <- permutation_OLIG2[[i]] %>% rownames_to_column('Gene')
      myfiles[[i]] <- permutation_OLIG2[[i]][permutation_OLIG2[[i]]$Anova < 0.05,]
    }

df <- data.frame(Perm = sapply(myfiles, nrow))
obs <- read.table("OUTPUTS_EXONS_OLIG2/SPECIES_DGE.txt")

pdf("Permutation_Histogram_OLIG2.pdf",width=4,height=2,useDingbats=FALSE)
ggpubr::gghistogram(df, 
  x = "Perm", 
  fill = "magenta4",
   add = "mean", rug = FALSE,
   bins = 30)+
  ggplot2::geom_vline(xintercept=1000, colour="red",linetype="dashed") 
dev.off()



AnovaP <- do.call(cbind,permutation_OLIG2) %>% as.data.frame()
AnovaP$P_mean <- rowMeans(AnovaP)

obs <- read.table("OUTPUTS_EXONS_OLIG2/SPECIES_DGE.txt")

df <- data.frame(Obs = -log10(obs$Anova_Species), Boot = -log10(AnovaP$P_mean))

pdf("Permutation_Anova_OLIG2.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(df, x = "Boot", y = "Obs",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.x = 0.3, label.sep = "\n")
   )+
xlab("-log10(Permutation Anova)")+ 
ylab("-log10(Observed Anova)")+
theme(legend.position="none")
dev.off()


