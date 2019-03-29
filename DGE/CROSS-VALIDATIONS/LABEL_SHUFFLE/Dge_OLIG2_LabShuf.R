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
suppressPackageStartupMessages(library(ggpubr))
source("Utility_Functions.R")

####################################
# Load Inputs and filter for Neun+ #
####################################
load("OUTPUTS_EXONS_OLIG2/WGCNA_EXONS_OLIG2.RData")

exp <- as.data.frame(t(expQuant))
pd <- pd

mod <- model.matrix(~Species+Sex+HumAge+RIN, data =pd)
mod0 <- model.matrix(~Sex+HumAge+RIN, data = pd)
svaobj <- sva(as.matrix(expQuant),mod,mod0,n.sv=NULL,B=100,method="two-step")
svs <- svaobj$sv
colnames(svs) <- paste("SV",1:ncol(svs),sep="")
pd <- cbind(pd,svs)

B=100  ## select number of times for leave-multiple-out method
new_pd <- vector("list", length = B)
Y.b <- vector("list",length = B)
for (i in 1:B)
  {
  new_pd[[i]] <- pd %>% # Subset the tables by 10 samples per each species and creating a 100 leave-multiple-out table
                  rownames_to_column('Sample') %>% 
                  group_by(Species) %>% 
                  mutate( Sample = sample(Sample)) %>% 
                  column_to_rownames('Sample') %>% 
                  as.data.frame() %>%
                  droplevels()
  Y.b[[i]]=exp[rownames(exp) %in% rownames(new_pd[[i]]),]
  Y.b[[i]]=Y.b[[i]][match(rownames(new_pd[[i]]),rownames(Y.b[[i]])),]
}
save(Y.b,new_pd,file = "Labels_Shuffled_OLIG2_Data.RData")

# ANOVA
Model <- vector("list", length = B)
LabShuf_OLIG2 <- vector("list",length = B)
for (i in 1:B)
    {
        LabShuf_OLIG2[[i]] <- data.frame(matrix(nrow=ncol(Y.b[[i]]), ncol=1, dimnames = list(colnames(Y.b[[i]]), c("Anova"))))
        {
        for (j in 1:ncol(Y.b[[i]])) 
            {   
            Model[[i]]=tryCatch(anova(lm(as.formula(paste("Y.b[[i]][,j] ~ ", paste(colnames(pd),collapse = " + "))), data = pd)),warning =  function(w) w)
                if (j %% 7560 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    LabShuf_OLIG2[[i]][j,"Anova"] = Model[[i]]$"Pr(>F)"[1]
                }
            }
        }
    }
save(LabShuf_OLIG2,file="LabShuf_OLIG2.RData")


myfiles <- vector("list", length = B)
for (i in 1:B)
    {
      LabShuf_OLIG2[[i]] <- LabShuf_OLIG2[[i]] %>% rownames_to_column('Gene')
      myfiles[[i]] <- LabShuf_OLIG2[[i]][LabShuf_OLIG2[[i]]$Anova < 0.05,]
    }

df <- data.frame(Perm = sapply(myfiles, nrow))
obs <- read.table("OUTPUTS_EXONS_OLIG2/SPECIES_DGE.txt")
intercept <- nrow(obs[obs$AnovaAdj_Species < 0.05,])

pdf("LabShuf_Histogram_OLIG2.pdf",width=4,height=2,useDingbats=FALSE)
gghistogram(df, 
  x = "Perm", 
  fill = "magenta4",
  add = "mean", rug = FALSE,
  bins = 30)+
  geom_vline(xintercept=intercept, colour="red",linetype="dashed") +
  ggtitle("Label Shuffle OLIG2")
dev.off()
