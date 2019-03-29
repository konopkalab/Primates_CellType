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
suppressPackageStartupMessages(library(ggplot2))
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
                  filter(!(Species == 'Hsap' & row_number() == sample(1:n(),1)))  %>% 
                  filter(!(Species == 'PanTro' & row_number() == sample(1:n(),1)))  %>% 
                  filter(!(Species == 'RheMac' & row_number() == sample(1:n(),1)))  %>% 
                  column_to_rownames('Sample') %>% 
                  as.data.frame() %>%
                  droplevels()
  Y.b[[i]]=exp[rownames(exp) %in% rownames(new_pd[[i]]),]
  Y.b[[i]]=Y.b[[i]][match(rownames(new_pd[[i]]),rownames(Y.b[[i]])),]
}


########################
# Anova across SPECIES #
########################

# ANOVA
# ANOVA
Model <- vector("list", length = B)
downsample_OLIG2 <- vector("list",length = B)
for (i in 1:B)
    {
        downsample_OLIG2[[i]] <- data.frame(matrix(nrow=ncol(Y.b[[i]]), ncol=1, dimnames = list(colnames(Y.b[[i]]), c("Anova"))))
        {
        for (j in 1:ncol(Y.b[[i]])) 
            {   
            Model[[i]]=tryCatch(anova(lm(as.formula(paste("Y.b[[i]][,j] ~ ", paste(colnames(new_pd[[i]]),collapse = " + "))), data = new_pd[[i]])),warning =  function(w) w)
                if (j %% 7500 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    downsample_OLIG2[[i]][j,"Anova"] = Model[[i]]$"Pr(>F)"[1]
                }
            }
        }
    }
save(downsample_OLIG2,file="LOO_OLIG2.RData")
tmp <- downsample_OLIG2

myfiles <- vector("list", length = B)
for (i in 1:B)
    {
      tmp[[i]] <- tmp[[i]] %>% rownames_to_column('Gene')
      tmp[[i]]$AnovaAdj <- p.adjust(tmp[[i]]$Anova,method="bonferroni")
      myfiles[[i]] <- tmp[[i]][tmp[[i]]$AnovaAdj < 0.05,]
    }


df <- data.frame(Perm = sapply(myfiles, nrow))
obs <- read.table("OUTPUTS_EXONS_OLIG2/SPECIES_DGE.txt")
intercept <- nrow(obs[obs$AnovaAdj_Species < 0.05,])
pdf("LOO_Histogram_OLIG2.pdf",width=4,height=2,useDingbats=FALSE)
gghistogram(df, 
  x = "Perm", 
  fill = "magenta4",
   add = "mean", rug = FALSE,
   bins = 30)+
  geom_vline(xintercept=intercept, colour="red",linetype="dashed")+
  xlab("")+
  xlim(2000,2600)
  dev.off()


AnovaP <- do.call(cbind,downsample_OLIG2) %>% as.data.frame()
AnovaP$P_mean <- rowMeans(AnovaP)

obs <- read.table("OUTPUTS_EXONS_OLIG2/SPECIES_DGE.txt")

df <- data.frame(Obs = -log10(obs$Anova_Species), Boot = -log10(AnovaP$P_mean))

pdf("LOO_Anova_OLIG2.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(df, x = "Boot", y = "Obs",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
   )+
xlab("-log10(Bootstrap Anova)")+ 
ylab("-log10(Observed Anova)")+
theme(legend.position="none")
dev.off()


