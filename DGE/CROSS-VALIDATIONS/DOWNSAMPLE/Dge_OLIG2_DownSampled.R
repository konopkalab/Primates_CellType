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
suppressPackageStartupMessages(library(purrr))
source("Utility_Functions.R")

# Create Output directories
folder_names <- c("OUTPUT_DOWNSAMPLE_OLIG2")
sapply(folder_names, dir.create)

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
                  sample_n(10) %>% 
                  column_to_rownames('Sample') %>% 
                  as.data.frame() %>%
                  droplevels()
  Y.b[[i]]=exp[rownames(exp) %in% rownames(new_pd[[i]]),]
  Y.b[[i]]=Y.b[[i]][match(rownames(new_pd[[i]]),rownames(Y.b[[i]])),]
}
save(Y.b,new_pd,file = "OUTPUT_DOWNSAMPLE_OLIG2/Downsampled_raw_data.RData")

########################
# Anova across SPECIES #
########################
Model <- vector("list", length = B)
SPECIES_DGE <- vector("list",length = B)
for (i in 1:B)
    {
        SPECIES_DGE[[i]] <- data.frame(matrix(nrow=ncol(Y.b[[i]]), ncol=2, dimnames = list(colnames(Y.b[[i]]), c("Fstat_Species","Anova_Species"))))
        {
        for (j in 1:ncol(Y.b[[i]])) 
            {   
            Model[[i]]=tryCatch(anova(lm(as.formula(paste("Y.b[[i]][,j] ~ ", paste(colnames(new_pd[[i]]),collapse = " + "))), data = new_pd[[i]])),warning =  function(w) w)
                if (j %% 8372 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    SPECIES_DGE[[i]][j,"Fstat_Species"] = Model[[i]]$"F value"[1]
                    SPECIES_DGE[[i]][j,"Anova_Species"] = Model[[i]]$"Pr(>F)"[1]
                    SPECIES_DGE[[i]]$AnovaAdj_Species <- p.adjust(SPECIES_DGE[[i]]$Anova_Species,"bonferroni")
                }
            }
        }
    }
save(SPECIES_DGE,file="OUTPUT_DOWNSAMPLE_OLIG2/SPECIES_DGE.RData")

##########
# H vs C #
##########
new_pd_HC <- vector("list", length = B)
Y.b_HC <- vector("list",length = B)
for (i in 1:B)
  {
  Y.b_HC[[i]] <- Y.b[[i]][grep("Hsap|PanTro",rownames(Y.b[[i]])),]
  new_pd_HC[[i]] <- droplevels(new_pd[[i]][grep("Hsap|PanTro",new_pd[[i]]$Species),])
  }

Model <- vector("list", length = B)
HsapVsPanTro_DGE <- vector("list",length = B)
coefs <- vector("list",length = B)
t_value <- vector("list",length = B)
for (i in 1:B)
    {
        HsapVsPanTro_DGE[[i]] <- data.frame(matrix(nrow=ncol(Y.b_HC[[i]]), ncol=2, dimnames = list(colnames(Y.b_HC[[i]]), c("Estimate_HsapVsPanTro", "Pval_HsapVsPanTro"))))
        {
        for (j in 1:ncol(Y.b_HC[[i]])) 
            {   
            Model[[i]]=tryCatch(lm(as.formula(paste("Y.b_HC[[i]][,j] ~ ", paste(colnames(new_pd_HC[[i]]),collapse = " + "))), data = new_pd_HC[[i]]),warning =  function(w) w)
                if (j %% 8372 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    coefs[[i]] = data.frame(coef(summary(Model[[i]])))
                    t_value[[i]] = coefs[[i]]["Species", "t.value"]
                    HsapVsPanTro_DGE[[i]][j,"Estimate_HsapVsPanTro"] = -1 * coefs[[i]]["Species", "Estimate"]
                    HsapVsPanTro_DGE[[i]][j,"Pval_HsapVsPanTro"] = 2 * (1 - pnorm(abs(t_value[[i]])))
                    HsapVsPanTro_DGE[[i]]$FDR_HsapVsPanTro <- p.adjust(HsapVsPanTro_DGE[[i]]$Pval_HsapVsPanTro,method="BH")
                }
            }
        }
    }
save(HsapVsPanTro_DGE,file="OUTPUT_DOWNSAMPLE_OLIG2/HsapVsPanTro_DGE.RData")

##########
# H vs R #
##########
new_pd_HR <- vector("list", length = B)
Y.b_HR <- vector("list",length = B)
for (i in 1:B)
  {
  Y.b_HR[[i]] <- Y.b[[i]][grep("Hsap|RheMac",rownames(Y.b[[i]])),]
  new_pd_HR[[i]] <- droplevels(new_pd[[i]][grep("Hsap|RheMac",new_pd[[i]]$Species),])
  }

Model <- vector("list", length = B)
HsapVsRheMac_DGE <- vector("list",length = B)
coefs <- vector("list",length = B)
t_value <- vector("list",length = B)
for (i in 1:B)
    {
        HsapVsRheMac_DGE[[i]] <- data.frame(matrix(nrow=ncol(Y.b_HR[[i]]), ncol=2, dimnames = list(colnames(Y.b_HR[[i]]), c("Estimate_HsapVsRheMac", "Pval_HsapVsRheMac"))))
        {
        for (j in 1:ncol(Y.b_HR[[i]])) 
            {   
            Model[[i]]=tryCatch(lm(as.formula(paste("Y.b_HR[[i]][,j] ~ ", paste(colnames(new_pd_HR[[i]]),collapse = " + "))), data = new_pd_HR[[i]]),warning =  function(w) w)
                if (j %% 8372 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    coefs[[i]] = data.frame(coef(summary(Model[[i]])))
                    t_value[[i]] = coefs[[i]]["Species", "t.value"]
                    HsapVsRheMac_DGE[[i]][j,"Estimate_HsapVsRheMac"] = -1 * coefs[[i]]["Species", "Estimate"]
                    HsapVsRheMac_DGE[[i]][j,"Pval_HsapVsRheMac"] = 2 * (1 - pnorm(abs(t_value[[i]])))
                    HsapVsRheMac_DGE[[i]]$FDR_HsapVsRheMac <- p.adjust(HsapVsRheMac_DGE[[i]]$Pval_HsapVsRheMac,method="BH")
                }
            }
        }
    }
save(HsapVsRheMac_DGE,file="OUTPUT_DOWNSAMPLE_OLIG2/HsapVsRheMac_DGE.RData")

##########
# C vs R #
##########
new_pd_CR <- vector("list", length = B)
Y.b_CR <- vector("list",length = B)
for (i in 1:B)
  {
  Y.b_CR[[i]] <- Y.b[[i]][grep("PanTro|RheMac",rownames(Y.b[[i]])),]
  new_pd_CR[[i]] <- droplevels(new_pd[[i]][grep("PanTro|RheMac",new_pd[[i]]$Species),])
  }

Model <- vector("list", length = B)
PanTroVsRheMac_DGE <- vector("list",length = B)
coefs <- vector("list",length = B)
t_value <- vector("list",length = B)
for (i in 1:B)
    {
        PanTroVsRheMac_DGE[[i]] <- data.frame(matrix(nrow=ncol(Y.b_CR[[i]]), ncol=2, dimnames = list(colnames(Y.b_CR[[i]]), c("Estimate_PanTroVsRheMac", "Pval_PanTroVsRheMac"))))
        {
        for (j in 1:ncol(Y.b_CR[[i]])) 
            {   
            Model[[i]]=tryCatch(lm(as.formula(paste("Y.b_CR[[i]][,j] ~ ", paste(colnames(new_pd_CR[[i]]),collapse = " + "))), data = new_pd_CR[[i]]),warning =  function(w) w)
                if (j %% 8372 == 0) {cat(paste("Done on Boot ",i,"\n",sep=""))}
                if(typeof(Model[[i]]) == "list"){
                    coefs[[i]] = data.frame(coef(summary(Model[[i]])))
                    t_value[[i]] = coefs[[i]]["Species", "t.value"]
                    PanTroVsRheMac_DGE[[i]][j,"Estimate_PanTroVsRheMac"] = -1 * coefs[[i]]["Species", "Estimate"]
                    PanTroVsRheMac_DGE[[i]][j,"Pval_PanTroVsRheMac"] = 2 * (1 - pnorm(abs(t_value[[i]])))
                    PanTroVsRheMac_DGE[[i]]$FDR_PanTroVsRheMac <- p.adjust(PanTroVsRheMac_DGE[[i]]$Pval_PanTroVsRheMac,method="BH")
                }
            }
        }
    }
save(PanTroVsRheMac_DGE,file="OUTPUT_DOWNSAMPLE_OLIG2/PanTroVsRheMac_DGE.RData")

# Combine data
for (i in 1:B)
    {
      HsapVsPanTro_DGE[[i]] <- HsapVsPanTro_DGE[[i]] %>% rownames_to_column('Gene')
      HsapVsRheMac_DGE[[i]] <- HsapVsRheMac_DGE[[i]] %>% rownames_to_column('Gene')
      PanTroVsRheMac_DGE[[i]] <- PanTroVsRheMac_DGE[[i]] %>% rownames_to_column('Gene')
      SPECIES_DGE[[i]] <- SPECIES_DGE[[i]] %>% rownames_to_column('Gene')
    }

myfiles <- vector("list", length = B)
for (i in 1:B)
    {
    myfiles[[i]] <- Reduce(function(x, y) merge(x, y,by="Gene",all=TRUE), list(HsapVsPanTro_DGE[[i]], HsapVsRheMac_DGE[[i]], PanTroVsRheMac_DGE[[i]],SPECIES_DGE[[i]]))
    }
save(myfiles,file="OUTPUT_DOWNSAMPLE_OLIG2/Combined_DGE.RData")


# Human Specific 
posH <- vector("list", length = B)
negH <- vector("list", length = B)

for (i in 1:B)
    {
      posH[[i]] <- myfiles[[i]] %>% 
      filter(FDR_PanTroVsRheMac > 0.1 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro > 0 & 
        Estimate_HsapVsRheMac  > 0) %>%
        as.data.frame()

      negH[[i]] <- myfiles[[i]] %>%
      filter(FDR_PanTroVsRheMac > 0.1 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro < 0 & 
        Estimate_HsapVsRheMac < 0) %>%
        as.data.frame()
}


human_specific <- Map(rbind, posH, negH)
save(human_specific,file="OUTPUT_DOWNSAMPLE_OLIG2/Hsap_Specific_DGE.RData")


# Chimp Specific 
posC <- vector("list", length = B)
negC <- vector("list", length = B)

for (i in 1:B)
    {
      posC[[i]] <- myfiles[[i]] %>% 
      filter(FDR_PanTroVsRheMac > 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.1 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro < 0 & 
        Estimate_PanTroVsRheMac  > 0) %>%
        as.data.frame()

      negC[[i]] <- myfiles[[i]] %>%
      filter(FDR_PanTroVsRheMac > 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.1 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro > 0 & 
        Estimate_PanTroVsRheMac  < 0) %>%
        as.data.frame()
}


pantro_specific <- Map(rbind, posC, negC)
save(pantro_specific,file="OUTPUT_DOWNSAMPLE_OLIG2/PanTro_Specific_DGE.RData")


# RheMac Specific 
posR <- vector("list", length = B)
negR <- vector("list", length = B)

for (i in 1:B)
    {
      posR[[i]] <- myfiles[[i]] %>% 
      filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro > 0.1 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_PanTroVsRheMac < 0 & 
        Estimate_HsapVsRheMac < 0) %>%
        as.data.frame()

      negR[[i]] <- myfiles[[i]] %>%
      filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro > 0.1 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_PanTroVsRheMac > 0 & 
        Estimate_HsapVsRheMac  > 0) %>%
        as.data.frame()
}


rhemac_specific <- Map(rbind, posR, negR)
save(rhemac_specific,file="OUTPUT_DOWNSAMPLE_OLIG2/RheMac_Specific_DGE.RData")
































