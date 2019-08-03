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
source("Utility_Functions.R")

####################################
# Load Inputs and filter for Neun+ #
####################################
dir.create("OUTPUTS_EXONS_NEUN/")
load(list.files(pattern = "NeuN",path = "PRIMATE_EXPRESSION_INPUT",full.names = TRUE))

data <- HCR_NeuN_Exons
pd <- HCR_NeuN_pd

# Convert Age into Factor
pd$HumAge <- ifelse(pd$HumAge <= 40,"1", ifelse (pd$HumAge > 40 & pd$HumAge <= 60 ,"2", ifelse(pd$HumAge > 60,"3","3")))
pd$HumAge <- as.factor(pd$HumAge)


##############################
# Low expressed gene removal #
##############################
#class <- levels(pd$Species)
#filter <- apply(data, 1, function(x) (all(x[grep(paste(class[1]),names(x))] > 0) | all(x[grep(paste(class[2]),names(x))] > 0) | all(x[grep(paste(class[3]),names(x))] > 0)))
#data <- data[filter,]

##########################
#    normalization       #
##########################
dat <- log2(cpm(data)+1)

##########################
# Quantile normalization #
##########################
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
write.table(p, "OUTPUTS_EXONS_NEUN/NeuN_Primates_Exons_CPM.txt",sep="\t",quote=F)

#########################
# Variance explained    #
#########################
var <- VarExp(p,pd,5,FALSE)
pdf("OUTPUTS_EXONS_NEUN/Variance_Explained_NeuN.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

#########################
#         PCA           #
#########################
pdf("OUTPUTS_EXONS_NEUN/NeuN_Primates_PCA.pdf",width=5,height=5,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Species=pd$Species)
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",color = "Species",palette=c("steelblue","darkgrey","green"), shape = 21, size = 3)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()

#########################
# Surrogates Variables  #
#########################
mod <- model.matrix(~Species+Sex+HumAge+RIN, data =pd)
mod0 <- model.matrix(~Sex+HumAge+RIN, data = pd)
svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100,method="two-step")
pdSv <- cbind(pd,svaobj$sv)

#########################
# Expression adjustment #
#########################
svs <- svaobj$sv
colnames(svs) <- paste("SV",1:ncol(svs),sep="")
pd_sva <- cbind(pd[c(-1)],svs)
pd_sva <- droplevels(pd_sva)
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~.,data=pd_sva)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals)<-rownames(p)
rownames(residuals) <- rownames(p)
write.table(adj.residuals, "OUTPUTS_EXONS_NEUN/NeuN_Primates_AdjExp.txt",sep="\t",quote=F)
pAdj <- adj.residuals

########################
# Anova across SPECIES #
########################
TRAITSr=cbind(pd,svs)
TRAITSr <- droplevels(TRAITSr)

# ANOVA
Data=t(p)
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("Fstat","Anova", "Warning"))))
output[,] <- NA

#' Looping on the Data
for (i in 1:ncol(Data)) {   
            Model=tryCatch(anova(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSr),collapse = " + "))), data = TRAITSr)),warning =  function(w) w)
                if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
                if(typeof(Model) == "list"){
                    output[i,"Fstat"] = Model$"F value"[1]
                    output[i,"Anova"] = Model$"Pr(>F)"[1]
                    } else {
                    output[i,"Warning"] = as.character(Model)
                    output[i, "Fstat"] = 0
                    output[i,"Anova"] = 1
            }
        }

SPECIES_DGE <- output
SPECIES_DGE$Warning <- NULL
SPECIES_DGE$FDR <- p.adjust(SPECIES_DGE$Anova,"bonferroni")
colnames(SPECIES_DGE) <- c("Fstat_Species","Anova_Species","AnovaAdj_Species")
write.table(SPECIES_DGE, "OUTPUTS_EXONS_NEUN/SPECIES_DGE.txt",sep="\t",quote=F)

########################
# Human vs Chimpanzee  #
########################
TRAITSfilt <- TRAITSr[grep("Hsap|PanTro",TRAITSr$Species),]
TRAITSfilt <- droplevels(TRAITSfilt)
Data=t(p[,grep("Hsap|PanTro",colnames(p))])
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("Estimate", "Pval", "Warning"))))
output[,] <- NA
for (i in 1:ncol(Data)) {   
            Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSfilt),collapse = " + "))), data = TRAITSfilt),warning =  function(w) w)
                if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
                if(typeof(Model) == "list"){
                    coefs = data.frame(coef(summary(Model)))
                    t_value = coefs["Species", "t.value"]
                    output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
                    output[i,"Estimate"]= -1 * coefs["Species", "Estimate"]
                    } else {
                    output[i,"Warning"] = as.character(Model)
                    output[i, "Estimate"] = 0
                    output[i,"Pval"] = 1
            }
        }

HsapVsPanTro_DGE <- output
HsapVsPanTro_DGE$Warning <- NULL
HsapVsPanTro_DGE$FDR <- p.adjust(HsapVsPanTro_DGE$Pval,"BH")
colnames(HsapVsPanTro_DGE) <- c("Estimate_HsapVsPanTro","Pval_HsapVsPanTro","FDR_HsapVsPanTro")
write.table(HsapVsPanTro_DGE, "OUTPUTS_EXONS_NEUN/HsapVsPanTro_DGE.txt",sep="\t",quote=F)

########################
# Human vs Macaque     #
########################
TRAITSfilt <- TRAITSr[grep("Hsap|RheMac",TRAITSr$Species),]
TRAITSfilt <- droplevels(TRAITSfilt)
Data=t(p[,grep("Hsap|RheMac",colnames(p))])
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("Estimate", "Pval", "Warning"))))
output[,] <- NA

#' Looping on the Data
for (i in 1:ncol(Data)) {   
            Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSfilt),collapse = " + "))), data = TRAITSfilt),warning =  function(w) w)
                if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
                if(typeof(Model) == "list"){
                    coefs = data.frame(coef(summary(Model)))
                    t_value = coefs["Species", "t.value"]
                    output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
                    output[i,"Estimate"]= -1 * coefs["Species", "Estimate"]
                    } else {
                    output[i,"Warning"] = as.character(Model)
                    output[i, "Estimate"] = 0
                    output[i,"Pval"] = 1
            }
        }

HsapVsRheMac_DGE <- output
HsapVsRheMac_DGE$Warning <- NULL
HsapVsRheMac_DGE$FDR <- p.adjust(HsapVsRheMac_DGE$Pval,"BH")
colnames(HsapVsRheMac_DGE) <- c("Estimate_HsapVsRheMac","Pval_HsapVsRheMac","FDR_HsapVsRheMac")
write.table(HsapVsRheMac_DGE, "OUTPUTS_EXONS_NEUN/HsapVsRheMac_DGE.txt",sep="\t",quote=F)

#########################
# Chimpanzee vs Macaque #
#########################
TRAITSfilt <- TRAITSr[grep("PanTro|RheMac",TRAITSr$Species),]
TRAITSfilt <- droplevels(TRAITSfilt)
Data=t(p[,grep("PanTro|RheMac",colnames(p))])
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("Estimate", "Pval", "Warning"))))
output[,] <- NA

#' Looping on the Data
for (i in 1:ncol(Data)) {   
            Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSfilt),collapse = " + "))), data = TRAITSfilt),warning =  function(w) w)
                if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
                if(typeof(Model) == "list"){
                    coefs = data.frame(coef(summary(Model)))
                    t_value = coefs["Species", "t.value"]
                    output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
                    output[i,"Estimate"]= -1 * coefs["Species", "Estimate"]
                    } else {
                    output[i,"Warning"] = as.character(Model)
                    output[i, "Estimate"] = 0
                    output[i,"Pval"] = 1
            }
        }

PanTroVsRheMac_DGE <- output
PanTroVsRheMac_DGE$Warning <- NULL
PanTroVsRheMac_DGE$FDR <- p.adjust(PanTroVsRheMac_DGE$Pval,"BH")
colnames(PanTroVsRheMac_DGE) <- c("Estimate_PanTroVsRheMac","Pval_PanTroVsRheMac","FDR_PanTroVsRheMac")
write.table(PanTroVsRheMac_DGE, "OUTPUTS_EXONS_NEUN/PanTroVsRheMac_DGE.txt",sep="\t",quote=F)

################
# DGE database #
################
files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = 'DGE.txt',full.names = TRUE)
myfiles = lapply(files, read.table,header=T,sep="\t")

for(i in 1:length(myfiles)){
    myfiles[[i]]$Gene <- rownames(myfiles[[i]])
    rownames(myfiles[[i]]) <- NULL
}

res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, myfiles)
write.table(res, "OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE.txt",sep="\t",quote=F)
write.xlsx(res, file="OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE.xlsx",sheetName = "GeneSetsDb_NEUN",row.names=FALSE)

#########################
# Species-specific DGE  #
#########################
files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_PRIMATES_DGE.txt',full.names = TRUE)
dge <- read.table(files,header=T,sep="\t")

#Human
posH <- dge %>% 
    filter(FDR_PanTroVsRheMac > 0.1 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro > 0.3 & 
        Estimate_HsapVsRheMac  > 0.3) %>%
    as.data.frame()

negH <- dge %>%
    filter(FDR_PanTroVsRheMac > 0.1 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro < -0.3 & 
        Estimate_HsapVsRheMac  < -0.3) %>%
    as.data.frame()
dfH <- rbind(posH,negH)

files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_CPM.txt',full.names = TRUE)
exp = lapply(files, read.table,header=T,sep="\t")[[1]]

mat <- exp[rownames(exp)%in%dfH$Gene,]
pdf("OUTPUTS_EXONS_NEUN/Human_Genes_Heatmap.pdf",width=4,height=4)
pheatmap(mat,cluster_cols=T,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),annotation=pd, show_rownames = F,show_colnames=F)
dev.off()

write.table(dfH, "OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt",sep="\t",quote=F)
write.xlsx(dfH, file="OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_SIGN.xlsx",sheetName = "HUMAN SPECIFIC",row.names=FALSE,append=TRUE, showNA=FALSE)

# Primate Accelerated
posACC <- dge %>% 
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro > 0.3 & 
        Estimate_HsapVsRheMac  > 0.3 &
        Estimate_PanTroVsRheMac  > 0.3 ) %>%
    as.data.frame()

negACC <- dge %>%
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro < -0.3 & 
        Estimate_HsapVsRheMac  < -0.3 &
        Estimate_PanTroVsRheMac < -0.3) %>%
    as.data.frame()
dfACC <- rbind(posACC,negACC)
write.table(dfACC, "OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_PRIMATE_ACCELERATED.txt",sep="\t",quote=F)
write.xlsx(dfACC, file="OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_SIGN.xlsx",sheetName = "PRIMATE ACCELERATED",row.names=FALSE,append=TRUE, showNA=FALSE)

# Chimp
posC <- dge %>% 
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac > 0.1 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro < -0.3 & 
        Estimate_PanTroVsRheMac > 0.3) %>%
    as.data.frame()
negC <- dge %>%
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro < 0.05 & 
        FDR_HsapVsRheMac > 0.1 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_HsapVsPanTro > 0.3 & 
        Estimate_PanTroVsRheMac  < -0.3) %>%
    as.data.frame()
dfC <- rbind(posC,negC)

files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_CPM.txt',full.names = TRUE)
exp = lapply(files, read.table,header=T,sep="\t")[[1]]

mat <- exp[rownames(exp)%in%dfC$Gene,]
pdf("OUTPUTS_EXONS_NEUN/Chimp_Genes_Heatmap.pdf",width=4,height=4)
pheatmap(mat,cluster_cols=T,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),annotation=pd, show_rownames = F,show_colnames=F)
dev.off()


write.table(dfC, "OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_CHIMP_SPECIFIC.txt",sep="\t",quote=F)
write.xlsx(dfC, file="OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_SIGN.xlsx",sheetName = "CHIMP SPECIFIC",row.names=FALSE,append=TRUE, showNA=FALSE)

# Macaque
posR <- dge %>% 
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro > 0.1 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_PanTroVsRheMac < -0.3 & 
        Estimate_HsapVsRheMac < -0.3) %>%
    as.data.frame()
negR <- dge %>%
    filter(FDR_PanTroVsRheMac < 0.05 & 
        FDR_HsapVsPanTro > 0.1 & 
        FDR_HsapVsRheMac < 0.05 & 
        AnovaAdj_Species < 0.05 & 
        Estimate_PanTroVsRheMac > 0.3 & 
        Estimate_HsapVsRheMac  > 0.3) %>%
    as.data.frame()
dfR <- rbind(posR,negR)

files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_CPM.txt',full.names = TRUE)
exp = lapply(files, read.table,header=T,sep="\t")[[1]]

mat <- exp[rownames(exp)%in%dfR$Gene,]
pdf("OUTPUTS_EXONS_NEUN/Rhesus_Genes_Heatmap.pdf",width=4,height=4)
pheatmap(mat,cluster_cols=T,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),annotation=pd, show_rownames = F,show_colnames=F)
dev.off()

write.table(dfR, "OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_MACAQUE_SPECIFIC.txt",sep="\t",quote=F)
write.xlsx(dfR, file="OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_SIGN.xlsx",sheetName = "MACAQUE SPECIFIC",row.names=FALSE,append=TRUE, showNA=FALSE)

####################
# WGCNA input data #
####################
files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_AdjExp.txt',full.names = TRUE)
expAdj = lapply(files, read.table,header=T,sep="\t")[[1]]
files = list.files(path = "OUTPUTS_EXONS_NEUN/",pattern = '_CPM.txt',full.names = TRUE)
expQuant = lapply(files, read.table,header=T,sep="\t")[[1]]

tmp <- pd
tmp$Hsap <- ifelse(tmp$Species == "Hsap",1,0)
tmp$PanTro <- ifelse(tmp$Species == "PanTro",1,0)
tmp$RheMac <- ifelse(tmp$Species == "RheMac",1,0)
tmp$Sex <- ifelse(tmp$Sex == "M",1,0)
tmp$Species <- NULL

pheno <- tmp
save(expAdj,expQuant,pheno,pd,file = "OUTPUTS_EXONS_NEUN/WGCNA_EXONS_NEUN.RData")

