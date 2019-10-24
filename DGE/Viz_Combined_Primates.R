
rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(graphics))

# Boxplot difference estimate size
n <- read.table(here("OUTPUTS_EXONS_NEUN","NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt"),header=T,sep="\t")
o <- read.table(here("OUTPUTS_EXONS_OLIG2","OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt"),header=T,sep="\t")

n$Class <- rep("NeuN",nrow(n))
o$Class <- rep("OLIG2",nrow(o))

df <- rbind(n,o)

df$absHC <- abs(df$Estimate_HsapVsPanTro)
df$absHR <- abs(df$Estimate_HsapVsRheMac)

A <- ggboxplot(df, 
      		x = "Class", 
      		y = "Estimate_HsapVsPanTro",
          	color = "Class", 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="H vs C: DGE",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)

C <- ggboxplot(df, 
      		x = "Class", 
      		y = "Estimate_HsapVsRheMac",
          	color = "Class", 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="H vs M: DGE",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)


# Boxplot difference estimate size
n <- read.table(here("OUTPUTS_EXONS_NEUN","NEUN_PRIMATES_DGE.txt"),header=T,sep="\t")
o <- read.table(here("OUTPUTS_EXONS_OLIG2","OLIG2_PRIMATES_DGE.txt"),header=T,sep="\t")

n$Class <- rep("NeuN",nrow(n))
o$Class <- rep("OLIG2",nrow(o))

df2 <- rbind(n,o)

df2$absHC <- abs(df2$Estimate_HsapVsPanTro)
df2$absHR <- abs(df2$Estimate_HsapVsRheMac)

B <- ggboxplot(df2, 
      		x = "Class", 
      		y = "Estimate_HsapVsPanTro",
          	color = "Class", 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="H vs C: All Genes",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)

D <- ggboxplot(df2, 
      		x = "Class", 
      		y = "Estimate_HsapVsRheMac",
          	color = "Class", 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="H vs M: All Genes",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)


plot2by2 <- cowplot::plot_grid(A,B,C,D,labels=c("A", "B","C","D"), ncol = 2)
cowplot::save_plot("Boxplots_Estimates.pdf", plot2by2, ncol = 2,base_height=5,base_width=2.5)

# Barplot number of DGE
data <- data.frame(Species = c("Hsap","PanTro","RheMac","Hsap","PanTro","RheMac"),CellType = c(rep("NeuN",3),rep("OLIG2",3)), Up=c(313,141,437,354,99,337),Down=c(181,134,457,179,64,293))

df <- melt(data)
df$Species = factor(df$Species, levels=c("RheMac","PanTro","Hsap"))
df$variable = factor(df$variable, levels=c("Down","Up"))

pdf("BarPlot_DGE_Both_Oriz.pdf",width=5,height=5,useDingbats=FALSE)
ggbarplot(df, "Species", "value",
  fill = "variable", color = "variable", palette = c("steelblue","red"),
  label = TRUE,
  position = position_dodge(0.8),
  orientation = "horiz",
  lab.pos = "out", lab.col = "black",
  xlab = "", ylab = "", facet.by = "CellType")+
theme_classic()+
ylim(0,700)
dev.off()

# unrooted tree
library(ape)
n <- read.tree(text="(RheMac:19.9,(PanTro:45.8,Hsap:82.4):1);")

pdf("ChangesXyear_Tree_NeuN_New.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(n), type = "unrooted", cex = 0.6,
     edge.color = "cyan4", edge.width = 2, edge.lty = 2,
     tip.color = "cyan4")
edgelabels(n$edge.length, bg = "white",col="black", font=2)
dev.off()

o <- read.tree(text="(RheMac:14,(PanTro:27.2,Hsap:88.8):1);")

pdf("ChangesXyear_Tree_OLIG2_New.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(o), type = "unrooted", cex = 0.6,
     edge.color = "magenta4", edge.width = 2, edge.lty = 2,
     tip.color = "magenta4")
edgelabels(o$edge.length, bg = "white",col="black", font=2)
dev.off()

# Matrix 
mat <- as.data.frame(matrix(c(82.4,45.8, 19.9,88.8,27.2,14),nc=2))
rownames(mat) <- c("H","P","M")
colnames(mat) <- c("NeuN","OLIG2")
chisq.test(mat)

# Density of gene expression
load("DGE_DATA_INPUTS/DGE_Primates_NeuN.RData")
load("OUTPUTS_EXONS_NEUN/WGCNA_EXONS_NEUN.RData")

dge <- AllStat
adj <- expAdj

adj$Gene <- rownames(adj)
NEUN <- melt(adj)
NEUN$Species <- do.call(rbind,strsplit(as.character(NEUN$variable),"_"))[,2]
NEUN$Class <- rep("NeuN",nrow(NEUN))

load("DGE_DATA_INPUTS/DGE_Primates_OLIG2.RData")
load("OUTPUTS_EXONS_OLIG2/WGCNA_EXONS_OLIG2.RData")

dge <- AllStat
adj <- expAdj

adj$Gene <- rownames(adj)
OLIG2 <- melt(adj)
OLIG2$Species <- do.call(rbind,strsplit(as.character(OLIG2$variable),"_"))[,2]
OLIG2$Class <- rep("OLIG2",nrow(OLIG2))

df <- rbind(NEUN,OLIG2)

pdf("Density_Expression_Values.pdf",width=4,height=2,useDingbats=FALSE)
ggdensity(df, x = "value",
   add = "mean", rug = FALSE,
   color = "Class", palette = c("cyan4", "magenta4"),
   facet.by="Species",
   alpha=0.5) +
theme_classic()+
xlab("Adjusted Expression")+
theme(legend.position="none")
dev.off()

# Mosaic Plot and Fisher exact
M <- as.table(rbind(c(494,533), c(7878,7027)))
dimnames(M) <- list(Class = c("DEG", "NotDEG"),
                    Cell = c("NeuN","OLIG2"))

fisher.test(M,alternative="l")

pdf("MosaicPlot_Hsap.pdf",width=4,height=4,useDingbats=FALSE)
mosaicplot(M, shade = TRUE, las=2,main = "Hsap")
dev.off()


M <- as.table(rbind(c(275,163), c(8097,7397)))
dimnames(M) <- list(Class = c("DEG", "NotDEG"),
                    Cell = c("NeuN","OLIG2"))

fisher.test(M,alternative="l")

pdf("MosaicPlot_PanTro.pdf",width=4,height=4,useDingbats=FALSE)
mosaicplot(M, shade = TRUE, las=2,main = "PanTro")
dev.off()

M <- as.table(rbind(c(894,630), c(7478,6930)))
dimnames(M) <- list(Class = c("DEG", "NotDEG"),
                    Cell = c("NeuN","OLIG2"))

pdf("MosaicPlot_MacMul.pdf",width=4,height=4,useDingbats=FALSE)
mosaicplot(M, shade = TRUE, las=2,main = "PanTro")
dev.off()

