suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
library(ape)


load("OUTPUT_DOWNSAMPLE_NeuN/Hsap_Specific_DGE.RData")
load("OUTPUT_DOWNSAMPLE_NeuN/PanTro_Specific_DGE.RData")
load("OUTPUT_DOWNSAMPLE_NeuN/RheMac_Specific_DGE.RData")


round(mean(sapply(human_specific, nrow)))/6 #25.5
round(mean(sapply(pantro_specific, nrow)))/6 #11.3
round(mean(sapply(rhemac_specific, nrow)))/45 #10.3

# TREE
n <- read.tree(text="(RheMac:10.3,(PanTro:11.3,Hsap:25.5):1);")

pdf("ChangedXyear_Tree_NeuN_Downsample_Mean.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(n), type = "unrooted", cex = 0.6,
     edge.color = "cyan4", edge.width = 2, edge.lty = 2,
     tip.color = "cyan4")
edgelabels(n$edge.length, bg = "white",col="black", font=2)
dev.off()

rm(list=ls())

load("OUTPUT_DOWNSAMPLE_OLIG2/Hsap_Specific_DGE.RData")
load("OUTPUT_DOWNSAMPLE_OLIG2/PanTro_Specific_DGE.RData")
load("OUTPUT_DOWNSAMPLE_OLIG2/RheMac_Specific_DGE.RData")


round(mean(sapply(human_specific, nrow)))/6 #37.8
round(mean(sapply(pantro_specific, nrow)))/6 #17.5
round(mean(sapply(rhemac_specific, nrow)))/45 #10.7

o <- read.tree(text="(RheMac:10.7,(PanTro:17.5,Hsap:37.8):1);")

pdf("ChangedXyear_Tree_OLIG2_Downsample_Mean.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(o), type = "unrooted", cex = 0.6,
     edge.color = "magenta4", edge.width = 2, edge.lty = 2,
     tip.color = "magenta4")
edgelabels(o$edge.length, bg = "white",col="black", font=2)
dev.off()


# Isabel Tree
n <- read.tree(text="(RheMac:3.9,(PanTro:3.5,Hsap:9.5):1);")

pdf("ChangedXyear_Tree_NeuN_Downsample_Isabel.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(n), type = "unrooted", cex = 0.6,
     edge.color = "cyan4", edge.width = 2, edge.lty = 2,
     tip.color = "cyan4")
edgelabels(n$edge.length, bg = "white",col="black", font=2)
dev.off()


o <- read.tree(text="(RheMac:3.2,(PanTro:3,Hsap:12.7):1);")

pdf("ChangedXyear_Tree_OLIG2_Downsample_Isabel.pdf",width=5,height=5,useDingbats=FALSE)
plot(as.phylo(o), type = "unrooted", cex = 0.6,
     edge.color = "magenta4", edge.width = 2, edge.lty = 2,
     tip.color = "magenta4")
edgelabels(o$edge.length, bg = "white",col="black", font=2)
dev.off()
