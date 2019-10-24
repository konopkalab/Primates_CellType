library(gplots)
library(graphics)
library(ape)
library(vcd)

rm(list=ls())
# unrooted tree
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
