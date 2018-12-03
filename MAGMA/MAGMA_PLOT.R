library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)

files=list.files(pattern="STATISTICS")
tmp=as.data.frame(lapply(files,read.table,sep="\t",header=T)[[1]])

tmp <- tmp[c(1,7,8)]
df <- melt(tmp)
df$log <- -log10(df$value)

pdf("MAGMA_ENRICH_BarPlot.pdf",width=8,height=5,useDingbats=FALSE)
ggbarplot(df, "Sample", "log",
  fill = "VARIABLE", color = "VARIABLE", palette = c("cyan4","magenta4"),
  label = FALSE,
  position = position_dodge(0.9))+
geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1) +
ylim(0,5)+ 
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
ylab("-log10(FDR)")+
xlab("")+
theme(legend.position = c(0.15,0.8)) # remove legend
dev.off()







