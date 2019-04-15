
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

# Boxplot difference estimate size
n <- read.table(here("OUTPUTS_EXONS_NEUN","NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt"),header=T,sep="\t")
o <- read.table(here("OUTPUTS_EXONS_OLIG2","OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt"),header=T,sep="\t")

n$Class <- rep("NeuN",nrow(n))
o$Class <- rep("OLIG2",nrow(o))

df <- rbind(n,o)

df$absHC <- abs(df$Estimate_HsapVsPanTro)
df$absHR <- abs(df$Estimate_HsapVsRheMac)

A <- ggviolin(df, 
      		x = "Class", 
      		y = "Estimate_HsapVsPanTro",
          	fill = "Class", 
          	add = "boxplot", add.params = list(fill = "white"),
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="Human vs Chimpanzee", subtitle= "Human Specific DGE",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)

C <- ggviolin(df, 
      		x = "Class", 
      		y = "Estimate_HsapVsRheMac",
          	fill = "Class", 
          	add = "boxplot", add.params = list(fill = "white"), 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="Human vs Macaque",subtitle= "Human Specific DGE",x="", y = "log2(Fold Change)")+
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

B <- ggviolin(df2, 
      		x = "Class", 
      		y = "Estimate_HsapVsPanTro",
          	fill = "Class", 
          	add = "boxplot", add.params = list(fill = "white"), 
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="Human vs Chimpanzee", subtitle= "All Genes",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)

D <- ggviolin(df2, 
      		x = "Class", 
      		y = "Estimate_HsapVsRheMac",
          	fill = "Class", 
          	add = "boxplot", add.params = list(fill = "white"),
          	palette = c("cyan4","magenta4"),
          	short.panel.labs = FALSE) + 
          	stat_compare_means(label = "p.format",method = "t.test",label.x = 1.5)+
          	theme_classic()+
        		labs(title="Human vs Macaque",subtitle= "All Genes",x="", y = "log2(Fold Change)")+
        		theme(legend.position="none")+
          	rotate_x_text(angle = 45)


plot2by2 <- cowplot::plot_grid(A,B,C,D,labels=c("A", "B","C","D"), ncol = 2)
cowplot::save_plot("Boxplots_Estimates.pdf", plot2by2, ncol = 2,base_height=6,base_width=3)

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

