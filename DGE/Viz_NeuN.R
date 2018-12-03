# Viz Hsap
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))


dir.create("NeuN_Viz")
load("DGE_DATA_INPUTS/DGE_Primates_NeuN.RData")
load("OUTPUTS_EXONS_NEUN/WGCNA_EXONS_NEUN.RData")

hsp <- Hspecific
hsp$Direction <- ifelse(hsp$Estimate_HsapVsPanTro > 0,"UP","FALSE")
top_labelled <- tbl_df(hsp) %>% group_by(Direction) %>% top_n(n = 10, wt = abs(Estimate_HsapVsPanTro))

pdf("NeuN_Viz/SCATTER_Estimates_Hsap_NeuN.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(hsp, x = "Estimate_HsapVsPanTro", y = "Estimate_HsapVsRheMac",
   color = "Direction",palette=c("steelblue","red"),size = 0.5,
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
xlab("log2FC Hsap - PanTro")+ 
ylab("log2FC Hsap - RheMac")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
geom_text_repel(data = top_labelled, mapping = aes(label = Gene), size = 4,box.padding = unit(0.4, "lines"),point.padding = unit(0.4, "lines"))+
ylim(-6,6) +
xlim(-6,6) +
theme(legend.position="none")
dev.off()

pdf("NeuN_Viz/SCATTER_Estimates_Hsap_NeuN_NoLabels.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(hsp, x = "Estimate_HsapVsPanTro", y = "Estimate_HsapVsRheMac",
   color = "Direction",palette=c("steelblue","red"),size = 0.5,
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
theme(legend.position="none")+
xlab("log2FC Hsap - PanTro")+ 
ylab("log2FC Hsap - RheMac")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
ylim(-6,6) +
xlim(-6,6) 
dev.off()

# Rendered PCA
PC<-prcomp(t(expQuant),scale=F)
eig <- (PC$sdev)^2
variance <- eig*100/sum(eig)

# Get the name
PCi<-data.frame(PC$x[,1:3],Class = pd$Species)
pdf("NeuN_Viz/NeuN_Primates_PCA.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(PCi, x = "PC1", y = "PC2",color = "Class",
palette = c("steelblue", "grey60","green"),
ellipse = TRUE, mean.point = FALSE,star.plot = FALSE,rug = FALSE)+
ggtitle("Principal Component Analysis: NeuN")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )")) +
theme(legend.position="none")+
border(color = "#008b8b")
dev.off()


# Avarege expression line plot
dge <- AllStat
adj <- expQuant

pos <- dge %>% 
filter(FDR_PanTroVsRheMac > 0.1 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro > 0.3 & Estimate_HsapVsRheMac > 0.3) %>%
as.data.frame()

neg <- dge %>%
filter(FDR_PanTroVsRheMac > 0.1 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro < -0.3 & Estimate_HsapVsRheMac  < -0.3) %>%
as.data.frame()

df <- rbind(pos,neg)

tmp1 <- adj[rownames(adj)%in%pos$Gene,]
tmp2 <- adj[rownames(adj)%in%neg$Gene,]

tmp1$Gene <- rownames(tmp1)
tmp1$Direction<- rep("Positive",nrow(tmp1))

tmp2$Gene <- rownames(tmp2)
tmp2$Direction<- rep("Negative",nrow(tmp2))

tmp <- as.data.frame(rbind(tmp1,tmp2))
tmp.m <- melt(tmp)
tmp.m$Species <- do.call(rbind,strsplit(as.character(tmp.m$variable),"_"))[,2]


pdf("NeuN_Viz/NEUN_PRIMATES_DGE_MEANEXP.pdf",width=4,height=3,useDingbats=FALSE)
ggline(tmp.m, x = "Species", y = "value", color = "Direction",
add = c("mean_se"), palette = c("steelblue", "red"))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
theme(legend.position="right")+
border(color = "#008b8b")+
ggtitle("DGE Species: NeuN")+
xlab("")+ 
ylab("Expression: Mean")
dev.off()

# RADAR plot
dge <- AllStat

#Human
posH <- dge %>% 
    filter(FDR_PanTroVsRheMac > 0.1 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro > 0.3 & Estimate_HsapVsRheMac  > 0.3) %>%
    as.data.frame()

negH <- dge %>%
    filter(FDR_PanTroVsRheMac > 0.1 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro < -0.3 & Estimate_HsapVsRheMac  < -0.3) %>%
    as.data.frame()

# Chimp
posC <- dge %>% 
    filter(FDR_PanTroVsRheMac < 0.05 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac > 0.1 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro < -0.3 & Estimate_PanTroVsRheMac > 0.3) %>%
    as.data.frame()
negC <- dge %>%
    filter(FDR_PanTroVsRheMac < 0.05 & FDR_HsapVsPanTro < 0.05 & FDR_HsapVsRheMac > 0.1 & AnovaAdj_Species < 0.05 & Estimate_HsapVsPanTro > 0.3 & Estimate_PanTroVsRheMac  < -0.3) %>%
    as.data.frame()

# Macaque
posR <- dge %>% 
    filter(FDR_PanTroVsRheMac < 0.05 & FDR_HsapVsPanTro > 0.1 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_PanTroVsRheMac < -0.3 & Estimate_HsapVsRheMac < -0.3) %>%
    as.data.frame()
negR <- dge %>%
    filter(FDR_PanTroVsRheMac < 0.05 & FDR_HsapVsPanTro > 0.1 & FDR_HsapVsRheMac < 0.05 & AnovaAdj_Species < 0.05 & Estimate_PanTroVsRheMac > 0.3 & Estimate_HsapVsRheMac  > 0.3) %>%
    as.data.frame()

df <- data.frame(row.names = c("Human","Chimpanzee","Macaque"), Up = c(nrow(posH),nrow(posC),nrow(posR)),Down = c(nrow(negH),nrow(negC),nrow(negR)))

df$Dge <- df$Up + df$Down
df$Dummy <- rep(0,3)
tab <- (df*100)/nrow(dge)

plot.data <- 
  tab %>%
  as.data.frame %>%
  mutate(Species=rownames(tab))  %>%
  melt(id.vars='Species')

# create new coord : inherit coord_polar
coord_radar <- 
  function(theta='x', start=0, direction=1){
    # input parameter sanity check
    match.arg(theta, c('x','y'))

    ggproto(
      NULL, CoordPolar, 
      theta=theta, r=ifelse(theta=='x','y','x'),
      start=start, direction=sign(direction),
      is_linear=function() TRUE)
  }

pdf("NeuN_Viz/RADAR_NeuN.pdf",width=5,height=5,useDingbats=FALSE)
ggplot(plot.data,aes(x=Species, y=value, group=variable, colour=variable)) + 
geom_polygon(fill=NA) +
geom_point(size=2) +
geom_text_repel(data = subset(plot.data, value > 0),aes(label=paste(round(value,1),"%")),fontface = 'bold',box.padding = 0.25, point.padding = 0.25,segment.color = 'grey50')+
coord_radar() + 
theme_minimal()+
theme(axis.ticks =element_blank(), axis.text.y =element_blank(),axis.title=element_blank(), axis.text.x=element_text())+
scale_color_manual(values=c("red", "steelblue","darkgrey","white"))+
ggtitle("Species Specific DGE: NeuN")+
border(color = "#008b8b")+
theme(legend.position = c(0.9, 0.1))
dev.off()

# with no DGE
plot.data2 <- plot.data[plot.data$variable != "Dge",]
plot.data2 <- droplevels(plot.data2)
pdf("NeuN_Viz/RADAR_NeuN_NoDGE.pdf",width=5,height=5,useDingbats=FALSE)
ggplot(plot.data2,aes(x=Species, y=value, group=variable, colour=variable)) + 
geom_polygon(fill=NA) +
geom_point(size=2) +
geom_text_repel(data = subset(plot.data2, value > 0),aes(label=paste(round(value,1),"%")),fontface = 'bold',box.padding = 0.25, point.padding = 0.25,segment.color = 'grey50')+
coord_radar() + 
theme_minimal()+
theme(axis.ticks =element_blank(), axis.text.y =element_blank(),axis.title=element_blank(), axis.text.x=element_text())+
scale_color_manual(values=c("red", "steelblue","white"))+
ggtitle("Species Specific DGE: NeuN")+
border(color = "#008b8b")+
theme(legend.position = c(0.9, 0.1))
dev.off()

# Boxplot CANDIDATES
load("WGCNA_EXONS_NEUN.RData")

cool <- c("STK39","HYDIN","AIDA","KMT2A")
mat <- expAdj[rownames(expAdj) %in% cool,]
mat$Rows <- rownames(mat)
df <- melt(mat)

df$Species <- as.factor(do.call(rbind,strsplit(as.character(df$variable),"_"))[,2])

pdf("NeuN_Viz/Boxplot_Coolguys.pdf",width=7,height=3,useDingbats=FALSE)
tocompare <- list( c("Hsap", "PanTro"), c("Hsap", "RheMac"),c("PanTro", "RheMac"))
ggboxplot(df, "Species", "value",
color = "Species", palette =c("steelblue","grey60","green"),
add = "jitter", shape = "Species")+
facet_wrap(~Rows,ncol=4,nrow=1,scales="free")+
ylab("Adjusted Expression")+
theme_classic()+
theme(legend.position="none")+
xlab("")+
stat_compare_means(comparisons = tocompare, label = "p.signif")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()











