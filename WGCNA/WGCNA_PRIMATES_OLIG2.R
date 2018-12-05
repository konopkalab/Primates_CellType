# WGCNA

library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Load tables 
load("WGCNA_EXONS_OLIG2.RData")

tab <- expAdj
datExpr <- as.data.frame(t(tab));
names(datExpr) <- rownames(tab);
rownames(datExpr) <- names(tab);

# Loading Covariate
datTraits=pheno
datTraits$HumAge <- as.numeric(datTraits$HumAge) 

## Powers analysis
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize= 14000, networkType = "signed",RsquaredCut = 0.85) 
pdf("SoftThresholdingPower_signed.pdf")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######################################################################################################################
############################################ MODULE construction #####################################################
############################################    signed network   #####################################################
PWR=sft$powerEstimate
net = blockwiseModules(datExpr,
corType="bicor",
maxBlockSize = 15000,
networkType="signed",
minCoreKME = 0.4, 
minKMEtoStay = 0.5,
power=PWR, 
checkMissingData = TRUE,
minModuleSize=35,
nThreads=15,
TOMType = "signed",
TOMDenom = "mean",
deepSplit=4,
verbose=1,
mergeCutHeight=0.10,
reassignThreshold = 1e-10,
numericLabels=TRUE)

moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, "OLIG2_colors.txt",sep="\t",quote=F)
save(net,file="OLIG2_net.RData")

#KMEs
KMEs<-signedKME(datExpr, net$MEs,corFnc = "bicor")
kme=data.frame(rownames(tab), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"KME_OLIG2.txt",sep="\t",quote=F)

# Dat Trait Heatmap
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic)$eigengenes
MEsTemp = MEs0
MEsTemp$MEgrey=NULL
modTraitCor= cor(MEsTemp, datTraits,method="pearson")
write.table(modTraitCor,"modTraitCor_OLIG2.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"modTraitP_OLIG2.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor) 
par(mar = c(6, 8.5, 3, 3))
pdf("Heatmap_DatTraits.pdf",width=7,height=7)
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.4, zlim = c(-1,1), main = paste("Module Association"))
dev.off()

#Dendro Plot
pdf("NetworkDendrogram_OLIG2.pdf",width=10,height=3)
x=cor(datTraits,datExpr,method="pearson")
x[x > 0.25]=6
x[x < -0.25]=30
x[ x > -0.25 & x < 0.25] = 27
listLabels2 <- t(labels2colors(x))
rownames(listLabels2) = colnames(x) 
colnames(listLabels2) = rownames(x)
plotColors=as.data.frame(cbind(moduleColorsAutomatic,listLabels2));
colnames(plotColors)=c("Module",names(datTraits)) 
plotDendroAndColors(net$dendrograms[[1]], plotColors,colorHeight = 0.1, dendroLabels=FALSE,hang = 0.03, addGuide = TRUE, guideHang = 0.05 )
dev.off()

#GGplotInput
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic,softPower = PWR,impute = TRUE)$eigengenes
MEs0$Rows=colnames(tab)
write.table(MEs0, "Matrix_module_correlation.txt",sep="\t",quote=F)

# Adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "bicor")
moduleOutput <- data.frame(rownames(tab))
moduleOutput[,2]<- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
write.table(moduleOutput, "ModuleOutput_OLIG2.txt", sep="\t", quote=F)

# TO connectivity
TOM = TOMsimilarityFromExpr(datExpr, power= PWR,corType = "bicor",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(TOM)=rownames(TOM)=colnames(datExpr)
save(TOM,file="TOM_OLIG2_SIGNED.RData")
Connectivity=apply(TOM,1,sum)
save(Connectivity,file="Connectivity.RData")

# CytoScape output
dir.create("Cyto")
for(module in unique(moduleColorsAutomatic)){
inModule <- is.finite(match(moduleColorsAutomatic, module))
modTOM <- TOM[inModule, inModule]
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("Cyto/CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile=paste("Cyto/CytoNode",paste(module,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColorsAutomatic[inModule], nodeNames = names(datExpr)[inModule])
}

# WGCNA barplot output
library(ggplot2)
library(reshape2)
library(RColorBrewer)
df=melt(MEs0)
dir.create("Plot")
df$Rows=factor(df$Rows, levels = colnames(tab))
df=df[!df$variable == "MEgrey", ]
# Function to make the barplot and saving as pdf according to the module name
cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- c("orange","steelblue")
doPlot = function(sel_name) 
{
    df = subset(df, variable == sel_name)
    PLOT=ggplot(data=df, aes(x=Rows, y=value)) +
 geom_bar(aes(fill=Class),stat="identity",position = "dodge")+
 scale_y_continuous(limits=c(-1,+1))+
 theme_bw()+
 theme(strip.text.x = element_text(size=12, face="bold"),
 strip.background = element_rect(colour="black", fill="#CCCCFF"))+
 scale_fill_manual(values = cols)+
 theme(axis.title.x = element_blank(),
            axis.text.x  = element_text(face="bold", size=6,angle = 45, hjust = 1))+
 theme(axis.title.y = element_blank(),
            axis.text.y  = element_text(face="bold", size=6))
    print(PLOT)
    ggsave(sprintf("Plot/%s.pdf", sel_name))
 }

lapply(unique(df$variable), doPlot)




