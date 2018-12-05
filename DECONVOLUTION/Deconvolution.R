suppressPackageStartupMessages(library(xbioc))
suppressPackageStartupMessages(library(MuSiC))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))

# Load NeuN data and create Expset
expr <- read.table("HCR_brain.txt",header=T,sep="\t",row.names=1)
expr$Symbol <- mapIds(EnsDb.Hsapiens.v86,keys=rownames(expr), column="GENENAME",keytype="GENEID",multiVals="first")
expr <- na.omit(expr)
rownames(expr) <- expr$Symbol
expr$Symbol <- NULL
expr <- as.matrix(log2(cpm(expr)+1))

pheno <- data.frame(Sample=colnames(expr),Species = c(rep("Hsap",5),rep("PanTro",5),rep("RheMac",5)))
rownames(pheno) <- pheno$Sample
phenoData <- new("AnnotatedDataFrame",data=pheno)
bertoSet <- ExpressionSet(assayData=expr, phenoData=phenoData)

# Load Allen MTG scRNA-seq data and create ExpSet
load("Allen_MTG.RData")
phenoSC <- allenMetaData
rownames(phenoSC) <- phenoSC$Sample_ID
phenoDataSC <- new("AnnotatedDataFrame",data=phenoSC)
exprSC <- as.matrix(allenData)
exprSC <- log2(cpm(exprSC)+1)
scSet <- ExpressionSet(assayData=exprSC, phenoData=phenoDataSC)

# Deconvolution
# NeuN
BERTO_Deconv = music_prop(bulk.eset = bertoSet, sc.eset = scSet,clusters = 'Cell', samples = 'Sample_ID',iter.max = 1000,nu = 1e-10)
save(BERTO_Deconv,file = "BERTO_Deconv.RData")
BERTO_Prob <- as.data.frame(BERTO_Deconv$Est.prop.weighted)
BERTO_Prob$Rows <- rownames(BERTO_Prob)
tmp <- melt(BERTO_Prob)
dfBERTO <- merge(tmp,pheno,by.x="Rows",by.y="row.names",all=T)

pdf("DECONVOLUTION_BERTO.pdf",width=5,height=4,useDingbats=FALSE)
ggboxplot(dfBERTO, "variable", "value",
color = "Species", palette =c("steelblue", "grey60","green"))+
ggtitle("Cell Type Proportion: BertoEtAl DFC")+
xlab("")+ 
ylab("Proportion")+
theme(legend.position="right")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
ylim(0,1)
dev.off()




