dir.create("DGE_DATA_INPUTS")

Hspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
Cspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_CHIMP_SPECIFIC.txt")
Mspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")
PrimateAcc <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_PRIMATE_ACCELERATED.txt")
AllStat <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE.txt")

save(Hspecific,Cspecific,Mspecific,PrimateAcc,AllStat,file="DGE_DATA_INPUTS/DGE_Primates_NeuN.RData")

Hspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
Cspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_CHIMP_SPECIFIC.txt")
Mspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")
PrimateAcc <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_PRIMATE_ACCELERATED.txt")
AllStat <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE.txt")

save(Hspecific,Cspecific,Mspecific,PrimateAcc,AllStat,file="DGE_DATA_INPUTS/DGE_Primates_OLIG2.RData")

# Make input for enrichment
Hspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
Cspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_CHIMP_SPECIFIC.txt")
Mspecific <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")
PrimateAcc <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_PRIMATE_ACCELERATED.txt")

h <- data.frame(Gene = Hspecific$Gene,Class = ifelse(Hspecific$Estimate_HsapVsPanTro > 0,"Hsap_UpReg","Hsap_DownReg"))
c <- data.frame(Gene = Cspecific$Gene,Class = ifelse(Cspecific$Estimate_HsapVsPanTro < 0,"PanTro_UpReg","PanTro_DownReg"))
r <- data.frame(Gene = Mspecific$Gene,Class = ifelse(Mspecific$Estimate_HsapVsRheMac < 0,"RheMac_UpReg","RheMac_DownReg"))
Pacc <- data.frame(Gene = PrimateAcc$Gene,Class = ifelse(PrimateAcc$Estimate_HsapVsPanTro > 0,"Primate_Accel","Primate_Decel"))

tmp <- rbind(h,c,r)
write.table(tmp, "DGE_DATA_INPUTS/NeuN_Primates_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(h, "DGE_DATA_INPUTS/NeuN_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(c, "DGE_DATA_INPUTS/NeuN_PanTro_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(r, "DGE_DATA_INPUTS/NeuN_RheMac_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(Pacc, "DGE_DATA_INPUTS/NeuN_Primate_Acceleration.txt",sep="\t",quote=F,row.names=F)

# Make input for enrichment
Hspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
Cspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_CHIMP_SPECIFIC.txt")
Mspecific <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")
PrimateAcc <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_PRIMATE_ACCELERATED.txt")

h <- data.frame(Gene = Hspecific$Gene,Class = ifelse(Hspecific$Estimate_HsapVsPanTro > 0,"Hsap_UpReg","Hsap_DownReg"))
c <- data.frame(Gene = Cspecific$Gene,Class = ifelse(Cspecific$Estimate_HsapVsPanTro < 0,"PanTro_UpReg","PanTro_DownReg"))
r <- data.frame(Gene = Mspecific$Gene,Class = ifelse(Mspecific$Estimate_HsapVsRheMac < 0,"RheMac_UpReg","RheMac_DownReg"))
Pacc <- data.frame(Gene = PrimateAcc$Gene,Class = ifelse(PrimateAcc$Estimate_HsapVsPanTro > 0,"Primate_Accel","Primate_Decel"))

tmp <- rbind(h,c,r)
write.table(tmp, "DGE_DATA_INPUTS/OLIG2_Primates_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(h, "DGE_DATA_INPUTS/OLIG2_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(c, "DGE_DATA_INPUTS/OLIG2_PanTro_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(r, "DGE_DATA_INPUTS/OLIG2_RheMac_Dge_Input.txt",sep="\t",quote=F,row.names=F)
write.table(Pacc, "DGE_DATA_INPUTS/OLIG2_Primate_Acceleration.txt",sep="\t",quote=F,row.names=F)

# Make input for enrichment with only shared
n <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
o <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt")


nUp <- n[n$Estimate_HsapVsPanTro > 0,]
nDown <- n[n$Estimate_HsapVsPanTro < 0,]
oUp <- o[o$Estimate_HsapVsPanTro > 0,]
oDown <- o[o$Estimate_HsapVsPanTro < 0,]

intUp <- intersect(nUp$Gene, oUp$Gene)
intDown <- intersect(nDown$Gene, oDown$Gene)

nUp <- nUp[!(nUp$Gene %in% intUp),]
oUp <- oUp[!(oUp$Gene %in% intUp),]

nDown <- nDown[!(nDown$Gene %in% intDown),]
oDown <- oDown[!(oDown$Gene %in% intDown),]


dfN <- rbind(nUp,nDown)
dfO <- rbind(oUp,oDown)

h <- data.frame(Gene = dfN$Gene,Class = ifelse(dfN$Estimate_HsapVsPanTro > 0,"Hsap_Spec_UpReg","Hsap_Spec_DownReg"))
write.table(h, "DGE_DATA_INPUTS/NeuN_Specific_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)
h <- data.frame(Gene = dfO$Gene,Class = ifelse(dfO$Estimate_HsapVsPanTro > 0,"Hsap_Spec_UpReg","Hsap_Spec_DownReg"))
write.table(h, "DGE_DATA_INPUTS/OLIG2_Specific_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)

Common <- data.frame(Gene = c(intUp,intDown),Class=c(rep("Hsap_Common_UpReg",length(intUp)),rep("Hsap_Common_DownReg",length(intDown))))
write.table(Common, "DGE_DATA_INPUTS/CommonDEGs_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)

## Both 
n <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_HUMAN_SPECIFIC.txt")
o <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_HUMAN_SPECIFIC.txt")

tmp <- data.frame(Gene = c(as.character(n$Gene), as.character(o$Gene)),Class = c(rep("NeuN",nrow(n)),rep("OLIG2",nrow(o))))
write.table(tmp, "DGE_DATA_INPUTS/BOTH_Hsap_Dge_Input.txt",sep="\t",quote=F,row.names=F)

n <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_CHIMP_SPECIFIC.txt")
o <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_CHIMP_SPECIFIC.txt")

tmp <- data.frame(Gene = c(as.character(n$Gene), as.character(o$Gene)),Class = c(rep("NeuN",nrow(n)),rep("OLIG2",nrow(o))))
write.table(tmp, "DGE_DATA_INPUTS/BOTH_PanTro_Dge_Input.txt",sep="\t",quote=F,row.names=F)

n <- read.table("LM/OUTPUTS_EXONS_NEUN/NEUN_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")
o <- read.table("LM/OUTPUTS_EXONS_OLIG2/OLIG2_PRIMATES_DGE_MACAQUE_SPECIFIC.txt")

tmp <- data.frame(Gene = c(as.character(n$Gene), as.character(o$Gene)),Class = c(rep("NeuN",nrow(n)),rep("OLIG2",nrow(o))))
write.table(tmp, "DGE_DATA_INPUTS/BOTH_RheMac_Dge_Input.txt",sep="\t",quote=F,row.names=F)









