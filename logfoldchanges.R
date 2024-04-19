library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(reshape)

#ULA
rawcountData_ULA <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\RHB_rawcount_tissuevsULA.csv", header = TRUE, sep = ",", row.names = "Column1")
condition <-rep(c("ULA","tissue"), times=3)
my_colData <-as.data.frame(condition)
rownames(my_colData) <-colnames(rawcountData_ULA)
dds <- DESeqDataSetFromMatrix(countData = rawcountData_ULA,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)

res <- results(dds, tidy = TRUE)
res <- res[order(res$log2FoldChange),]
my_annotation <- read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
res <- left_join(res, my_annotation, by = c("row" = "Gene.stable.ID"))
colnames(res)[colnames(res) == "row"]<-"ensembl_id"
ULA_res<-res

#737
rawcountData_73 <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\RHB_rawcount_tissuevsRGD73.csv", header = TRUE, sep = ",", row.names = "Column1")
condition <-rep(c("73.7","tissue"), times=3)
my_colData <-as.data.frame(condition)
rownames(my_colData) <-colnames(rawcountData_73)
dds <- DESeqDataSetFromMatrix(countData = rawcountData_73,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)
res <- res[order(res$log2FoldChange),]
res <- left_join(res, my_annotation, by = c("row" = "Gene.stable.ID"))
colnames(res)[colnames(res) == "row"]<-"ensembl_id"
RGD737_res<-res

#1053
rawcountData_105 <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\RHB_rawcount_tissuevsRGD105.csv", header = TRUE, sep = ",", row.names = "Column1")
condition <-rep(c("105.3","tissue"), times=3)
my_colData <-as.data.frame(condition)
rownames(my_colData) <-colnames(rawcountData_105)
dds <- DESeqDataSetFromMatrix(countData = rawcountData_105,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)
res <- res[order(res$log2FoldChange),]
res <- left_join(res, my_annotation, by = c("row" = "Gene.stable.ID"))
colnames(res)[colnames(res) == "row"]<-"ensembl_id"
RGD1053_res<-res

#2D
rawcountData_2D <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\RHB_rawcount_tissuevs2D.csv", header = TRUE, sep = ",", row.names = "Column1")
condition <-rep(c("2D","tissue"), times=3)
my_colData <-as.data.frame(condition)
rownames(my_colData) <-colnames(rawcountData_2D)
dds <- DESeqDataSetFromMatrix(countData = rawcountData_2D,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)
res <- res[order(res$log2FoldChange),]
res <- left_join(res, my_annotation, by = c("row" = "Gene.stable.ID"))
colnames(res)[colnames(res) == "row"]<-"ensembl_id"
res_2D<-res

#inflammatory response
inflam<-c("ABCA1","ABI1","ACVR1B","ACVR2A", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "ADGRE1", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "CXCL8", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELL", "SELENOS", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP")
ULA_res$InPathway <-ifelse(ULA_res$Gene.name%in%inflam, "Pathway", "Other")
ULA_res_inflam <- ULA_res[ULA_res$InPathway == "Pathway", ]
write.csv(ULA_res_inflam,file ="USA_res_inflam.csv")
inflam_log2FC<-read.csv("USA_res_inflam_p005.csv",header = TRUE, sep = ",")
inflam_log2FC <- merge(inflam_log2FC,inflam_log2FC_a, by = "Gene.name",all = TRUE)

RGD737_res$InPathway <-ifelse(RGD737_res$Gene.name%in%inflam, "Pathway", "Other")
RGD737_res_inflam <- RGD737_res[RGD737_res$InPathway == "Pathway", ]
write.csv(RGD737_res_inflam,file ="RGD737_res_inflam.csv")
inflam_log2FC_a<-read.csv("RGD737_res_inflam_p005.csv",header = TRUE, sep = ",")

RGD1053_res$InPathway <-ifelse(RGD1053_res$Gene.name%in%inflam, "Pathway", "Other")
RGD1053_res_inflam <- RGD1053_res[RGD1053_res$InPathway == "Pathway", ]
write.csv(RGD1053_res_inflam,file ="RGD1053_res_inflam.csv")
inflam_log2FC_a<-read.csv("RGD1053_res_inflam_p005.csv",header = TRUE, sep = ",")

res_2D$InPathway <-ifelse(res_2D$Gene.name%in%inflam, "Pathway", "Other")
res_2D_inflam <- res_2D[res_2D$InPathway == "Pathway", ]
write.csv(res_2D_inflam,file ="res_2D_inflam.csv")
inflam_log2FC_a<-read.csv("res_2D_inflam_p005.csv",header = TRUE, sep = ",")

colnames(inflam_log2FC)[colnames(inflam_log2FC) == "X2D"]<-"2D"
row.names(inflam_log2FC)<-inflam_log2FC[,1]
columns_to_delete <- c("Gene.name")
column_indices <- which(names(inflam_log2FC) %in% columns_to_delete)
inflam_log2FC <- inflam_log2FC[, -column_indices]
inflam_log2FC<-na.omit(inflam_log2FC)

df<-melt(inflam_log2FC)
row.names(inflam_log2FC)<-inflam_log2FC[,1]
colnames(df) <- c("Gene.name", "condition", "log2FC")
ggplot(df, aes(x = condition, y = Gene.name, fill = log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF",
                                     mid = "#FFFFCC",
                                     high = "#FF0000") + 
  ggtitle("Expression of genes involved in inflammatory response")

#myogenesis
myogen<-c("ABLIM1","ACHE","ACSL1","ACTA1","ACTC1","ACTN2","ACTN3","ADAM12","ADCY9","AEBP1","AGL","AGRN","AK1","AKT2","ANKRD2","APLNR","APOD","APP","ATP2A1","ATP6AP1","BAG1","BDKRB2","BHLHE40","BIN1","CACNA1H","CACNG1","CAMK2B","CASQ1","CASQ2","CAV3","CD36","CDH13","CDKN1A","CFD","CHRNA1","CHRNB1","CHRNG","CKB","CKM","CKMT2","CLU","CNN3","COL15A1","COL1A1","COL3A1","COL4A2","COL6A2","COL6A3","COX6A2","COX7A1","CRAT","CRYAB","CSRP3","CTF1","DAPK2","DES","DMD","DMPK","DTNA","EFS","EIF4A2","ENO3","EPHB3","ERBB3","FABP3","FDPS","FGF2","FHL1","FKBP1B","FLII","FOXO4","FST","FXYD1","GAA","GABARAPL2","GADD45B","GJA5","GNAO1","GPX3","GSN","HBEGF","HDAC5","HRC","HSPB2","HSPB8","IFRD1","IGF1","IGFBP3","IGFBP7","ITGA7","ITGB1","ITGB4","ITGB5","KCNH1","KCNH2","KIFC3","KLF5","LAMA2","LARGE1","LDB3","LPIN1","LSP1","MAPK12","MAPRE3","MB","MEF2A","MEF2C","MEF2D","MRAS","MYBPC3","MYBPH","MYF6","MYH1","MYH11","MYH2","MYH3","MYH4","MYH7","MYH8","MYH9","MYL1","MYL2","MYL3","MYL4","MYL6B","MYL7","MYLK","MYL11","MYO1C","MYOG","MYOM1","MYOM2","MYOZ1","NAV2","NCAM1","NOS1","NOTCH1","NQO1","OCEL1","PC","PDE4DIP","PDLIM7","PFKM","PGAM2","PICK1","PKIA","PLXNB2","PPFIA4","PPP1R3C","PRNP","PSEN2","PTGIS","PTP4A3","PVALB","PYGM","RB1","REEP1","RIT1","RYR1","SCD","SCHIP1","SGCA","SGCD","SGCG","SH2B1","SH3BGR","SIRT2","SLC6A8","SLN","SMTN","SOD3","SORBS1","SORBS3","SPARC","SPDEF","SPEG","SPHK1","SPTAN1","SSPN","DENND2B","STC2","SVIL","SYNGR2","TAGLN","TCAP","TEAD4","TGFB1","TNNC1","TNNC2","TNNI1","TNNI2","TNNT1","TNNT2","TNNT3","TPD52L1","TPM2","TPM3","TSC2","VIPR1","WWTR1")
ULA_res$InPathway <-ifelse(ULA_res$Gene.name%in%myogen, "Pathway", "Other")
ULA_res_myogen <- ULA_res[ULA_res$InPathway == "Pathway", ]
write.csv(ULA_res_myogen,file ="USA_res_myogen.csv")
myogen_log2FC<-read.csv("USA_res_myogen_p005.csv",header = TRUE, sep = ",")
myogen_log2FC <- merge(myogen_log2FC,myogen_log2FC_a, by = "Gene.name",all = TRUE)

RGD737_res$InPathway <-ifelse(RGD737_res$Gene.name%in%myogen, "Pathway", "Other")
RGD737_res_myogen <- RGD737_res[RGD737_res$InPathway == "Pathway", ]
write.csv(RGD737_res_myogen,file ="RGD737_res_myogen.csv")
myogen_log2FC_a<-read.csv("RGD737_res_myogen_p005.csv",header = TRUE, sep = ",")

RGD1053_res$InPathway <-ifelse(RGD1053_res$Gene.name%in%myogen, "Pathway", "Other")
RGD1053_res_myogen <- RGD1053_res[RGD1053_res$InPathway == "Pathway", ]
write.csv(RGD1053_res_myogen,file ="RGD1053_res_myogen.csv")
myogen_log2FC_a<-read.csv("RGD1053_res_myogen_p005.csv",header = TRUE, sep = ",")

res_2D$InPathway <-ifelse(res_2D$Gene.name%in%myogen, "Pathway", "Other")
res_2D_myogen <- res_2D[res_2D$InPathway == "Pathway", ]
write.csv(res_2D_myogen,file ="res_2D_myogen.csv")
myogen_log2FC_a<-read.csv("res_2D_myogen_p005.csv",header = TRUE, sep = ",")

colnames(myogen_log2FC)[colnames(myogen_log2FC) == "X2D"]<-"2D"
row.names(myogen_log2FC)<-myogen_log2FC[,1]
columns_to_delete <- c("Gene.name")
column_indices <- which(names(myogen_log2FC) %in% columns_to_delete)
myogen_log2FC <- myogen_log2FC[, -column_indices]
myogen_log2FC<-na.omit(myogen_log2FC)

df_myo<-melt(myogen_log2FC)
row.names(myogen_log2FC)<-myogen_log2FC[,1]
colnames(df_myo) <- c("Gene.name", "condition", "log2FC")
ggplot(df_myo, aes(x = condition, y = Gene.name, fill = log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF",
                                     mid = "#FFFFCC",
                                     high = "#FF0000") + 
  ggtitle("Expression of genes involved in myogenesis")

#TNFA_SIGNALING_VIA_NFKB
TNFA<-c("ABCA1","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","ACKR3","CCN1","RIGI","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PMEPA1","PNRC1","PLPP3","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")
ULA_res$InPathway <-ifelse(ULA_res$Gene.name%in%TNFA, "Pathway", "Other")
ULA_res_TNFA <- ULA_res[ULA_res$InPathway == "Pathway", ]
write.csv(ULA_res_TNFA,file ="USA_res_TNFA.csv")
TNFA_log2FC<-read.csv("USA_res_TNFA_p005.csv",header = TRUE, sep = ",")
TNFA_log2FC <- merge(TNFA_log2FC,TNFA_log2FC_a, by = "Gene.name",all = TRUE)

RGD737_res$InPathway <-ifelse(RGD737_res$Gene.name%in%TNFA, "Pathway", "Other")
RGD737_res_TNFA <- RGD737_res[RGD737_res$InPathway == "Pathway", ]
write.csv(RGD737_res_TNFA,file ="RGD737_res_TNFA.csv")
TNFA_log2FC_a<-read.csv("RGD737_res_TNFA_p005.csv",header = TRUE, sep = ",")

RGD1053_res$InPathway <-ifelse(RGD1053_res$Gene.name%in%TNFA, "Pathway", "Other")
RGD1053_res_TNFA <- RGD1053_res[RGD1053_res$InPathway == "Pathway", ]
write.csv(RGD1053_res_TNFA,file ="RGD1053_res_TNFA.csv")
TNFA_log2FC_a<-read.csv("RGD1053_res_TNFA_p005.csv",header = TRUE, sep = ",")

res_2D$InPathway <-ifelse(res_2D$Gene.name%in%TNFA, "Pathway", "Other")
res_2D_TNFA <- res_2D[res_2D$InPathway == "Pathway", ]
write.csv(res_2D_TNFA,file ="res_2D_TNFA.csv")
TNFA_log2FC_a<-read.csv("res_2D_TNFA_p005.csv",header = TRUE, sep = ",")

colnames(TNFA_log2FC)[colnames(TNFA_log2FC) == "X2D"]<-"2D"
TNFA_log2FC<-na.omit(TNFA_log2FC)

df_TNFA<-melt(TNFA_log2FC)
row.names(TNFA_log2FC)<-TNFA_log2FC[,1]
colnames(df_TNFA) <- c("Gene.name", "condition", "log2FC")
ggplot(df_TNFA, aes(x = condition, y = Gene.name, fill = log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF",
                                     mid = "#FFFFCC",
                                     high = "#FF0000") + 
  ggtitle("Expression of genes involved in TNFA signalling via NFKB")

#hypoxia
hypo<-c("ADM","ADORA2B","AK4","AKAP12","ALDOA","ALDOB","ALDOC","AMPD3","ANGPTL4","ANKZF1","ANXA2","ATF3","ATP7A","B3GALT6","B4GALNT2","BCAN","BCL2","BGN","BHLHE40","BNIP3L","BRS3","BTG1","CA12","CASP6","CAV1","CCNG2","NOCT","CDKN1A","CDKN1B","CDKN1C","CHST2","CHST3","CITED2","COL5A1","CP","CSRP2","CCN2","CXCR4","ACKR3","CCN1","DCN","DDIT3","DDIT4","DPYSL4","DTNA","DUSP1","EDN2","EFNA1","EFNA3","EGFR","ENO1","ENO2","ENO3","ERO1A","ERRFI1","ETS1","EXT1","F3","FAM162A","FBP1","FOS","FOSL2","FOXO3","GAA","GALK1","GAPDH","GAPDHS","GBE1","GCK","GCNT2","GLRX","GPC1","GPC3","GPC4","GPI","GRHPR","GYS1","HAS1","HDLBP","HEXA","HK1","HK2","HMOX1","HOXB9","HS3ST1","HSPA5","IDS","IER3","IGFBP1","IGFBP3","IL6","ILVBL","INHA","IRS2","ISG20","JMJD6","JUN","KDELR3","KDM3A","KIF5A","KLF6","KLF7","KLHL24","LALBA","LARGE1","LDHA","LDHC","LOX","LXN","MAFF","MAP3K1","MIF","MT1E","MT2A","MXI1","MYH9","NAGK","NCAN","NDRG1","NDST1","NDST2","NEDD4L","NFIL3","NR3C1","P4HA1","P4HA2","PAM","PCK1","PDGFB","PDK1","PDK3","PFKFB3","PFKL","PFKP","PGAM2","PGF","PGK1","PGM1","PGM2","PHKG1","PIM1","PKLR","PKP1","PLAC8","PLAUR","PLIN2","PNRC1","PPARGC1A","PPFIA4","PPP1R15A","PPP1R3C","PRDX5","PRKCA","CAVIN3","CAVIN1","PYGM","RBPJ","RORA","RRAGD","S100A4","SAP30","SCARB1","SDC2","SDC3","SDC4","SELENBP1","SERPINE1","SIAH2","SLC25A1","SLC2A1","SLC2A3","SLC2A5","SLC37A4","SLC6A6","SRPX","STBD1","STC1","STC2","SULT2B1","TES","TGFB3","TGFBI","TGM2","TIPARP","TKTL1","TMEM45A","TNFAIP3","TPBG","TPD52","TPI1","TPST2","UGP2","VEGFA","VHL","VLDLR","CCN5","WSB1","XPNPEP1","ZFP36","ZNF292")
ULA_res$InPathway <-ifelse(ULA_res$Gene.name%in%hypo, "Pathway", "Other")
ULA_res_hypo <- ULA_res[ULA_res$InPathway == "Pathway", ]
write.csv(ULA_res_hypo,file ="USA_res_hypo.csv")
hypo_log2FC<-read.csv("USA_res_hypo_p005.csv",header = TRUE, sep = ",")
hypo_log2FC <- merge(hypo_log2FC,hypo_log2FC_a, by = "Gene.name",all = TRUE)

RGD737_res$InPathway <-ifelse(RGD737_res$Gene.name%in%hypo, "Pathway", "Other")
RGD737_res_hypo <- RGD737_res[RGD737_res$InPathway == "Pathway", ]
write.csv(RGD737_res_hypo,file ="RGD737_res_hypo.csv")
hypo_log2FC_a<-read.csv("RGD737_res_hypo_p005.csv",header = TRUE, sep = ",")

RGD1053_res$InPathway <-ifelse(RGD1053_res$Gene.name%in%hypo, "Pathway", "Other")
RGD1053_res_hypo <- RGD1053_res[RGD1053_res$InPathway == "Pathway", ]
write.csv(RGD1053_res_hypo,file ="RGD1053_res_hypo.csv")
hypo_log2FC_a<-read.csv("RGD1053_res_hypo_p005.csv",header = TRUE, sep = ",")

res_2D$InPathway <-ifelse(res_2D$Gene.name%in%hypo, "Pathway", "Other")
res_2D_hypo <- res_2D[res_2D$InPathway == "Pathway", ]
write.csv(res_2D_hypo,file ="res_2D_hypo.csv")
hypo_log2FC_a<-read.csv("res_2D_hypo_p005.csv",header = TRUE, sep = ",")

colnames(hypo_log2FC)[colnames(hypo_log2FC) == "X2D"]<-"2D"
hypo_log2FC<-na.omit(hypo_log2FC)

df_hypo<-melt(hypo_log2FC)
row.names(hypo_log2FC)<-hypo_log2FC[,1]
colnames(df_hypo) <- c("Gene.name", "condition", "log2FC")
ggplot(df_hypo, aes(x = condition, y = Gene.name, fill = log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF",
                                     mid = "#FFFFCC",
                                     high = "#FF0000") + 
  ggtitle("Expression of genes involved in hypoxia")

#muscle cell state gene signatures
musc<-read.csv("Musclestategenesignature.csv",header = FALSE, sep = ",")
musc<-musc[["V1"]]
ULA_res$InPathway <-ifelse(ULA_res$Gene.name%in%musc, "Pathway", "Other")
ULA_res_musc<- ULA_res[ULA_res$InPathway == "Pathway", ]
write.csv(ULA_res_musc,file ="USA_res_musc.csv")
musc_log2FC<-read.csv("USA_res_musc_p005.csv",header = TRUE, sep = ",")
musc_log2FC <- merge(musc_log2FC,musc_log2FC_a, by = "Gene.name",all = TRUE)

RGD737_res$InPathway <-ifelse(RGD737_res$Gene.name%in%musc, "Pathway", "Other")
RGD737_res_musc <- RGD737_res[RGD737_res$InPathway == "Pathway", ]
write.csv(RGD737_res_musc,file ="RGD737_res_musc.csv")
musc_log2FC_a<-read.csv("RGD737_res_musc_p005.csv",header = TRUE, sep = ",")

RGD1053_res$InPathway <-ifelse(RGD1053_res$Gene.name%in%musc, "Pathway", "Other")
RGD1053_res_musc <- RGD1053_res[RGD1053_res$InPathway == "Pathway", ]
write.csv(RGD1053_res_musc,file ="RGD1053_res_musc.csv")
musc_log2FC_a<-read.csv("RGD1053_res_musc_p005.csv",header = TRUE, sep = ",")

res_2D$InPathway <-ifelse(res_2D$Gene.name%in%musc, "Pathway", "Other")
res_2D_musc <- res_2D[res_2D$InPathway == "Pathway", ]
write.csv(res_2D_musc,file ="res_2D_musc.csv")
musc_log2FC_a<-read.csv("res_2D_musc_p005.csv",header = TRUE, sep = ",")

colnames(musc_log2FC)[colnames(musc_log2FC) == "X2D"]<-"2D"
musc_log2FC<-na.omit(musc_log2FC)

musc_log2FC_005<-read.csv("musc_res_p005.csv",header = TRUE, sep = ",")
musc_log2FC_005<-musc_log2FC_005[,2:5]
colnames(musc_log2FC_005)[colnames(musc_log2FC_005) == "X2D"]<-"2D"
pheatmap(musc_log2FC_005, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, breaks = seq(-5, 5, length.out = 100))

df_musc<-melt(musc_log2FC)
row.names(musc_log2FC)<-musc_log2FC[,1]
colnames(df_musc) <- c("Gene.name", "condition", "log2FC")
ggplot(df_musc, aes(x = condition, y = Gene.name, fill = log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF",
                                     mid = "#FFFFCC",
                                     high = "#FF0000") + 
  ggtitle("Expression of genes")
