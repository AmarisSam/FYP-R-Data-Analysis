library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tibble)

gsea_2D <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\fgsearesults_2D.csv")
columns_to_delete <- c("pval", "padj","ES")
column_indices <- which(names(gsea_2D) %in% columns_to_delete)
gsea_2D <- gsea_2D[, -column_indices]

gsea_ULA <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\fgsearesults_ULA.csv")
columns_to_delete <- c("X","pval", "padj","ES")
column_indices <- which(names(gsea_ULA) %in% columns_to_delete)
gsea_ULA <- gsea_ULA[, -column_indices]

gsea_73 <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\fgsearesults_73.csv")
columns_to_delete <- c("X","pval", "padj","ES")
column_indices <- which(names(gsea_73) %in% columns_to_delete)
gsea_73 <- gsea_73[, -column_indices]

gsea_105 <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\fgsearesults_105.csv")
columns_to_delete <- c("X","pval", "padj","ES")
column_indices <- which(names(gsea_105) %in% columns_to_delete)
gsea_105 <- gsea_105[, -column_indices]

merged_df <- merge(gsea_2D, gsea_ULA, by = "pathway", all = TRUE)
colnames(merged_df)[colnames(merged_df) == "NES.x"]<-"2D"
colnames(merged_df)[colnames(merged_df) == "NES.y"]<-"ULA"
merged_df <- merge(merged_df, gsea_73, by = "pathway", all = TRUE)

#without 2D
merged_df <- merge(gsea_ULA, gsea_73, by = "pathway", all = TRUE)
colnames(merged_df)[colnames(merged_df) == "NES.x"]<-"ULA"
colnames(merged_df)[colnames(merged_df) == "NES.y"]<-"73.7"
merged_df <- merge(merged_df, gsea_105, by = "pathway", all = TRUE)
colnames(merged_df)[colnames(merged_df) == "NES"]<-"105.3"

row.names(merged_df)<-merged_df[,1]
columns_to_delete <- c("pathway")
column_indices <- which(names(merged_df) %in% columns_to_delete)
merged_df <- merged_df[, -column_indices]
merged_df<-na.omit(merged_df)

pheatmap(merged_df, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)

NES_scores <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\NES_scores.csv")
colnames(NES_scores)[colnames(NES_scores) == "X2D"]<-"2D"
colnames(NES_scores)[colnames(NES_scores) == "X73.7"]<-"73.7"
colnames(NES_scores)[colnames(NES_scores) == "X105.3"]<-"105.3"
NES_scores<-column_to_rownames(NES_scores,var="pathway")
pheatmap(NES_scores, breaks = seq(1.0, 2.0, length.out = 100), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, legend = TRUE, legend_title = "NES",
         main = "GSEA NES")

NES_scores <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\NES_scores_no2D.csv")
colnames(NES_scores)[colnames(NES_scores) == "X73.7"]<-"73.7"
colnames(NES_scores)[colnames(NES_scores) == "X105.3"]<-"105.3"
NES_scores<-column_to_rownames(NES_scores,var="X")
pheatmap(NES_scores, breaks = seq(1.0, 2.0, length.out = 100), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, legend = TRUE, legend_title = "NES",
         main = "GSEA NES")
