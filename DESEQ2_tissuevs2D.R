#tissue vs 2D
library("DESeq2")
library(ggplot2)
library(pheatmap)
library(fgsea)
library(dplyr)
library(tibble)
library(tidyverse)

rawcountData_2D <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\RHB_rawcount_tissuevs2D.csv", header = TRUE, sep = ",", row.names = "Column1")
condition <-rep(c("2D","tissue"), times=3)
my_colData <-as.data.frame(condition)
rownames(my_colData) <-colnames(rawcountData_2D)

dds <- DESeqDataSetFromMatrix(countData = rawcountData_2D,
                              colData = my_colData,
                              design = ~condition)


dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = T)
res <- results(dds, tidy = TRUE)
res <- res[order(res$log2FoldChange),]

annotation <- read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
annotated_data <- right_join(annotation, normalized_counts, by = c("Gene.stable.ID" = "ensembl_id"))

#ez volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
}
vsd <- vst(dds, blind = TRUE)
plotDists(vsd)


generate_DE_results <- function (dds, comparisons, padjcutoff = 0.05, log2cutoff = 1.0, cpmcutoff = 2) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "ensembl_id")
  # read in the annotation and append it to the data
  my_annotation <- read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "Gene.stable.ID"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  
  # write output to files
  write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
                     nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
                     nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
                     nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
                     "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}


plot_volcano <- function (res, padj_cutoff, nlabel = 10, label.by = "padj", title){
  # assign significance to results based on padj
  res <- mutate(res, significance=ifelse(res$padj<padj_cutoff, paste0("padj < ", padj_cutoff), paste0("padj > ", padj_cutoff)))
  res = res[!is.na(res$significance),]
  significant_genes <- res %>% filter(significance == paste0("padj < ", padj_cutoff))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    top_genes <- significant_genes %>% arrange(padj) %>% head(nlabel)
    bottom_genes <- significant_genes %>% filter (log2FoldChange < 0) %>% arrange(padj) %>% head (nlabel)
  } else if (label.by == "log2FoldChange") {
    top_genes <- head(arrange(significant_genes, desc(log2FoldChange)),nlabel)
    bottom_genes <- head(arrange(significant_genes, log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
  ggplot(res, aes(log2FoldChange, -log(padj))) +
    geom_point(aes(col=significance)) + ggtitle(title) +
    scale_color_manual(values=c("red", "black")) + 
    ggrepel::geom_text_repel(data=top_genes, aes(label=head(Gene.name,nlabel)), size = 3)+
    ggrepel::geom_text_repel(data=bottom_genes, aes(label=head(Gene.name,nlabel)), color = "#619CFF", size = 3)+
    labs ( x = "Log2FoldChange", y = "-(Log normalized p-value)")+
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme_minimal()
}

tissuevs2D_output <- generate_DE_results(dds, c("2D","tissue" ))
res <-read.csv("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\GSEA_Rfiles\\2D_vs_tissue_allgenes.csv", header = T)
plot_volcano(res, 0.05, nlabel = 15, label.by = "padj","2D vs tissue")

#gsea
library(fgsea)
# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))
ranked_list <- read.table("C:\\Users\\SamGa\\Downloads\\FYP Stuff\\GSEA_Rfiles\\2D_vs_tissue_rank.rnk", header = T, stringsAsFactors = F)

prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$log2FoldChange, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

ranked_list <- prepare_ranked_list(ranked_list)

fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = ranked_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm= 1000)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}

waterfall_plot(fgsea_results, "Hallmark pathways altered in 2D culture")
