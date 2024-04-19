1) R Scripts used in RNA bulk sequencing
  a) DE analysis
   - files:
    1. DESEQ2_tissuevs2D.R
    2. DESEQ2_tissuevsRGD105.R
    3. DESEQ2_tissuevsRGD73.R
   4. DESEQ2_tissuevsULA.R
  - DE analyses was done for all conditions using DESEQ2 & tumour tissue sample as reference
  - output (Figures): 2-5, 6-9
  - output: Differential Expression Results, GSEA Data

  b) NES Heatmap Plot
   - files: Heatmap_Plot_NESscores.R
   - input: NES scores.csv
       - location: Supplementary Data/GSEA Data
   - output (Figure): 10

  c) Pathway-specific Gene Expression Heatmap Plot
   - files: logfoldchanges.R
   - output (Figure): 11-15

2) Differential Expression Results
- output: List of differentially expressed genes (p = 0.05) in conditions using tumour tissue sample as reference
  - location: Supplementary Data/Differential Expression Data
- files: 
  1. 2D_vs_tissue_padj_cutoff.csv
  2. ULA_vs_tissue_padj_cutoff.csv
  3. 73.7_vs_tissue_padj_cutoff.csv
  4. 105.3_vs_tissue_padj_cutoff.csv

3) GSEA Data
- input: List of ordered ranked of only protein coding genes for GSEA from DESEQ2-processed data
  - location: Supplementary Data/GSEA Data
- files: 
  1. 2D_vs_tissue_rank.csv
  2. ULA_vs_tissue_rank.csv
  3. 73.7_vs_tissue_rank.csv
  4. 105.3_vs_tissue_rank.csv

- output: NES score and p-adjusted value of hallmark pathways
  - location: Supplementary Data/GSEA Data
- files: 
  1. fgsearesults_2D.csv
  2. fgsearesults_ULA.csv
  3. fgsearesults_73.7.csv
  4. fgsearesults_105.3.csv


