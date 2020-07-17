library(SingleCellExperiment)
library(SC3)
library(scater)
# use SC3 to get the similarity matrix
# use the feature selected data
HSMM_expr_matrix <- read.table("SCC3.txt",header = TRUE,row.names = 1,sep = "\t")
HSMM_sample_sheet <- read.delim("CC3_cellname.txt")
gene_short_name <- read.delim("SCC3_gene_annotations.txt")
gene_short_name <- gene_short_name[,1]
HSMM_gene_annotation <- data.frame(gene_short_name,row.names = gene_short_name)
ann <- HSMM_sample_sheet[,1]
matrix <- as.matrix(HSMM_expr_matrix)
sce <- SingleCellExperiment(
  assays = list(
    counts = matrix,
    logcounts = log2(matrix+ 1)
  ), 
  colData = ann
)
rowData(sce)$feature_symbol <- rownames(sce)
#explore clustering of the data in the range of ks (the number of clusters) from 3 to 6
sce <- sc3(sce, ks = 3:6, biology = TRUE, gene_filter = FALSE)
M_3 = (sce@metadata$sc3$consensus$`3`$consensus)
M_4 = (sce@metadata$sc3$consensus$`4`$consensus)
M_5 = (sce@metadata$sc3$consensus$`5`$consensus)
M_6 = (sce@metadata$sc3$consensus$`6`$consensus)
# take average of all the clustering results to get cell-cell similarity matrix M
M = (M_3+M_4+M_5+M_6)/4
