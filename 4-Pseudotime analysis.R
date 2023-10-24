library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(CytoTRACE)
rm(list=ls())
maize_sub <- readRDS("maize_sub.rds")
## Extract a subset of cells
Cells_alse <- subset(maize_sub@meta.data, seurat_clusters==c("0","1","3","4","5","11")
maize_alse <- subset(maize_sub, cells=row.names(Cells_alse))
Cells_alse_dap7 <- subset(maize_alse@meta.data, day==c("7d") 
maize_alse_dap7 <- subset(maize_alse, cells=row.names(Cells_alse_dap7))  
##run monocle2
data <- as(as.matrix(maize_alse_dap7@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
###
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=64, relative_expr = TRUE)
###
##Highly variable genes selected using monocle2
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)

#dimensionality
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#order
mycds <- orderCells(mycds)
##plot
plot <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
##Extract monocle2 sort coordinates
monocle_dim<-plot$data[,c("data_dim_1","data_dim_2")]
rownames(monocle_dim)<-rownames(plot$data)
###run CytoTRACE
maize_alse_dap7_expr <- maize_alse_dap7@assays$RNA@counts                       
maize_alse_dap7_results <- CytoTRACE(maize_alse_dap7_expr, ncores = 64)
plotCytoTRACE(maize_alse_dap7_results, emb = monocle_dim)       
##
saveRDS(maize_alse_dap7,"maize_alse_dap7.rds")
                    
