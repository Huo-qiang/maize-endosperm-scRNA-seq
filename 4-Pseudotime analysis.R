library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(CytoTRACE)
rm(list=ls())
maize_sub <- readRDS("maize_sub.rds")
##提取细胞子集
Cells_alse <- subset(maize_sub@meta.data, seurat_clusters==c("0","1","3","4","5","11")
maize_alse <- subset(maize_sub, cells=row.names(Cells_alse))
Cells_alse_dap7 <- subset(maize_alse@meta.data, day==c("7d") 
maize_alse_dap7 <- subset(maize_alse, cells=row.names(Cells_alse_dap7))  
##monocle2拟时分析
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
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
##Cluster轨迹分布图
plot <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
##提取monocle2排序坐标
monocle_dim<-plot$data[,c("data_dim_1","data_dim_2")]
rownames(monocle_dim)<-rownames(plot$data)
###CytoTRACE拟时分析
maize_alse_dap7_expr <- maize_alse_dap7@assays$RNA@counts                       
maize_alse_dap7_results <- CytoTRACE(maize_alse_dap7_expr, ncores = 64)
plotCytoTRACE(maize_alse_dap7_results, emb = monocle_dim)       
##
saveRDS(maize_alse_dap7,"maize_alse_dap7.rds")
                    
