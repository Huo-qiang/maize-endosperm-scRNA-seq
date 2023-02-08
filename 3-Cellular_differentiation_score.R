library(Seurat)
library(ggplot2)
library(UCell)
library(viridis)
##读取6dap和8dap的marker基因
maize_sub <- readRDS("maize_sub.rds")
dap6markergene<-read.csv("day6marker.csv")
dap8markergene<-read.csv("day5marker.csv")
##标记基因与单细胞检测到的基因取交集
dap6markergene<-intersect(rownames(submaize),dap6markergene$day6)
dap8markergene<-intersect(rownames(submaize),dap6markergene$day8)
# 构建UCell需要输入的gene sets
markers <- list()
markers$dap6<-dap6markergene
markers$dap8<-dap8markergene
###
maize_sub_ucell <- AddModuleScore_UCell(maize_sub, features = markers)
cell_information <- maize_sub_ucell@meta.data
###
cell_information$differentiation_score<-cell_information$dap8_UCell/cell_information$dap6_UCell
cell_information$differentiation_score<- log2((1+cell_information[,22]))
maize_sub_ucell@meta.data<- cell_information
##
FeaturePlot(maize_sub_ucell, reduction = "umap", features = "differentiation_score", 
            ncol = 1, order = T,
            min.cutoff = "q03", max.cutoff = "q99", pt.size = 0.1)+scale_color_viridis(option = "A")
##
saveRDS(maize_sub_ucell,"maize_sub_ucell.rds")
