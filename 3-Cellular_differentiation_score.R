library(Seurat)
library(ggplot2)
library(UCell)
library(viridis)
##load 6dap and 7dap marker gene
maize_sub <- readRDS("maize_sub.rds")
dap6markergene<-read.csv("day6marker.csv")
dap7markergene<-read.csv("day7marker.csv")
##Marker genes take intersection with genes detected in single cells
dap6markergene<-intersect(rownames(submaize),dap6markergene$day6)
dap8markergene<-intersect(rownames(submaize),dap7markergene$day7)
# Gene sets that need input to build UCell
markers <- list()
markers$dap6<-dap6markergene
markers$dap8<-dap8markergene
###
submaize_ucell <- AddModuleScore_UCell(submaize, features = markers,maxRank = nrow(submaize))
cell_information <- submaize_ucell@meta.data
###
cell_information$differentiation_score<-cell_information$dap8_UCell/cell_information$dap6_UCell
cell_information$differentiation_score<- log2((1+cell_information[,22]))
submaize_ucell@meta.data<- cell_information
##
FeaturePlot(submaize_ucell, reduction = "umap", features = "differentiation_score", 
            ncol = 1, order = T,
            min.cutoff = "q03", max.cutoff = "q99", pt.size = 0.1)+scale_color_viridis(option = "A")
##
saveRDS(maize_sub_ucell,"submaize_ucell.rds")
