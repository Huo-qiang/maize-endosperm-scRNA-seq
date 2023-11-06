################Calculating the AUCell value for regulon#########################
#load GRN#
library(AUCell)
library（Seurat)
maize_sub <- readRDS("maize_sub.rds)
grn <- read.table("resultdap.txt")
##########################################regulon#####################
geneSets <- lapply(unique(grn$source),function(x){grn$target[grn$source==x]})
names(geneSets) <- unique(grn$source)
# Calculate enrichment scores
exprMatrix <- maize_sub@assays$RNA@counts
cells_AUC_0.01 <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.01)
cells_AUC_0.05 <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC_0.1 <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.1)
cells_AUC_0.15 <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.15)
cells_AUC_0.2 <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.2)
# Optional: Set the assignment thresholds
par(mfrow=c(3,3))
set.seed(123)
cells_assignment_0.01 <- AUCell_exploreThresholds(cells_AUC_0.01, plotHist=TRUE, nCores=64, assign=TRUE)
cells_assignment_0.05 <- AUCell_exploreThresholds(cells_AUC_0.05, plotHist=TRUE, nCores=64, assign=TRUE)
cells_assignment_0.1 <- AUCell_exploreThresholds(cells_AUC_0.1, plotHist=TRUE, nCores=64, assign=TRUE)
cells_assignment_0.15 <- AUCell_exploreThresholds(cells_AUC_0.15, plotHist=TRUE, nCores=64, assign=TRUE)
cells_assignment_0.2 <- AUCell_exploreThresholds(cells_AUC_0.2, plotHist=TRUE, nCores=64, assign=TRUE)
########################################################################################################################
##Calculating RAS###
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)
#Load AUCell 
rasMat <- fread("AUCell.txt", sep = "\t", header = T, data.table = F) # data.frame
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
saveRDS(rasMat, "rasMat.rds")
#Load cell type infortation
cell.info <- maize_sub@meta.data
cell.types <- names(table(cell.info$CellType))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$CellType == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)
##Calculating RSS
rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
saveRDS(rssMat, "rssMat.rds")
rssMat <- readRDS("rssMat.rds")
binMat <- read.table("binary_mtx.txt", sep = "\t", header = T, row.names = 1, check.names = FALSE)
###Regulon Rank Plot
source("plotRegulonRank.R")
source("DimPlot.R")
all_Plot <- function(cell.type, regulon) {
  p.list <- list(
    PlotRegulonRank(rssMat, cell.type),
    DimPlot(cell.info, cell.type = cell.type),
    DimPlot(cell.info, regulon = regulon)
  )
  cowplot::plot_grid(plotlist = p.list, ncol = 3, 
                     rel_widths = c(3,5,5))
}
all_Plot("0", "WRKY95")
######################################################################################################################################
