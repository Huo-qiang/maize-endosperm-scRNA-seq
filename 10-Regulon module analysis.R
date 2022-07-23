###Load packages
library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)
###Load data
rasMat <- readRDS("rasMat.rds")
regulons <- read.table("regulons.txt", sep = "\t")
regulon.names <- sapply(strsplit(regulons$V1, split = "\\("), function(x) x[1])
regulon.sizes <- sapply(strsplit(regulons$V3, split = ","), function(x) length(x))
regulon.names <- regulon.names[regulon.sizes>=5]
rasMat <- rasMat[, regulon.names]
pccMat <- cor(rasMat)
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}
csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)
csiMat.binary <- matrix(as.numeric(csiMat >= 0.6), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
saveRDS(csiMat, "csiMat.rds")
write.table(csiMat.binary, "csiMat.binary.txt", sep = "\t")
####module
mat = readRDS("csiMat.rds")
h = 4.5
row_dend = as.dendrogram(hclust(dist(mat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)
##Heatmap
col_range = c(0.6, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

ht <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

lgd <- Legend(
  col_fun = col_fun, 
  title = "", 
  at = col_range, 
  labels = c("low", "high"), 
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
  )
  draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
decorate_heatmap_body("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    x1 = x1/length(ind)
    x2 = x2/length(ind)
    grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2), 
              hjust = 0, vjust = 0, default.units = "npc", 
              gp = gpar(fill=NA, col="#FCB800", lwd=3))
    grid.text(label = paste0("M",clusters),
              x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
              default.units = "npc",
              hjust = 1, vjust = 0.5,
              gp = gpar(fontsize=12, fontface="bold"))
})
decorate_column_dend("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
        default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
})
##umap
tree = column_dend(ht)
ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))

write.table(regulon.clusters, "regulon_clusters.txt", sep = "\t", quote = F, row.names = F)
k = length(clusters)
cell.info <- maize_sub@meta.data
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])
p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(UMAP_1, UMAP_2, color=get(module))) + 
    geom_point(size=0.05) + 
    theme_bw(base_size = 15) + 
    ggtitle(module) + 
    scale_color_gradientn(name = NULL, colors = expression.color) + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
          )
})
cowplot::plot_grid(plotlist = p.list, ncol = 2)
##########################################################################################################################################
