############
rm(list=ls())
memory.limit(204800)

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(clustree)
library(monocle)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
## 扩大future.global
options(future.globals.maxSize= 10000000000)

##Setting 10x files location
data6day1_dir <- "F:/PhD/single cell/combined/6DAP-1"
data7day1_dir <- "F:/PhD/single cell/combined/7DAP-1"
data7day2_dir <- "F:/PhD/single cell/combined/7DAP-2"
data7day3_dir <- "F:/PhD/single cell/combined/7DAP-3"
#data8day1_dir <- "F:/PhD/single cell/combined/8DAP-1"
#data8day2_dir <- "F:/PhD/single cell/combined/8DAP-2"
#data8day3_dir <- "F:/PhD/single cell/combined/8DAP-3"

##load 10x files
data6day1.data <- Read10X(data.dir = data6day1_dir)
data7day1.data <- Read10X(data.dir = data7day1_dir)
data7day2.data <- Read10X(data.dir = data7day2_dir)
data7day3.data <- Read10X(data.dir = data7day3_dir)
#data8day1.data <- Read10X(data.dir = data8day1_dir)
#data8day2.data <- Read10X(data.dir = data8day2_dir)
#data8day3.data <- Read10X(data.dir = data8day3_dir)



###Create seurat objects and further filter the data, with each feature occurring in at least 3 cells

data6day1 <- CreateSeuratObject(counts = data6day1.data, project = "6day1", min.cells = 3)
data7day1 <- CreateSeuratObject(counts = data7day1.data, project = "7day1", min.cells = 3)
data7day2 <- CreateSeuratObject(counts = data7day2.data, project = "7day2", min.cells = 3)
data7day3 <- CreateSeuratObject(counts = data7day3.data, project = "7day3", min.cells = 3)
#data8day1 <- CreateSeuratObject(counts = data8day1.data, project = "8day1", min.cells = 3)
#data8day2 <- CreateSeuratObject(counts = data8day2.data, project = "8day2", min.cells = 3)
#data8day3 <- CreateSeuratObject(counts = data8day3.data, project = "8day3", min.cells = 3)

##Marker mitochondrial genes
data6day1 <- PercentageFeatureSet(data6day1, pattern = "^Zeam", col.name = "percent.mt")
data7day1 <- PercentageFeatureSet(data7day1, pattern = "^Zeam", col.name = "percent.mt")
data7day2 <- PercentageFeatureSet(data7day2, pattern = "^Zeam", col.name = "percent.mt")
data7day3 <- PercentageFeatureSet(data7day3, pattern = "^Zeam", col.name = "percent.mt")
#data8day1 <- PercentageFeatureSet(data8day1, pattern = "^Zeam", col.name = "percent.mt")
#data8day2 <- PercentageFeatureSet(data8day2, pattern = "^Zeam", col.name = "percent.mt")
#data8day3 <- PercentageFeatureSet(data8day3, pattern = "^Zeam", col.name = "percent.mt")

## Calculate the correlation between UMI and gene
data6day1$log10GenesPerUMI <- log10(data6day1$nFeature_RNA) / log10(data6day1$nCount_RNA)
data7day1$log10GenesPerUMI <- log10(data7day1$nFeature_RNA) / log10(data7day1$nCount_RNA)
data7day2$log10GenesPerUMI <- log10(data7day2$nFeature_RNA) / log10(data7day2$nCount_RNA)
data7day3$log10GenesPerUMI <- log10(data7day3$nFeature_RNA) / log10(data7day3$nCount_RNA)
#data8day1$log10GenesPerUMI <- log10(data8day1$nFeature_RNA) / log10(data8day1$nCount_RNA)
#data8day2$log10GenesPerUMI <- log10(data8day2$nFeature_RNA) / log10(data8day2$nCount_RNA)
#data8day3$log10GenesPerUMI <- log10(data8day3$nFeature_RNA) / log10(data8day3$nCount_RNA)

##View the mitochondrial gene ratio distribution map
# VlnPlot(data6day, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## look feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot7day1 <- FeatureScatter(data7day, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot7day2 <- FeatureScatter(data7day, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot7day2
# plot7day1 + plot7day2

## 查看数据质量
# metadata <- data7day@meta.data

# metadata %>% 
#  ggplot(aes(x=nCount_RNA)) + 
#  geom_density(alpha = 0.2) + 
#  scale_x_log10() + 
#  theme_classic() +
#  ylab("Cell density") +
#  geom_vline(xintercept = 500)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
# metadata %>%
#  ggplot(aes(x=log10GenesPerUMI)) +
#  geom_density(alpha = 0.2) +
#  theme_classic() +
#  geom_vline(xintercept = 0.80)

###Further filtering of cells: removal of data containing more than 2 cells feature< 8000; retention of mitochondrial genes with less than 5% of transcripts, removal of dead cells

data6day1 <- subset(data6day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day1 <- subset(data7day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day2 <- subset(data7day2, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day3 <- subset(data7day3, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
#data8day1 <- subset(data8day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
#data8day2 <- subset(data8day2, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
#data8day3 <- subset(data8day3, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )

## Adding grouping information to seurat objects
data6day1$day <- "6d"
data7day1$day <- "7d"
data7day2$day <- "7d"
data7day3$day <- "7d"
#data8day1$day <- "8d"
#data8day2$day <- "8d"
d#ata8day3$day <- "8d"
#cell cycle marker
scRNA_endosperm <- CellCycleScoring(scRNA_endosperm, 
                                    s.features = g1sgene, 
                                    g2m.features = g2mgene, 
                                    set.ident = TRUE)
##save RDS
save(data6day1,data7day1,data7day2,data7day3,data8day1,data8day2,data9day3,file="scRNAlist.Rdata")

##==Integration of multiple samples use harmony==##
library(harmony)
mazie_scRNAlist <- load("scRNAlist.Rdata")
scRNA_harmony <- merge(data6day1, y=data7day1, data7day2, data7day3)
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress =c("percent.mt")) %>% RunPCA(verbose=FALSE)
system.time({scRNA_endosperm <- RunHarmony(scRNA_endosperm, group.by.vars = "orig.ident")})
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:50) %>% FindClusters(resulotion=0.75)
scRNA_harmony <- CellCycleScoring(scRNA_harmony, s.features = g1sgene,g2m.features = g2mgene, set.ident = TRUE)
cell_cycle_before<-scRNA_harmony <- RunPCA(scRNA_harmony, features = c(g1sgene, g2mgene))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData(features = rownames(scRNA_harmony),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
cell_cycle_after <- RunPCA(scRNA_harmony, features = c(g1sgene, g2mgene))
#group_by_cluster
DimPlot(scRNA_harmony, reduction = "umap", label=T) 
#group_by_sample
DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#group_by_DAP
DimPlot(scRNA_harmony, reduction = "umap", group.by='day') 
#save files
saveRDS(scRNA_harmony, 'maize_endosperm_integrated.rds')

## annotated cell population
maize <- readRDS("maize_endosperm_integrated.rds")
AL1.marker.list <- c("pyruvate, orthophosphate dikinase2","sucrose synthase1","aleurone9","Zm00001d017852","Zm00001d008620","defensin-like protein2")
AL2.marker.list <- c("umc1664","Zm00001d033909","subtilisin1","Zm00001d024120","Zm00001d025924","histidine-containing phosphotransfer protein3")
AL3.marker.list <- c("Zm00001d043709","Zm00001d016379","Zm00001d007800","hexokinase5","Zm00001d023964","Zm00001d015569")
ALmarker.list <- c("pyruvate, orthophosphate dikinase2","sucrose synthase1","aleurone9","Zm00001d017852","Zm00001d008620","defensin-like protein2","umc1664","Zm00001d033909","subtilisin1","Zm00001d024120","Zm00001d025924","histidine-containing phosphotransfer protein3","Zm00001d043709","Zm00001d016379","Zm00001d007800","hexokinase5","Zm00001d023964","Zm00001d015569")
BETL.marker.list <- c("TF-Myb related protein1","basal layer antifungal protein2","basal endosperm transfer layer9","Zm00001d053785","Zm00001d052759","Zm00001d053108")
CSE.marker.list <- c("TF-floury3","Zm00001d017286","cytokinin N-glucosyl transferase1","Zm00001d010434","Zm00001d038682","Zm00001d009415")
ESR.marker.list <- c("maternally expressed gene14","Zm00001d019032","embryo surrounding region2","Zm00001d026755","embryo surrounding region1","Zm00001d038758")
EMB.marker.list <- c("Zm00001d011342","Zm00001d011340","Zm00001d011345","Zm00001d022089","TF-Zea mays MADS6","sugars will eventually be exported transporter4a")

# How many cells are in each cluster
table(Idents(maize))

# How many cells are in each day?
table(maize$day)
table(maize$orig.ident)

# Frequency statistics
prop.table(table(Idents(maize), maize$day), margin = 2)

write.csv((prop.table(table(Idents(maize), maize$day), margin = 2)), "cell-propotions-by-day-res0.75.csv")

# library(RColorBrewer)
# display.brewer.all()
# colorset <- brewer.pal(10,"RdYlBu")
# colorset <- brewer.pal(9,"YlOrRd")


pdf("Dotplot-res.0.75-AL-markers-assay-RNA.pdf",paper = "letter")
AL=DotPlot(maize, features = AL.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.0.75-BETL-markers-assay-RNA.pdf",paper = "letter")
BETL=DotPlot(maize, features = BETL.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.0.75-ESR-markers-assay-RNA.pdf",paper = "letter")
ESR=DotPlot(maize, features = ESR.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.0.75-CSE-markers-assay-RNA.pdf",paper = "letter")
CSE=DotPlot(maize, features = CSE.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.0.75-EMB-markers-assay-RNA.pdf",paper = "letter")
EMB=DotPlot(maize, features = EMB.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

## find all markers
maize.res0.75.RNA.markers.wilcox <- FindAllMarkers(maize, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
maize.res0.75.RNA.markers.wilcox %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% print(n=Inf)
write.csv(mmaize.res0.75.RNA.markers.wilcox, "all.res0.75.RNA.markers.wilcox.csv")
## find day marker
Idents(maize) <- "day"
maize.DAP.RNA.markers.wilcox <- FindAllMarkers(maize, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
maize.res0.75.markers.wilcox %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% print(n=Inf)
write.csv(maize.DAP.RNA.markers.wilcox, "all.DAP.RNA.markers.wilcox.csv")
##Average expression level
Idents(maize) <- "seurat_clusters"
clusters.res.0.75.averagexpression.integreted.RNA.data <- AverageExpression(maize, assays = "RNA", slot = "data")
write.csv(as.matrix(clusters.res.0.75.averagexpression.integreted.RNA.data$RNA), "RNA-average.expression-res.0.75.csv")
clusters.res.0.75.averagexpression.split.RNA.data <- AverageExpression(maize, assays = "RNA", slot = "data", add.ident = ("day"))
write.csv(as.matrix(clusters.res.0.75.averagexpression.split.RNA.data$RNA), "RNA-average.expression-splitbyday-res.0.75.csv"
## Annotation of cell populations based on annotation and in situ hybridisation
maize$celltype <- NA
maize$celltype[WhichCells(object = maize, idents = c(0))] <- "aleurone Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(4))] <- "aleurone Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(5))] <- "aleurone Ⅲ"
maize$celltype[WhichCells(object = maize, idents = c(11))] <- "aleurone Ⅳ"
maize$celltype[WhichCells(object = maize, idents = c(16))] <- "endosperm adjacent to scutellum"
maize$celltype[WhichCells(object = maize, idents = c(18))] <- "embryo surrounding region"
maize$celltype[WhichCells(object = maize, idents = c(1))] <- "starch endosperm Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(3))] <- "starch endosperm Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(7))] <-"basal endosperm transfer layer Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(13))] <- "basal endosperm transfer layer Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(2))] <- "Synthesis phase Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(6))] <- "Synthesis phase Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(9))] <- "Synthesis phase Ⅲ"
maize$celltype[WhichCells(object = maize, idents = c(12))] <- "Mitosis phase Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(14))] <- "Mitosis phase Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(15))] <- "Mitosis phase Ⅲ"
maize$celltype[WhichCells(object = maize, idents = c(8))] <- "undefined Ⅰ"
maize$celltype[WhichCells(object = maize, idents = c(10))] <- "undefined Ⅱ"
maize$celltype[WhichCells(object = maize, idents = c(17))] <-"stimulus response"
maize$celltype[WhichCells(object = maize, idents = c(19))] <- "nucellus"
maize$celltype[WhichCells(object = maize, idents = c(20))] <- "embryo"
maize$celltype[WhichCells(object = maize, idents = c(21))] <- "mitochondrial RNA enriched"
DimPlot(maize, reduction = "umap", label = TRUE,group.by ="celltype")
DimPlot(maize, reduction = "tsne", label = TRUE,group.by ="celltype")
## subset data
maize_sub<- subset(x = maize, subset = seurat_cluster == c("0","3","4","7","6","1","8","5","11","15","10","12"))
## re-normalize
maize_sub <- NormalizeData(maize_sub) %>% ScaleData(features = rownames(maize_sub),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
maize_sub <- RunUMAP(maize_sub, reduction = "harmony", dims = 1:50)
maize_sub <- FindNeighbors(maize_sub, reduction = "harmony", dims = 1:50)
DimPlot(maize_sub,reduction = "umap",label = T)
## subset data
maize_used<- subset(x = maize, subset = seurat_cluster == c("0","3","4","7","6","1","8","5","11","15"))
maize_used <- NormalizeData(maize_used) %>% ScaleData(features = rownames(maize_used),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
maize_used <- RunUMAP(maize_used, reduction = "harmony", dims = 1:50)
maize_used <- FindNeighbors(maize_used, reduction = "harmony", dims = 1:50)
DimPlot(maize_used,reduction = "umap",label = T)
#
## choose the resolution for downstream analysis
Idents(maize_used) <- "RNA_snn_res.0.75"
levels(maize_used)
## subset clusters
clusters.SE <- subset(x = maize_used, idents = c("1","8"), invert = F)
clusters.SE <- NormalizeData(clusters.SE) %>% ScaleData(features = rownames(clusters.SE),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
clusters.SE.markers.wilcox <- FindAllMarkers(clusters.SE, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
write.csv(clusters.SE.markers.wilcox, "2023-10-20-clusters.SE.markers.wilcox.csv")

clusters.BETL <- subset(x = maize_used, idents = c("2","9"), invert = F)
clusters.BETL <- NormalizeData(clusters.BETL) %>% FindVariableFeatures() %>% ScaleData(features = rownames(clusters.BETL),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
clusters.BETL.markers.wilcox <- FindAllMarkers(clusters.BETL, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
write.csv(clusters.BETL.markers.wilcox, "2023-10-20-clusters.BETL.markers.wilcox.csv")

clusters.AL <- subset(x = maize_used, idents = c("0","3","4","7"), invert = F)
clusters.AL <- NormalizeData(clusters.AL) %>% FindVariableFeatures() %>% ScaleData(features = rownames(clusters.AL),vars.to.regress =c("percent.mt","S.Score","G2M.Score")) %>% RunPCA(verbose=FALSE)
clusters.AL.markers.wilcox <- FindAllMarkers(clusters.AL, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
write.csv(clusters.AL.markers.wilcox, "2023-10-20-clusters.AL.markers.wilcox.csv")


## choose markers
SE.sub.markers.p005.strict <-subset(clusters.SE.markers.wilcox, p_val_adj < 0.05 & pct.1 > pct.2) 

top.SE.sub.markers <- SE.sub.markers.p005.strict %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top.SE.sub.markers <- top.SE.sub.markers$gene
top.SE.sub.markers.deduplicate <- top.SE.sub.markers[!duplicated(top.SE.sub.markers)]
write.csv(top.SE.sub.markers.deduplicate, "SE.sub.marker.label.csv")

##
BETL.sub.markers.p005.strict <-subset(clusters.BETL.markers.wilcox, p_val_adj < 0.05 & pct.1 > pct.2) 

top.BETL.sub.markers <- BETL.sub.markers.p005.strict %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top.BETL.sub.markers <- top.BETL.sub.markers$gene
top.BETL.sub.markers.deduplicate <- top.BETL.sub.markers[!duplicated(top.BETL.sub.markers)]
write.csv(top.BETL.sub.markers.deduplicate, "BETL.sub.marker.label.csv")

##
AL.sub.markers.p005.strict <-subset(clusters.AL.markers.wilcox, p_val_adj < 0.05 & pct.1 > pct.2) 

top.AL.sub.markers <- AL.sub.markers.p005.strict %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top.AL.sub.markers <- top.AL.sub.markers$gene
top.AL.sub.markers.deduplicate <- top.AL.sub.markers[!duplicated(top.AL.sub.markers)]
write.csv(top.AL.sub.markers.deduplicate, "AL.sub.marker.label.csv")

##############################################################################################################

## Plot
DotPlot(clusters.SE, features = top.SE.sub.markers.deduplicate, assay = "RNA", cluster.idents = T, group.by = "RNA_snn_res.0.75", cols = c("#FFFFCC","#A50026")) + theme(axis.text.x = element_text(angle = 90,hjust=1)) + FontSize(2)

DotPlot(clusters.BETL, features = top.BETL.sub.markers.deduplicate, assay = "RNA", cluster.idents = T, group.by = "RNA_snn_res.0.75", cols = c("#FFFFCC","#A50026")) + theme(axis.text.x = element_text(angle = 90,hjust=1)) + FontSize(5)

DotPlot(clusters.AL, features = top.AL.sub.markers.deduplicate, assay = "RNA", cluster.idents = T, group.by = "RNA_snn_res.0.75", cols = c("#FFFFCC","#A50026")) + theme(axis.text.x = element_text(angle = 90,hjust=1)) + FontSize(2)


######################################################################################

###find day marker
bulk678.expression.rpm <- read.csv("D:/0-SingleCell/5-real-time-course-results/bulk678expression.rpm.csv",header = TRUE,row.names="full.name")

cluster_marker <- read.csv("D:/0-SingleCell/2023_10_18_cluster_marker.csv")

name.to.symbol <- read.csv("D:/0-SingleCell/fullname.to.symbol.csv",header = TRUE,row.names="full.name")

name.to.symbol$gene <- row.names(name.to.symbol)

## choose the resolution for downstream analysis
##
Idents(maize_day7) <- "orig.ident"
table(Idents(maize_day7))
data.cluster.sample.select <- subset(x = maize_day7, idents = c("6dap-1"), invert = F)

## set loop
# test.cluster.list <- c(0,3,4,7,1,8,2,9,11,15,5,6)
test.cluster.list <- c(15)

for (y in test.cluster.list) {
  
  id <- y
  
  cluster.id <- paste("cluster_C",id,sep="")
  
  ## extract cluster marker list
  c.marker <- cluster_marker[cluster_marker$If_cluster_marker == id,]
  
  ##
  Idents(data.cluster.sample.select) <- "RNA_snn_res.0.75"
  data.cluster <- subset(x = data.cluster.sample.select, idents = c(id), invert = F)
  
  ## re-normalize
  data.cluster <- NormalizeData(data.cluster) %>% ScaleData(vars.to.regress =c("percent.mt","S.Score","G2M.Score"))
  ## find DEG
  Idents(data.cluster) <- "day"
  data.cluster.markers.day.wilcox <- FindAllMarkers(data.cluster, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25, exact = F)
  
  data.cluster.markers.day.filtered <-subset(data.cluster.markers.day.wilcox, p_val_adj < 0.05 & pct.1 > pct.2)  # avg_log2FC > 0.6)
  
  
  ## data.cluster.markers.day.filtered <- data.cluster.markers.day.filtered %>% group_by(gene) %>% filter(n()==1)
  
  ## row.names(data.cluster.markers.day.filtered) <- data.cluster.markers.day.filtered$gene
  
  
  averagexpression.cluster.RNA.data <- AverageExpression(data.cluster, assays = "RNA", slot = "data", group.by = "day")
  
  cluster.mat <- as.matrix(averagexpression.cluster.RNA.data$RNA)
  
  cluster.day.filtered <- data.cluster.markers.day.filtered$gene
  
  cluster.day.filtered.deduplicate <- cluster.day.filtered[!duplicated(cluster.day.filtered)]
  
  cluster.mat.filtered.deduplicate <- as.matrix(cluster.mat[cluster.day.filtered.deduplicate,])
  
  ##scaledata(0.1)
  mat.filtered.scaled <- t(apply(cluster.mat.filtered.deduplicate, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  ## plot heatmap with row label
  pdf(paste(cluster.id,"Heatmap_byday_rowlabel.pdf",sep=""), # File name
      width = 6, height = 8, # Width and height in inches
  )
  draw(Heatmap(mat.filtered.scaled,
               column_title = cluster.id,
               col=colorRamp2(c(0,0.5,1),c("#728EB5","#F5EDAE","#CA3B32")),
               cluster_rows = T,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               cluster_columns = F,
               cluster_column_slices = F,
               show_column_names = T,
               show_row_names = T,
               row_names_gp = gpar(fontsize = 1),
               heatmap_width = unit(4, "cm"), heatmap_height = unit(15, "cm"), row_dend_width = unit(0.2, "cm"),
               heatmap_legend_param = list(title = "NA", legend_height = unit(3, "cm"))
  ), 
  
  heatmap_legend_side = "left")
  dev.off()
  
  ### heatmap with highlight
  
  mat.order <- data.frame(list=1:length(cluster.day.filtered.deduplicate))
  mat.order.list <- cbind(cluster.mat.filtered.deduplicate, mat.order)
  mat.order.list$gene <- row.names(mat.order.list)  
  merge.mat.list <- merge(mat.order.list, c.marker, by= c("gene"), all.x=TRUE)
  
  merge.mat.list.label <- merge(merge.mat.list, name.to.symbol, by= c("gene"), all.x=TRUE)
  
  write.csv(merge.mat.list.label, paste(cluster.id,"_Heatmap_byday_label.csv",sep="") )
  
  mark_gene_all <- merge.mat.list.label[!is.na(merge.mat.list.label$If_cluster_marker), ]
  mark_gene <- mark_gene_all[ , c("Symbol")]
  gene_pos <- mark_gene_all[ , c("list")]
  
  row_anno <- rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene,labels_gp = gpar(fontsize = 0.5)) )
  
  pdf(paste(cluster.id,"Heatmap_byday_markers.pdf",sep=""), # File name
      width = 6, height = 8, # Width and height in inches
  )
  draw(Heatmap(mat.filtered.scaled,
               column_title = cluster.id,
               col=colorRamp2(c(0,0.5,1),c("#728EB5","#F5EDAE","#CA3B32")),
               cluster_rows = T,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               cluster_columns = F,
               cluster_column_slices = F,
               show_column_names = T,
               show_row_names = F,
               right_annotation = row_anno,
               row_names_gp = gpar(fontsize = 1),
               heatmap_width = unit(4, "cm"), heatmap_height = unit(15, "cm"), row_dend_width = unit(0.2, "cm"),
               heatmap_legend_param = list(title = "NA", legend_height = unit(3, "cm"))
  ), 
  
  heatmap_legend_side = "left")
  dev.off()
  
  ##
  mark_gene_top <- mark_gene_all[mark_gene_all$pct_fold_cluster > 2, ]
  mark_gene_top <- mark_gene_top %>% top_n(n = 5, wt = avg_log2FC_cluster)
  mark_gene_2 <- mark_gene_top[ , c("Symbol")]
  gene_pos_2 <- mark_gene_top[ , c("list")]
  
  row_anno_2 <- rowAnnotation(mark_gene_2 = anno_mark(at = gene_pos_2,labels = mark_gene_2,labels_gp = gpar(fontsize = 30)) )
  
  pdf(paste(cluster.id,"Heatmap_byday_highlight.pdf",sep=""), # File name
      width = 6, height = 8, # Width and height in inches
  )
  draw(Heatmap(mat.filtered.scaled,
               column_title = cluster.id,
               col=colorRamp2(c(0,0.5,1),c("#728EB5","#F5EDAE","#CA3B32")),
               cluster_rows = T,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               cluster_columns = F,
               cluster_column_slices = F,
               show_column_names = T,
               show_row_names = F,
               right_annotation = row_anno_2,
               row_names_gp = gpar(fontsize = 0.5),
               heatmap_width = unit(8, "cm"), heatmap_height = unit(15, "cm"), row_dend_width = unit(0.2, "cm"),
               heatmap_legend_param = list(title = "NA", legend_height = unit(3, "cm"))
  ), 
  
  heatmap_legend_side = "left")
  dev.off()
  
  ## top5 by day 
  top_by_day <- data.cluster.markers.day.filtered %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  merge.top_by_day <- merge(mat.order.list, top_by_day, by= c("gene"), all.x=TRUE)
  
  merge.top_by_day.label <- merge(merge.top_by_day, name.to.symbol, by= c("gene"), all.x=TRUE)
  
  mark_gene_top_day  <- merge.top_by_day.label[!is.na(merge.top_by_day.label$cluster), ]
  
  mark_gene_3 <- mark_gene_top_day[ , c("Symbol")]
  gene_pos_3 <- mark_gene_top_day[ , c("list")]
  
  row_anno_3 <- rowAnnotation(mark_gene_3 = anno_mark(at = gene_pos_3,labels = mark_gene_3,labels_gp = gpar(fontsize = 20)) )
  
  pdf(paste("topbyday",cluster.id,"Heatmap_byday_highlight.pdf",sep=""), # File name
      width = 6, height = 8, # Width and height in inches
  )
  draw(Heatmap(mat.filtered.scaled,
               column_title = cluster.id,
               col=colorRamp2(c(0,0.5,1),c("#728EB5","#F5EDAE","#CA3B32")),
               cluster_rows = T,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               cluster_columns = F,
               cluster_column_slices = F,
               show_column_names = T,
               show_row_names = F,
               right_annotation = row_anno_3,
               row_names_gp = gpar(fontsize = 0.5),
               heatmap_width = unit(8, "cm"), heatmap_height = unit(15, "cm"), row_dend_width = unit(0.2, "cm"),
               heatmap_legend_param = list(title = "NA", legend_height = unit(3, "cm"))
  ), 
  
  heatmap_legend_side = "left")
  dev.off()
  
  
  ## export data
  filterd.marker.logFC <- data.cluster.markers.day.filtered[data.cluster.markers.day.filtered$gene %in% cluster.day.filtered.deduplicate, ]
  
  filterd.expression <- as.data.frame(cluster.mat)
  filterd.expression$gene <- row.names(filterd.expression)                     # Apply row.names function
  
  merge.c <- merge(filterd.marker.logFC, filterd.expression, by= c("gene"), all.x=TRUE)
  
  bulk678.expression.rpm$gene <- row.names(bulk678.expression.rpm)                     # Apply row.names function
  
  merge.bulk <- merge(merge.c, bulk678.expression.rpm, by= c("gene"), all.x=TRUE)
  
  merge.bulk.ifmarker <- merge(merge.bulk, c.marker, by= c("gene"), all.x=TRUE)
  
  write.csv(merge.bulk.ifmarker,paste("2023-10-22-",cluster.id,".DEGbyday.expression.csv",sep=""))
  
}
    
