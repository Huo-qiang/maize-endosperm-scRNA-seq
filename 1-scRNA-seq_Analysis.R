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

##设置文件存储位置
data6day1_dir <- "F:/PhD/single cell/combined/6DAP-1"
data7day1_dir <- "F:/PhD/single cell/combined/7DAP-1"
data7day2_dir <- "F:/PhD/single cell/combined/7DAP-2"
data7day3_dir <- "F:/PhD/single cell/combined/7DAP-3"
data8day1_dir <- "F:/PhD/single cell/combined/8DAP-1"
data8day2_dir <- "F:/PhD/single cell/combined/8DAP-2"
data8day3_dir <- "F:/PhD/single cell/combined/8DAP-3"

##加载10X数据-这里用clean数据
data6day1.data <- Read10X(data.dir = data6day1_dir)
data7day1.data <- Read10X(data.dir = data7day1_dir)
data7day2.data <- Read10X(data.dir = data7day2_dir)
data7day3.data <- Read10X(data.dir = data7day3_dir)
data8day1.data <- Read10X(data.dir = data8day1_dir)
data8day2.data <- Read10X(data.dir = data8day2_dir)
data8day3.data <- Read10X(data.dir = data8day3_dir)



###创建seurat对象并进一步过滤数据，每个feature至少在3个细胞里出现

data6day1 <- CreateSeuratObject(counts = data6day1.data, project = "6day1", min.cells = 3)
data7day1 <- CreateSeuratObject(counts = data7day1.data, project = "7day1", min.cells = 3)
data7day2 <- CreateSeuratObject(counts = data7day2.data, project = "7day2", min.cells = 3)
data7day3 <- CreateSeuratObject(counts = data7day3.data, project = "7day3", min.cells = 3)
data8day1 <- CreateSeuratObject(counts = data8day1.data, project = "8day1", min.cells = 3)
data8day2 <- CreateSeuratObject(counts = data8day2.data, project = "8day2", min.cells = 3)
data8day3 <- CreateSeuratObject(counts = data8day3.data, project = "8day3", min.cells = 3)

##标记线粒体基因
data6day1 <- PercentageFeatureSet(data6day1, pattern = "^Zeam", col.name = "percent.mt")
data7day1 <- PercentageFeatureSet(data7day1, pattern = "^Zeam", col.name = "percent.mt")
data7day2 <- PercentageFeatureSet(data7day2, pattern = "^Zeam", col.name = "percent.mt")
data7day3 <- PercentageFeatureSet(data7day3, pattern = "^Zeam", col.name = "percent.mt")
data8day1 <- PercentageFeatureSet(data8day1, pattern = "^Zeam", col.name = "percent.mt")
data8day2 <- PercentageFeatureSet(data8day2, pattern = "^Zeam", col.name = "percent.mt")
data8day3 <- PercentageFeatureSet(data8day3, pattern = "^Zeam", col.name = "percent.mt")

## 计算UMI和gene的相关性
data6day1$log10GenesPerUMI <- log10(data6day1$nFeature_RNA) / log10(data6day1$nCount_RNA)
data7day1$log10GenesPerUMI <- log10(data7day1$nFeature_RNA) / log10(data7day1$nCount_RNA)
data7day2$log10GenesPerUMI <- log10(data7day2$nFeature_RNA) / log10(data7day2$nCount_RNA)
data7day3$log10GenesPerUMI <- log10(data7day3$nFeature_RNA) / log10(data7day3$nCount_RNA)
data8day1$log10GenesPerUMI <- log10(data8day1$nFeature_RNA) / log10(data8day1$nCount_RNA)
data8day2$log10GenesPerUMI <- log10(data8day2$nFeature_RNA) / log10(data8day2$nCount_RNA)
data8day3$log10GenesPerUMI <- log10(data8day3$nFeature_RNA) / log10(data8day3$nCount_RNA)

##查看线粒体基因比例分布图
# VlnPlot(data6day, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## 查看 feature-feature relationships, but can be used
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

###进一步过滤细胞：去掉包含2个以上细胞的数据 feature< 8000；保留线粒体基因的转录本数低于5%的细胞，去掉死细胞

data6day1 <- subset(data6day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day1 <- subset(data7day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day2 <- subset(data7day2, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data7day3 <- subset(data7day3, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data8day1 <- subset(data8day1, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data8day2 <- subset(data8day2, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )
data8day3 <- subset(data8day3, subset = nFeature_RNA >500 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80 )

## 为seurat对象添加分组信息
data6day1$day <- "6d"
data7day1$day <- "7d"
data7day2$day <- "7d"
data7day3$day <- "7d"
data8day1$day <- "8d"
data8day2$day <- "8d"
data8day3$day <- "8d"
##存储上游分析文件
save(data6day1,data7day1,data7day2,data7day3,data8day1,data8day2,data9day3,file="scRNAlist.Rdata")

##==harmony整合多样本==##
library(harmony)
mazie_scRNAlist <- load("scRNAlist.Rdata")
##PCA降维
scRNA_endosperm <- merge(data6day1, y=data7day1, data7day2, data7day3, data8day1, 
                       data8day2, data9day3)
scRNA_endosperm <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
system.time({scRNA_endosperm <- RunHarmony(scRNA_endosperm, group.by.vars = "orig.ident")})
#降维聚类
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:50) %>% FindClusters(resulotion=1.2)
##作图
#group_by_cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
#group_by_sample
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#group_by_DAP
plot3 = DimPlot(scRNA_harmony, reduction = "umap", group.by='day') 
#combinate
plot_all <- plot1+plot2+plot3
ggsave("scRNA_endosperm.png", plot = plot_all, width = 10, height = 5)
saveRDS(scRNA_endosperm, 'maize_endosperm_integrated.rds')

## 存储marker gene 列表
maize <- readRDS("20201229dataintegrated.rds")
AL1.marker.list <- c("pyruvate, orthophosphate dikinase2","sucrose synthase1","aleurone9","Zm00001d017852","Zm00001d008620","defensin-like protein2")
AL2.marker.list <- c("umc1664","Zm00001d033909","subtilisin1","Zm00001d024120","Zm00001d025924","histidine-containing phosphotransfer protein3")
AL3.marker.list <- c("Zm00001d043709","Zm00001d016379","Zm00001d007800","hexokinase5","Zm00001d023964","Zm00001d015569")
ALmarker.list <- c("pyruvate, orthophosphate dikinase2","sucrose synthase1","aleurone9","Zm00001d017852","Zm00001d008620","defensin-like protein2","umc1664","Zm00001d033909","subtilisin1","Zm00001d024120","Zm00001d025924","histidine-containing phosphotransfer protein3","Zm00001d043709","Zm00001d016379","Zm00001d007800","hexokinase5","Zm00001d023964","Zm00001d015569")
BETL.marker.list <- c("TF-Myb related protein1","basal layer antifungal protein2","basal endosperm transfer layer9","Zm00001d053785","Zm00001d052759","Zm00001d053108")
CSE.marker.list <- c("TF-floury3","Zm00001d017286","cytokinin N-glucosyl transferase1","Zm00001d010434","Zm00001d038682","Zm00001d009415")
ESR.marker.list <- c("maternally expressed gene14","Zm00001d019032","embryo surrounding region2","Zm00001d026755","embryo surrounding region1","Zm00001d038758")
EMB.marker.list <- c("Zm00001d011342","Zm00001d011340","Zm00001d011345","Zm00001d022089","TF-Zea mays MADS6","sugars will eventually be exported transporter4a")
## 统计分群信息

# How many cells are in each cluster
table(Idents(maize))

# How many cells are in each day?
table(maize$day)
table(maize$orig.ident)

# 频率统计
prop.table(table(Idents(maize), maize$day), margin = 2)

write.csv((prop.table(table(Idents(maize), maize$day), margin = 2)), "cell-propotions-by-day-res1.2.csv")

## heatmap选颜色
# library(RColorBrewer)
# display.brewer.all()
# colorset <- brewer.pal(10,"RdYlBu")
# colorset <- brewer.pal(9,"YlOrRd")


## 先用Dotplot和featureplot看一下细胞类型

##DotPlot(maize, features = AL.marker.list, cols = c(至少三种颜色), split.by","= "day") + RotatedAxis()

pdf("Dotplot-res.1.2-AL-markers-assay-RNA.pdf",paper = "letter")
AL=DotPlot(maize, features = AL.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.1.2-BETL-markers-assay-RNA.pdf",paper = "letter")
BETL=DotPlot(maize, features = BETL.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.1.2-ESR-markers-assay-RNA.pdf",paper = "letter")
ESR=DotPlot(maize, features = ESR.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.1.2-CSE-markers-assay-RNA.pdf",paper = "letter")
CSE=DotPlot(maize, features = CSE.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()

pdf("Dotplot-res.1.2-EMB-markers-assay-RNA.pdf",paper = "letter")
EMB=DotPlot(maize, features = EMB.marker.list, cols = c("#FFFFCC","#A50026"), assay = "RNA") + RotatedAxis()
dev.off()


## cols = (DiscretePalette(22,"stepped")) 调整颜色参数

## 画图看看几个feature的分布情况
FeaturePlot(maize, features = c("nFeature_RNA", "nFeature_SCT","log10GenesPerUMI", "percent.mt"), min.cutoff = "q10", max.cutoff = "q90",label = TRUE, ncol = 2)
## 用subset函数选取子集
data.6day.only.all <- subset(x = maize, subset = day == c("6d"))
data.7day.only.all <- subset(x = maize, subset = day == c("7d"))
data.8day.only.all <- subset(x = maize, subset = day == c("8d"))
## 分群marker基因的鉴定
maize.res1.2.RNA.markers.wilcox <- FindAllMarkers(maize, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
maize.res1.2.markers.wilcox %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% print(n=Inf)
write.csv(maize.res1.2.RNA.markers.wilcox, "all.res1.2.RNA.markers.wilcox.csv")
## DAP的marker基因鉴定
Idents(maize) <- "day"
maize.DAP.RNA.markers.wilcox <- FindAllMarkers(maize, assay = "RNA", slot = "data", test.use = "wilcox" , only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
maize.res1.2.markers.wilcox %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% print(n=Inf)
write.csv(maize.DAP.RNA.markers.wilcox, "all.DAP.RNA.markers.wilcox.csv")
# 频率统计
prop.table(table(Idents(maize), maize$day), margin = 2)
write.csv((prop.table(table(Idents(maize), maize$day), margin = 2)), "20210114-cell-propotions-by-day-res1.2.csv")
##平均表达量统计
Idents(maize) <- "seurat_clusters"
clusters.res.1.2.averagexpression.integreted.RNA.data <- AverageExpression(maize, assays = "RNA", slot = "data")
write.csv(as.matrix(clusters.res.1.2.averagexpression.integreted.RNA.data$RNA), "RNA-average.expression-res.1.2.csv")
clusters.res.1.2.averagexpression.split.RNA.data <- AverageExpression(maize, assays = "RNA", slot = "data", add.ident = ("day"))
write.csv(as.matrix(clusters.res.1.2.averagexpression.split.RNA.data$RNA), "RNA-average.expression-splitbyday-res.1.2.csv"
## 根据注释和原位杂交注释细胞群
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
## 用subset函数选取子集
maize_sub<- subset(x = maize, subset = seurat_cluster == c("0","1","3","4","5","7","11","13","16","18"))
data.6day.only.sub <- subset(x = maize_sub, subset = day == c("6d"))
data.7day.only.sub <- subset(x = maize_sub, subset = day == c("7d"))
data.8day.only.sub <- subset(x = maize_sub, subset = day == c("8d"))
save(maize_sub,data.6day.only.sub,data.7day.only.sub,data.8day.only.sub,file="maize_sub.Rdata")
##重新plot
load("maize_sub.Rdata")
maize_sub <- FindNeighbors(maize_sub, dims = 1:50)
maize_sub <- RunUMAP(maize_sub, dims = 1:50)
DimPlot(maize_sub, group.by = "seurat_cluster", reduction = "umap", label = TRUE, split.by = "day", label.size = 10, repel = TRUE, pt.size = 0.5, ncol = 3)

## 统计分群信息
# How many cells are in each cluster
table(Idents(maize_sub))
# How many cells are in each day?
table(maize_sub$day)
table(maize_sub$orig.ident)

# 频率统计
prop.table(table(Idents(maize_sub), maize_sub$day), margin = 2)
write.csv((prop.table(table(Idents(maize_sub), maize_sub$day), margin = 2)), "maize_sub_cell-propotions-by-day-res1.2.csv")
#
saveRDS(maize_sub,"maize_sub.rds")


#############################################################################################################################################
