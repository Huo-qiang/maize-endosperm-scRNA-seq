## Empty history space
ls()
rm(list=ls())
memory.limit(10000000)

## load packages
library(GenomicFeatures)
library(ChIPseeker)
library(clusterProfiler)

options(ChIPseeker.downstreamDistance = 500)

## options(ChIPseeker.ignore_1st_exon = TRUE)
## options(ChIPseeker.ignore_1st_intron = TRUE)
## options(ChIPseeker.ignore_downstream = TRUE)
## options(ChIPseeker.ignore_promoter_subcategory = TRUE)

dir.create("D:/Chipseeker/plot/")
dir.create("D:/Chipseeker/txt/")
dir.create("D:/Chipseeker/plot/TSS-enrich/")
dir.create("D:/Chipseeker/plot/heatmap/")
dir.create("D:/Chipseeker/plot/AnnoBar/")
dir.create("D:/Chipseeker/plot/Annopie/")
dir.create("D:/Chipseeker/plot/vennpie/")
dir.create("D:/Chipseeker/plot/upsetplot/")
dir.create("D:/Chipseeker/plot/Distrib/")
dir.create("D:/Chipseeker/txt/rda/")


## load gtf

x <- Seqinfo(seqnames=c("chr1", "chr10", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrctg1", "chrctg10", "chrctg100", "chrctg101", "chrctg102", "chrctg103", "chrctg104", "chrctg105", "chrctg106", "chrctg107", "chrctg108", "chrctg109", "chrctg11", "chrctg110", "chrctg111", "chrctg112", "chrctg113", "chrctg114", "chrctg115", "chrctg116", "chrctg117", "chrctg118", "chrctg119", "chrctg12", "chrctg120", "chrctg121", "chrctg122", "chrctg123", "chrctg124", "chrctg125", "chrctg126", "chrctg127", "chrctg128", "chrctg129", "chrctg13", "chrctg130", "chrctg131", "chrctg132", "chrctg133", "chrctg134", "chrctg135", "chrctg136", "chrctg137", "chrctg138", "chrctg139", "chrctg14", "chrctg140", "chrctg141", "chrctg142", "chrctg143", "chrctg144", "chrctg145", "chrctg146", "chrctg147", "chrctg148", "chrctg149", "chrctg15", "chrctg150", "chrctg151", "chrctg152", "chrctg153", "chrctg154", "chrctg155", "chrctg156", "chrctg157", "chrctg158", "chrctg159", "chrctg16", "chrctg160", "chrctg161", "chrctg162", "chrctg163", "chrctg164", "chrctg165", "chrctg166", "chrctg167", "chrctg168", "chrctg169", "chrctg17", "chrctg170", "chrctg171", "chrctg172", "chrctg173", "chrctg174", "chrctg175", "chrctg176", "chrctg177", "chrctg178", "chrctg179", "chrctg18", "chrctg180", "chrctg181", "chrctg182", "chrctg183", "chrctg184", "chrctg185", "chrctg186", "chrctg187", "chrctg188", "chrctg189", "chrctg19", "chrctg190", "chrctg191", "chrctg192", "chrctg193", "chrctg194", "chrctg195", "chrctg196", "chrctg197", "chrctg198", "chrctg199", "chrctg2", "chrctg20", "chrctg200", "chrctg201", "chrctg202", "chrctg203", "chrctg204", "chrctg205", "chrctg206", "chrctg207", "chrctg208", "chrctg209", "chrctg21", "chrctg210", "chrctg211", "chrctg212", "chrctg213", "chrctg214", "chrctg215", "chrctg216", "chrctg217", "chrctg218", "chrctg219", "chrctg22", "chrctg220", "chrctg221", "chrctg222", "chrctg223", "chrctg224", "chrctg225", "chrctg226", "chrctg227", "chrctg228", "chrctg229", "chrctg23", "chrctg230", "chrctg231", "chrctg232", "chrctg233", "chrctg234", "chrctg235", "chrctg236", "chrctg237", "chrctg238", "chrctg239", "chrctg24", "chrctg240", "chrctg241", "chrctg242", "chrctg243", "chrctg244", "chrctg245", "chrctg246", "chrctg247", "chrctg248", "chrctg249", "chrctg25", "chrctg250", "chrctg251", "chrctg252", "chrctg253", "chrctg254", "chrctg255", "chrctg26", "chrctg27", "chrctg2729", "chrctg28", "chrctg29", "chrctg3", "chrctg30", "chrctg31", "chrctg32", "chrctg33", "chrctg34", "chrctg35", "chrctg36", "chrctg37", "chrctg38", "chrctg39", "chrctg4", "chrctg40", "chrctg41", "chrctg42", "chrctg43", "chrctg44", "chrctg45", "chrctg46", "chrctg47", "chrctg48", "chrctg49", "chrctg5", "chrctg50", "chrctg51", "chrctg52", "chrctg53", "chrctg54", "chrctg56", "chrctg57", "chrctg58", "chrctg59", "chrctg6", "chrctg60", "chrctg61", "chrctg62", "chrctg63", "chrctg64", "chrctg65", "chrctg66", "chrctg67", "chrctg68", "chrctg69", "chrctg7", "chrctg70", "chrctg71", "chrctg72", "chrctg73", "chrctg74", "chrctg75", "chrctg76", "chrctg77", "chrctg78", "chrctg79", "chrctg8", "chrctg80", "chrctg81", "chrctg82", "chrctg83", "chrctg84", "chrctg85", "chrctg86", "chrctg87", "chrctg88", "chrctg89", "chrctg9", "chrctg90", "chrctg91", "chrctg92", "chrctg93", "chrctg94", "chrctg95", "chrctg96", "chrctg97", "chrctg98", "chrctg99", "chrMt", "chrPt"),
             seqlengths=c(307041717, 150982314, 244442276, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 50531, 992858, 87733, 76596, 101603, 83681, 85134, 61773, 57130, 83503, 57822, 134067, 100537, 57973, 107933, 139473, 62853, 67958, 71809, 9476, 79679, 62538, 43926, 34930, 55817, 63567, 90796, 95420, 44508, 79882, 9160, 46935, 86301, 254428, 70946, 73561, 40576, 83845, 83688, 86441, 57064, 60070, 78370, 40594, 56161, 84781, 344981, 63556, 739442, 51967, 53538, 66893, 69976, 70601, 73006, 58809, 55297, 1738914, 63917, 69060, 65477, 58645, 128268, 60365, 55665, 48311, 77382, 94873, 31799, 50526, 40634, 50620, 66428, 67904, 59352, 53643, 54946, 47089, 70750, 105189, 72536, 43047, 66538, 66228, 69465, 59418, 41679, 340736, 85784, 70926, 80324, 351693, 1278389, 96695, 52129, 56966, 64245, 6454, 82080, 71282, 76840, 59524, 72401, 66116, 526805, 41461, 58299, 138544, 45447, 46041, 60303, 60109, 89036, 54466, 73251, 53313, 82450, 28410, 55164, 73968, 46344, 67134, 58608, 65749, 60995, 55641, 43787, 74799, 74080, 55962, 61539, 58830, 53522, 61449, 24080, 62284, 44889, 80070, 63290, 55301, 64997, 53614, 41417, 17651, 85258, 77423, 65674, 73141, 52477, 94458, 36081, 70135, 54117, 7104, 61153, 95256, 70597, 45166, 52360, 56713, 51365, 48860, 76559, 47028, 81363, 66018, 30460, 41524, 63031, 51883, 52549, 82972, 69595, 46945, 215635, 11339, 124116, 65819, 66393, 59657, 110255, 846561, 54684, 108212, 39815, 74491, 70486, 48545, 91381, 29460, 66261, 84148, 66656, 67687, 45711, 101180, 57348, 65659, 33120, 72581, 100579, 83265, 381258, 83120, 756422, 458256, 82446, 695048, 62662, 65407, 53907, 255484, 117437, 76354, 92234, 63530, 51910, 666287, 69959, 60926, 60205, 65268, 71433, 54868, 77242, 342161, 83778, 81536, 516616, 5568, 87707, 93364, 73283, 64470, 65470, 80451, 83476, 67579, 98601, 109991, 47381, 57632, 58721, 82793, 110418, 52427, 795190, 94431, 58809, 55397, 76262, 44176, 90109, 103641, 36676, 569630, 140384),
             isCircular=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
             genome="zeamaysV4")

cornseqinfo <- makeTxDbFromGFF("E:/index/ZeamaysB73.gtf", format=c("gtf"),chrominfo= x)
readRDS("E:/index/maizev4gtfchr.rds")

## Profile of ChIP peaks binding to TSS regions

promoter <- getPromoters(TxDb=cornseqinfo, upstream=3000, downstream=3000)

##recycle peak-files

files <- list.files(path="E:/DAP", pattern="*.GEM_events.narrowPeak", full.names=TRUE, recursive=FALSE)

for (file in files) {

peak <- readPeakFile(file) ## 读入bed文件
fn   <- basename(file)
try( {
tagMatrix <- getTagMatrix(peak, windows=promoter) ## Profile of ChIP peaks binding to TSS regions


## plot

pdf(paste("E:/DAP/plot/TSS-enrich/",fn,".TSS-enrich.pdf",sep=""))
print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic (5'->3')", ylab = "Read Count Frequency"))
dev.off()

pdf(paste("E:/DAP/plot/heatmap/",fn,".TSS-heatmap.pdf",sep=""))
print(tagHeatmap(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic (5'->3')", ylab = "Read Count Frequency", color="blue4"))
dev.off()

## annotate peak

peakannotate <- annotatePeak(peak, tssRegion=c(-3000, 0), TxDb=cornseqinfo, addFlankGeneInfo=TRUE, flankDistance=50000) ## peak的注释
as.GRanges(peakannotate)
head(as.GRanges(peakannotate))

## plot

pdf(paste("E:/DAP/plot/AnnoBar/",fn,".AnnoBar.pdf",sep=""))
print(plotAnnoBar(peakannotate)) ## 可视化 Pie and Bar plot
dev.off()

pdf(paste("E:/DAP/plot/Annopie/",fn,".Annopie.pdf",sep=""))
print(plotAnnoPie(peakannotate)) ## 可视化 Pie and Bar plot
dev.off()

pdf(paste("E:/DAP/plot/vennpie/",fn,".vennpie.pdf",sep=""), width =13, height=7)
print(vennpie(peakannotate)) ## 可视化 Pie and Bar plot
dev.off()

pdf(paste("E:/DAP/plot/upsetplot/",fn,".upsetplot.pdf",sep=""))
print(upsetplot(peakannotate)) ## 可视化 Pie and Bar plot
dev.off()

pdf(paste("E:/DAP/plot/Distrib/",fn,".Distrib.pdf",sep=""))
print(plotDistToTSS(peakannotate, title="Distribution of binding loci relative to TSS"))
dev.off()

## save annotate files
save(peakannotate,file=(paste("E:/DAP/txt/rda/",fn,".rda",sep=""))) # Output peakAnnolist file

write.table(as.data.frame(peakannotate),file=paste("E:/DAP/txt/",fn,".txt",sep=""),sep='\t',quote = F,row.names=FALSE) # Output peakAnnolist file
})

}
##完成

