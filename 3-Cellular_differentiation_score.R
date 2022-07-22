library(Seurat)
library(ggplot2)
library(UCell)
##读取6dap和8dap的marker基因
maize_sub <- readRDS("maize_sub.rds")
dap6markergene<-read.csv("day6marker.csv")
dap8markergene<-read.csv("day5marker.csv")
