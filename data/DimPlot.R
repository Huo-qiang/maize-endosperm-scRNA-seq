DimPlot <- function(cell.info, dim.1="UMAP_1", dim.2="UMAP_2", cell.type=NULL, regulon=NULL) {
  if (!is.null(cell.type)){
    data <- cell.info[, c(dim.1, dim.2, "CellType")]
    data$pt.col <- ifelse(data$CellType == cell.type, "red", "#DFDFDF")
    data$pt.size <- ifelse(data$CellType == cell.type, 0.2, 0.1)
    title <- paste0(cell.type)
    col.title <- "red"
  } else {
    data <- cell.info[, c(dim.1, dim.2, regulon)]
    data$pt.col <- ifelse(data[, regulon], "#006464", "#DFDFDF")
    data$pt.size <- ifelse(data[, regulon], 0.2, 0.1)
    title <- paste0("Regulon: ", regulon)
    col.title = "#006464"
  }
  ggplot(data, aes(get(dim.1), get(dim.2))) + 
    geom_point(size=data$pt.size, color=data$pt.col) + 
    theme_bw(base_size = 12) + 
    ggtitle("") + 
    xlab(dim.1) + ylab(dim.2) + 
    annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,label=title,color=col.title,size=6) + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
          )
}