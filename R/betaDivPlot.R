
#' Title
#'
#' @param otuTab    otu table of your sample
#' @param metaData  design file
#' @param classToPlot   which column you want to plot
#' @param varAxes   whether to plot the contributive variables (TRUE/FALSE)
#' @param col   colour palette: including all the types of the "display.brewer.all()" in the RColorBrewer package
#'
#' @return
#' @export
#' @examples otu_table_L2.txt <- system.file("extdata", "otu_table_L2.txt", package = "microVisu")
#' @examples design.txt <- system.file("extdata", "design.txt", package = "microVisu")
#' @examples betaDivPlot(otuTab = otu_table_L2.txt, metaData = design.txt,
#' classToPlot = "status", col = "Set3")
betaDivPlot <- function(otuTab, metaData, classToPlot, varAxes = FALSE, col = "Set3") {
    # load packages needed
    library("ggbiplot")
    otuTab <- read.delim(otuTab, header = TRUE, sep = "\t") # Import otu table
    design <- read.table(metaData, header = T, row.names = 1, sep = "\t")
    idx <- intersect(colnames(otuTab), rownames(design))
    sub_design = design[idx, ]
    otuTab <- otuTab[, idx]
    otu.pca <- prcomp(t(otuTab), scale. = TRUE)
    groups <- sub_design[,intersect(classToPlot, colnames(sub_design))]
    ggbiplot(otu.pca, obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, var.axes = varAxes)+
        scale_color_brewer(palette = col, name = classToPlot)
}
