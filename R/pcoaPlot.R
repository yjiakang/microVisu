

#' Title
#'
#' @param distTab distance table: bray curtis, weighted/unweighted unifrac distance table
#' @param metaData design file
#' @param distType which kind of distance you use for plot
#' @param classForColor collor group
#' @param classForShape shape group
#' @param col colour palette: including all the types of the "display.brewer.all()" in the RColorBrewer package
#'
#' @return
#' @export
#'
#' @examples bray <- system.file("extdata", "bray_curtis_otu_table_css", package = "microVisu")
#' @examples design.txt <- system.file("extdata", "design.txt", package = "microVisu")
#' @examples PCoAPlot(bray, meta, distType = "Bray Curtis", classForColor = "climate",classForShape = "status",
#'  col = "Paired")
pcoaPlot <- function(distTab,
                     metaData,
                     distType,
                     classForColor,
                     classForShape = "None",
                     col = "Set3") {
    library(vegan)
    library(ggplot2)
    design <- read.table(metaData, header = T, row.names = 1, sep = "\t")
    distTab <- read.table(distTab, sep = "\t", header = T, check.names = F)
    idx <- rownames(design) %in% colnames(distTab)
    sub_design <- design[idx,]
    distTab <- distTab[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix
    pcoa <- cmdscale(distTab, k = 3, eig = T) # k is dimension, 3 is recommended
    points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
    colnames(points) <- c("x", "y", "z")
    eig <- pcoa$eig
    points <- cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
    if(classForShape == "None") {
        ggplot(points, aes(x = x, y = y, color = !!sym(classForColor))) +
            geom_point(alpha = .7, size = 3) +
            labs(x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
                 y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""),
                 title = paste(distType, "PCoA")) +
            scale_color_brewer(palette = col)
    } else {
        ggplot(points, aes(x = x, y = y, color = !!sym(classForColor), shape = !!sym(classForShape))) +
            geom_point(alpha = .7, size = 3) +
            labs(x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
                 y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""),
                 title = paste(distType, "PCoA")) +
            scale_color_brewer(palette = col)
    }
}

