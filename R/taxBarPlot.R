#' Visualize the amplicon data
#'
#' @param otuTab otu table of your sample
#' @param metaData design file
#' @param classToPlot which column you want to plot
#' @param topNum  top n taxa to plot
#' @param col colour palette: including all the types of the "display.brewer.all()" in the RColorBrewer package
#' @return
#'
#' @export 
#'
#' @examples otu_table_L2.txt <- system.file("extdata", "otu_table_L2.txt", package = "microVisu")
#' @examples design.txt <- system.file("extdata", "design.txt", package = "microVisu")
#' @examples taxBarPlot(otuTab = otu_table_L2.txt, metaData = design.txt,
#'  classToPlot = "status", topNum = 10, col = "Set3")
taxBarPlot  <- function(otuTab, metaData, classToPlot, topNum, col) {
    # load packages needed
    library("tidyr")
    library("ggplot2")
    otuTab <- read.delim(otuTab, header = TRUE, sep = "\t") # Import otu table
    otuTab <- as.data.frame(t(t(otuTab)/colSums(otuTab)*100)) # Tranfer to percent
    metaData <- read.table(metaData, header = TRUE, row.names = 1, sep = "\t") # Import metadata table
    idx <- intersect(rownames(metaData),colnames(otuTab)) # Find the common samples both in metadata and otu table
    metaData <- metaData[idx,]
    otuTab <- otuTab[,idx]
    samFile <- as.data.frame(metaData[,classToPlot], row.names = row.names(metaData))
    colnames(samFile)[1] <- classToPlot
    otuTab <- merge(samFile, t(otuTab), by = "row.names")[, -1]
    otuTabMean <- aggregate(otuTab[,-1], by = otuTab[1], FUN = mean) # Calculate the mean of the same group
    otuTabMeanFinal <- do.call(rbind, otuTabMean)[-1, ]
    colnames(otuTabMeanFinal) <- otuTabMean[, classToPlot]
    otuTabMeanFinal <- as.data.frame(otuTabMeanFinal)
    otuTabMeanFinal$total <- apply(otuTabMeanFinal, 1, sum)
    otuTabMeanFinal$taxa <- rownames(otuTabMeanFinal)
    otuTabMeanFinal <- dplyr::arrange(otuTabMeanFinal, desc(total)) # Sort based on the total counts using the imported pkg
    otuTabMeanFinal <- subset(head(otuTabMeanFinal, n = topNum), select = -total)
    dataForPlot <- otuTabMeanFinal %>% gather(classToPlot, abundance, -taxa) # Change into long data
    ggplot(dataForPlot, aes(x = classToPlot, y = abundance, fill = taxa)) +
        geom_bar(stat = "identity", width = 0.5) +
        scale_fill_brewer(palette = col) +
        xlab(NULL) +
        theme(axis.title = element_text(size = 10, face = "bold"),
              axis.text.x= element_text(size = 10, face = "bold"))+
        labs(fill = "Taxonomy") +
        ylab("Abundance(%)")
}
