#' Visualize the amplicon data
#'
#' @param otuTab otu table of your sample
#' @param metaData design file
#' @param classToPlot which column you want to plot
#' @param topNum  top n taxa to plot
#' @param col colour palette: including all the types of the "display.brewer.all()" in the RColorBrewer package
#' @param classToFacet which class you want for facet, default none
#' @param legCol column number of legend, default 2. If text is too long, suggestion is 1
#' @return
#'
#' @export 
#'
#' @examples otu_table_L2.txt <- system.file("extdata", "otu_table_L2.txt", package = "microVisu")
#' @examples design.txt <- system.file("extdata", "design.txt", package = "microVisu")
#' @examples taxBarPlot(otuTab = otu_table_L2.txt, metaData = design.txt,
#'  classToPlot = "status", topNum = 10, col = "Set3", classToFacet = "knownseverity")
taxBarPlot  <- function(otuTab, metaData, classToPlot, topNum, col, classToFacet = FALSE, legCol = 2) {
    # load packages needed
    library("tidyr")
    library("ggplot2")
    library("reshape2")
    library("RColorBrewer")
    library("dplyr")
    # color num in the palette
    if(col %in% (c("Set2", "Pastel2", "Dark2", "Accent"))) {
        ncolor = 8
    } else {
        ncolor = 9
    }
    otuTab <- read.delim(otuTab, header = TRUE, sep = "\t") # Import otu table
    otuTab <- as.data.frame(t(t(otuTab)/colSums(otuTab)*100)) # Tranfer to percent
    metaData <- read.table(metaData, header = TRUE, row.names = 1, sep = "\t") # Import metadata table
    idx <- intersect(rownames(metaData),colnames(otuTab)) # Find the common samples both in metadata and otu table
    metaData <- metaData[idx,]
    otuTab <- otuTab[,idx]
    if(classToFacet == FALSE) {
        samFile <- as.data.frame(metaData[,classToPlot], row.names = row.names(metaData))
        colnames(samFile)[1] <- classToPlot
        otuTab <- merge(samFile, t(otuTab), by = "row.names")[, -1]
        otuTabMean <- aggregate(otuTab[,-1], by = otuTab[1], FUN = mean) # Calculate the mean of the same group
        otuTabMeanFinal <- do.call(rbind, otuTabMean)[-1, ]
        colnames(otuTabMeanFinal) <- otuTabMean[, classToPlot]
        otuRowName <- row.names(otuTabMeanFinal)
        otuTabMeanFinal <- as.data.frame(otuTabMeanFinal)
        otuTabMeanFinal <- as.data.frame(lapply(otuTabMeanFinal, as.numeric))
        rownames(otuTabMeanFinal) <- otuRowName
        otuTabMeanFinal$total <- apply(otuTabMeanFinal, 1, sum)
        otuTabMeanFinal <- dplyr::arrange(otuTabMeanFinal, desc(total)) # Sort based on the total counts using the imported pkg
        otuTabMeanFinal <- subset(head(otuTabMeanFinal, n = topNum), select = -total)
        otuTabMeanFinal %<>% 
            t %>% 
            as.data.frame %>% 
            mutate(Others = 100 - rowSums(.)) %>% 
            t %>% 
            as.data.frame %>% 
            mutate(taxa = rownames(.))
        dataForPlot <- otuTabMeanFinal %>% gather(!!sym(classToPlot), abundance, -taxa) # Change into long data
        ggplot(dataForPlot, aes(x = !!sym(classToPlot), y = abundance, fill = taxa)) +
            geom_bar(stat = "identity", width = 0.5) +
            xlab(NULL) +
            scale_fill_manual(values =  colorRampPalette(brewer.pal(ncolor, col))(topNum+1)) +
            theme_bw() +
            theme(axis.title = element_text(size = 10, face = "bold"),
                  axis.text.x = element_text(size = 10, face = "bold"),
                  legend.position = "bottom") +
            guides(fill = guide_legend(ncol = legCol)) + 
            labs(fill = "Taxonomy") +
            ylab("Abundance(%)")
    } else {
        otuTab$sum <- rowSums(otuTab)
        otuTab <- arrange(otuTab, desc(sum))
        otuMeta <- t(subset(head(otuTab, n = topNum), select = -sum)) %>% 
            as.data.frame %>% 
            mutate(Others = 100 - rowSums(.)) %>% 
            bind_cols(metaData)
        otuGp <- aggregate(otuMeta[1:(topNum + 1)], 
                           by = otuMeta[c(which(colnames(otuMeta)==classToPlot), 
                                          which(colnames(otuMeta)==classToFacet))],
                           FUN = mean)
        otuLong <- melt(otuGp)
        ggplot(otuLong, aes(x = !!sym(classToPlot), y = value, fill = variable)) +
            geom_bar(stat = "identity", width = 0.5) +
            xlab(NULL) +
            scale_fill_manual(values =  colorRampPalette(brewer.pal(ncolor, col))(topNum + 1)) +
            theme_bw() + 
            theme(axis.title = element_text(size = 10, face = "bold"),
                  axis.text.x= element_text(size = 10, face = "bold"),
                  legend.position = "bottom")+
            guides(fill = guide_legend(ncol = legCol)) +
            labs(fill = "Taxonomy") +
            facet_grid(as.formula(paste(".", "~", classToFacet))) + 
            ylab("Abundance(%)")
    }    
}

