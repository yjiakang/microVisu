

#' Title PERMANOVA
#'
#' @param otuTab otu table of your sample
#' @param metaData design file
#' @param var variables
#' @param method "manhattan", "euclidean", "canberra", "clark",
#'  "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita",
#'   "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
#'    "chisq" or "chord".
#' @param permutations permutation times
#'
#' @return
#' @export
#'
#' @examples otu_table_L2.txt <- system.file("extdata", "otu_table_L2.txt", package = "microVisu")
#' @examples design.txt <- system.file("extdata", "design.txt", package = "microVisu")
#' @examples permanovaTest(otu_table_L2.txt, design.txt, "status*cultivars")
permanovaTest <- function(otuTab, metaData, var,
                         method = "bray", permutations = 999) {
    otuTab <- read.delim(otuTab, header = TRUE, sep = "\t") # Import otu table
    design <- read.table(metaData, header = T, row.names = 1, sep = "\t") # Import metadata info
    idx <- intersect(colnames(otuTab), rownames(design)) # Intersect design and otu
    sub_design <-  design[idx, ]
    otuTab <- otuTab[, idx]
    adonis_result <- vegan::adonis2(as.formula(paste("t(otuTab)", "~", var)), sub_design, method = method, permutations = 999) # Adonis test
    adonis_result
}

