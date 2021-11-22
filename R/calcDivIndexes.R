#' calcDivIndexes
#'
#' Calculates common diversity indices per sample in a dataframe of counts
#'
#' @title
#' Calculate diversity indices for barcode samples
#'
#' @description
#' Takes a dataframe of barcode counts and computes shannon,
#' simpson, inverse simpson and gini coefficients for each
#' sample
#'
#' @param dat dataframe of raw or normalised counts with samples as columns and observations/barcodes as rows
#'
#' @return Returns a data-frame containing the calculated diversity indices per sample
#' @importFrom vegan diversity
#' @importFrom ineq Gini
#' @export
#' @examples
#' indexes <- calcDivIndexes(df = counts.df)

dat <- dge.collapsed$counts

calcDivIndexes <- function(dat){
  colnames <- colnames(dat)
  df <- data.frame(name = c(), shannon = as.numeric(), simpson = as.numeric(),
                   invsimpson = as.numeric(), gini = as.numeric())
  for (i in 1:length(colnames(dat))) {
    name <- colnames(dat)[i]
    shannon <- diversity(dat[, i], index = "shannon")
    simpson <- diversity(dat[, i], index = "simpson")
    invsimpson <- diversity(dat[, i], index = "invsimpson")
    gini <- Gini(dat[, i])
    x <- data.frame(name, shannon, simpson, invsimpson, gini)
    df <- rbind(df, x)
  }
  return(df)
}
