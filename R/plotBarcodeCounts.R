#' Plot total counts per barcode in library
#'
#' Simple plot of total read counts per barcode in library
#'
#' @param counts data.frame of barcode count x sample
#'
#' @return Returns a plot of the read counts per barcode (row) in a data frame
#' @export
#'
#' @examples
#' plotBarcodeCounts(counts)
#'

plotBarcodeCounts <- function(counts, order = F, log10 = F){
  rowsums <- rowSums(counts)

  if(isTRUE(log10)){
    rowsums <- log10(rowsums)
  }

  if(isTRUE(order)){
    ordered <- sort(rowsums, decreasing = T)
    barplot(ordered,
            las=2,
            main="Total counts per barcode",
            axisnames = F,
            cex.axis=0.8,
            xlab = "Barcode - Descending order by total read count across samples", ylab = "Barcode total read count")
  } else {
    barplot(rowsums,
            las=2,
            main="Total counts per barcode",
            axisnames = F,
            cex.axis=0.8,
            xlab = "Barcode - Descending order by frequency in reference library", ylab = "Barcode total read count")

  }
}
