#' Plot barcode frequency distribution
#'
#' Plot barcode frequency distribution from a two column data frame containing barcode name and count.
#'
#' @param barcodes DGEList or DESeqDataSet object containing raw or normalised barcode counts with replicate grouping info in object metadata
#' @param samplename a character vector of containing grouping information. Must be as long as the columns of object. Can pass a metadata column from object.
#'
#' @return Returns a plot of the barcode frequency distribution for a sample
#' @export
#'
#' @examples
#' plotBarcodeFrequency(barcodes, samplename = "test sample")
#'

plotBarcodeFrequency <- function(barcodes, samplename){
  colnames(barcodes) <- c("Barcode", "Raw_count")
  barcodes$Raw_count <- as.numeric(barcodes$Raw_count)

  p <- ggplot(barcodes, aes(y=Raw_count, x=seq(1,length(rownames(barcodes))))) +
    geom_point(stat = "identity", show.legend = F) + scale_y_continuous(trans ='log10') +
    theme_bw() +
    scale_size_manual(values=c(2,2)) +
    geom_hline(yintercept = 100, color = "grey40") +
    xlab("Barcode") +
    ylab("Count") +
    ggtitle(paste(samplename, ": Barcode frequency distribution"))

  return(p)
}
