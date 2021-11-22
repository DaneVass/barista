#' Plot sample read counts
#'
#' Simple plot of total read counts per sample
#'
#' @param counts data.frame of barcode count x sample
#' @param group a character vector of containing grouping information. Must be equal the number of columns of counts. Can pass a metadata column from DGEList object.
#' @param log10 Boolean. log10 transform output?
#' @param legend Boolean. Include legend?
#' @param order Boolean. Order samples by group?
#'
#' @return Returns a plot of the read counts per column (sample) in a data frame
#' @export
#'
#' @examples
#' plotBarcodeFrequency(barcodes, samplename = "test sample")
#'

plotReadCounts <- function(counts, group, log10 = FALSE, legend = TRUE, order = TRUE){

  cols <- as.factor(group)
  dat <- colSums(counts)
  dat <- data.frame(sample = names(dat), counts = dat, group = factor(cols), row.names = NULL)

  if(isTRUE(order)){
    dat$sample <- factor(dat$sample, levels = dat$sample[order(dat$group, decreasing = T)])
  }

  # plot data
  p <- ggplot(data=dat, aes(y=sample, x=counts, fill = group)) +
    geom_bar(stat="identity", width = .75) +
    theme_bw() +
    theme(axis.text.y = element_text(angle=0, vjust=0, hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 5)) +
    ggtitle("Total Read Counts") +
    xlim(0,max(dat$counts+100)) +
    scale_fill_manual(values = rev(ggpubr::get_palette("npg", length(unique(dat$group)))))

  if(isTRUE(log10)){
    p <- p + scale_x_log10()
  }

  if(legend == FALSE){
    p <- p + theme(legend.position = "none")
  }

  print(p)
}
