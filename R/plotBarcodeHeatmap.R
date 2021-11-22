#' plot_proportional_heatmap
#'
#' Generate proportional heatmap from raw count / normalised count object.
#'
#' @param counts.obj dataframe containing raw / normalised counts of barcodes
#' @param outfile desired name of output plot
#' @param name desired plot title
#' @param scale scaling to use for heatmap: can be none, row or column as per pheatmap
#' @return Returns a plot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' plot_proportional_heatmap(counts.obj, outfile, name = "Proportional Barcode Plot", scale = "row")

plot_proportional_heatmap <- function(counts.obj, outfile, name = "Proportional Barcode Plot", scale = "row"){
    require(ggplot2)
    require(reshape2)
    barcodes.proportional <- counts.obj

    # convert counts to proportions
    barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100
    str(barcodes.proportional)
    barcodes.proportional <- as.data.frame(barcodes.proportional)

    # set colors

    pheatmap(barcodes.proportional, scale = scale)
            #beside = F, horiz = T, border = T, col = barcodes.proportional$color,
            #las=2, xlab = "Barcode proportion (%)", main = name)

    #dev.off()
}
