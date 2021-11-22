#' plot_proportional_timeseries
#'
#' Generate proportional timeseries plot from raw / normalised barcode count object.
#'
#' @param counts.obj dataframe containing raw counts of barcodes
#' @param outfile desired name of output plot
#' @param name desired plot title
#' @param seed RNG seed
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' plot_proportional_timeseries(counts.obj, name = "Proportional Bubble Plot", seed = 5)

plot_proportional_timeseries <- function(counts.obj, outfile, name = "", seed = 5){
    require(ggplot2)
    require(reshape2)
    barcodes.proportional <- counts.obj
    barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100

    barcodes.proportional.melted <- melt(barcodes.proportional)
    head(barcodes.proportional.melted)
    colnames(barcodes.proportional.melted) <- c("Barcode", "Sample", "Proportion")
    barcodes.proportional.melted$Barcode <- as.factor(barcodes.proportional.melted$Barcode)
    barcodes.proportional.melted$Proportion <- as.numeric(barcodes.proportional.melted$Proportion)

    timepoints <- unique(barcodes.proportional.melted$Sample)

    colors <- c(RColorBrewer::brewer.pal(12, "Set3"),RColorBrewer::brewer.pal(9, "Set1"))
    set.seed(seed) # set custom seed to get same color order every time
    colors <- sample(colors, length(rownames(barcodes.proportional.melted)), replace = TRUE)
    barcodes.proportional.melted$color <- colors

    #head(barcodes.proportional.melted)
    # inspired by genBaRcode package
    timeseries.plot <- ggplot(barcodes.proportional.melted,
                              aes_string(group = "Barcode",fill = "Barcode", x = "Sample", y = "Proportion")) +
        theme_bw() +
        geom_area(alpha=0.9) +
        scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
        scale_x_discrete(breaks = timepoints, labels = timepoints, limits = timepoints) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(label = paste("Proportional Timeseries Plot:", name))

    timeseries.plot
}
