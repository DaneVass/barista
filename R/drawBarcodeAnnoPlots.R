#' Draw Barcode Annotation Plots
#'
#' Plot summary statistics resulting from 10X barcode annotation
#'
#' @param bam.parsed.df Dataframe from mapBarcodeReads() containing annotated barcode reads
#' @param breaks Optional. Bin size for histograms
#' @param outdir Optional. Output directory. Will generate a plots directory inside here to store plots
#' @param prefix Optional. Prefix to add to filename for each plot
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence mapped to the lintraceR DNA barcode ID.
#' @export
#' @examples
#' drawBarcodeAnnoPlots(bam.parsed.df, breaks = 50, outdir = "./", prefix = "barcode_anno_plots")
drawBarcodeAnnoPlots <- function(bam.parsed.df, breaks = 50, outdir = "./", prefix = "barcode_anno_plots"){
    plots.dir <- file.path(outdir,"plots")
    dir.create(plots.dir, showWarnings = F)

    # plot cell ID frequency
    pdf(file.path(plots.dir, "barcodes_per_10X_cellID.pdf"))
    hist(table(bam.parsed.df$Cell.10X.Barcode),
         col = "steelblue",
         breaks = breaks,
         main = paste("Distribution of 10X cell ID tags in BAM entries", start, "to", end))
    dev.off()

    pdf(file.path(plots.dir, "barcode_match_lengths.pdf"))
    lengths <- sapply(X = bam.parsed.df$DNA.Barcode, FUN = nchar, simplify = T)
    hist(lengths,
         col = "steelblue",
         breaks = breaks,
         main = paste("Distribution of barcode match lengths in BAM entries", start, "to", end))
    dev.off()

}

