#' nthPercentile_barcodes
#'
#' Calculate barcodes comprising the Nth percentile for each sample & generate cumulative sum plots
#'
#' @param counts.obj dataframe containing raw or normalised barcode counts
#' @param percentile desired percentile value
#' @return Returns a cumulative sum & percentile plot per sample.
#' @export
#' @examples
#' nthPercentile_barcodes(counts, percentile = .95)

nthPercentile_barcodes <- function(counts, percentile = .95){
  #source("scripts/findElbow.R")

  counts.obj <- as.data.frame(counts)
  dim(counts.obj)
  samples <- colnames(counts.obj)

  # setup empty dataframes to collect cumulative sum data
  percentile.df <- data.frame()
  elbow.df <- data.frame()

  # setup colors for plotting
  colors <- c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
  colors <- rep(colors,times=length(rownames(counts.obj))/length(colors)+1)
  counts.obj$color <- as.character(colors[1:length(rownames(counts.obj))])

  iter = 1
  for (i in samples){

    if(i == 'color'){
      message("skipping color column")
    }

    else {
      print(i)
      # sort dataset
      sorted <- counts.obj[,i, drop = F]
      sorted <- sorted[order(sorted, decreasing = T),i, drop = F]
      sorted <- sorted[sorted > 0, , drop = F]

      dim(sorted)
      head(sorted, 10)

      colsum <- sum(sorted)
      percentile.cutoff <- sum(sorted)*percentile

      # generate cumulative sum of sorted dataset
      cumsum <- cumsum(sorted)
      max(cumsum) == colsum

      # find number of barcodes that make up Nth percentile
      len <- length(which(cumsum <= percentile.cutoff))
      print(paste(i, " ", percentile,"th percentile: ", len))

      # fill Nth percentile data frame for plotting below
      d <- data.frame(Sample=factor(i),Barcodes=len)
      percentile.df <- rbind(percentile.df, d)
      cum.sum <- cumsum(sorted[,1])

      # plot cumulative sum showing 90th percentile and elbow points
      #pdf(paste("plots/cumulative_plots/",i,"_Cumulative_plot.pdf", sep=''))
      plot(cum.sum, main = paste("Cumulative Plot: ", i))
      #abline(v=len, col = "grey60")
      #abline(v=elbow.point, col = "grey40")
      abline(h = percentile.cutoff, col = "grey60")
      elbow.point <- findElbow(cum.sum, plot = F, returnIndex = T)
      points(x = elbow.point ,y = cum.sum[elbow.point], col = "tomato3", pch = 19)

      points(x = len ,y = percentile.cutoff, col = "dodgerblue", pch = 19)
      #dev.off()

      # fill under elbow data frame for plotting below
      under.elbow <- length(which(cum.sum <= cum.sum[elbow.point]))
      d <- data.frame(Sample=factor(i),Barcodes=under.elbow)
      elbow.df <- rbind(elbow.df, d)
      print(paste(i," Under elbow: ", under.elbow))

      # extract barcodes in 90th percentile
      rows <- rownames(sorted)[which(cumsum <= percentile.cutoff)]
      length(rows)
      sorted.top <- counts.obj[rows,c(i,'color'),drop = F]
      top.bc.sum <- sum(sorted.top[,1])
      sorted.top$proportions <- (sorted.top[,1]/top.bc.sum)*100
      #print(head(sorted.top))

      #barcodes.90.allsamps.Rep1[[iter]] <- rows
      #names(barcodes.90.allsamps.Rep1[iter]) <- i

      # plot proportional plot for the barcodes in 90th percentile
      #pdf(paste("plots/barcode_plots/",i,"_barcode_plot.pdf", sep=''), height = 10, width = 3)
      p <- barplot(as.matrix(rev(sorted.top$proportions)), beside = F, horiz = F, border = T,
                   col = rev(sorted.top$color), las=2, ylab = "Barcode proportion (%)", main = paste("Barcode plot :", i))
      #dev.off()
      iter = iter + 1
    }
  }
}







