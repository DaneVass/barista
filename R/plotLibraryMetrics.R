

# Barcode frequency distribution
plot.barcode.dist <- function(barcodes, samplename){
  ggplot(barcodes, aes(y=Raw_count, x=seq(1,length(rownames(barcodes))))) +
    geom_point(stat = "identity", show.legend = F) + scale_y_continuous(trans='log10') +
    theme_bw() +
    scale_size_manual(values=c(2,2)) +
    geom_hline(yintercept = 100, color = "grey40") +
    xlab("Barcode") +
    ylab("Count") +
    ggtitle(paste(samplename, ": Barcode frequency distribution"))
  ggsave(paste("plots/", samplename, "_barcode_frequency_distribution.png", sep = ''), dpi = 300)
}

# Barcode cumulative sum plot
plot.barcode.cumsum <- function(barcodes, samplename){
  sorted <- sort(barcodes$Raw_count, decreasing = T)
  barcodes <- barcodes[sorted,]
  colsum <- sum(sorted)
  colsum.90 <- sum(sorted)*0.9
  colsum.99 <- sum(sorted)*0.99
  cumsum <- cumsum(sorted)
  cumsum <- cumsum/max(cumsum)
  length(cumsum)
  b <- seq(1,length(cumsum),1)
  b <- b/length(cumsum)
  d <- data.frame(barcode=b, proportion=cumsum)
  dim(d)

  ggplot(d, aes(y=proportion, x=barcode)) +
    geom_point(stat = "identity", show.legend = F) +
    theme_bw() +
    scale_size_manual(values=c(2,2)) +
    xlab("Barcode rank (Proportion)") +
    ylab("Cumulative Sum (Proportion)") +
    ggtitle(paste(samplename, ": Barcode cumulative sum"))
  ggsave(paste("plots/",samplename, "_cumulative_sum_distribution.png", sep = ''), dpi = 300)
}
