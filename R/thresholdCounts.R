#' Threshold counts
#'
#' Threshold dataframe to a given level and return number of barcodes meeting threshold in each
#' sample.
#'
#' @param df Dataframe to be thresholded.
#' @param threshold The threshold to use. Rows with count below this will be removed.
#' @param plot Logical. Draw plots of dataset?
#' @return Returns a thresholded data-frame
#'
#' @export
#' @examples
#' thresholdCounts(df, threshold = 20, plot = FALSE)

thresholdCounts <- function(df, threshold = 20, plot = FALSE){

  if(!isnull(threshold)){
    threshold <- threshold
  } else {
    message("No threshold given, defaulting to 20.")
    threshold <- 20
  }

  # setup output dataframe
  above.threshold.counts <- data.frame(Sample=factor(), Count=c())

  for(sample in colnames(df)){
    above.threshold = length(which(df$counts[,sample] >= threshold))
    d <- data.frame(Sample=factor(sample),Count=above.threshold)
    above.threshold.counts <- rbind(above.threshold.counts, d)
    print(paste(sample, sum(above.threshold)))
  }

  if(plot == FALSE){
    return(above.threshold.counts)
  } else {
    g <- ggplot(above.threshold.counts[which(grepl("BR1", above.threshold.counts$Sample)==TRUE),],
                aes(x=Sample,y=Count))
    g + geom_bar(stat = "identity") +
      theme(panel.grid.major.x=element_line(colour="grey70")) +
      labs(title = paste("Number of barcodes present per conditon. Threshold:", threshold)) +
             xlab("Condition") + ylab("Total number of barcodes") +
             theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}
