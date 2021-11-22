#' proportionalBubbleplot
#'
#' Generate proportional bubbleplots from raw count object with barcodes labelled above a specified threshold
#'
#' @param counts.obj dataframe containing raw counts of barcodes. assumes barcodes as rownames.
#' @param labels Boolean. print barcode labels?
#' @param proportion.cutoff barcodes represented at a percentage within any sample above this threshold will be labelled
#' @param name desired plot title
#'
#' @return Returns a bubbleplot of barcodes represented by proportion of total pool
#' @export
#' @examples
#' proportionalBubbleplot(counts.obj, name = "Proportional Bubble Plot", label = T, proportion.cutoff = 10)

proportionalBubbleplot <- function(counts.obj, labels = T, proportion.cutoff = 10, name = "Proportional Bubble Plot"){
  # transform CPM into percentage within sample
  barcodes.proportional <- as.data.frame(counts.obj)
  barcodes.proportional <- sweep(barcodes.proportional,2,colSums(barcodes.proportional),`/`) * 100

  # give each barcode a specific color
  colors <- scales::hue_pal()(30)
  colors <- sample(colors, length(rownames(barcodes.proportional)), replace = TRUE)
  names(colors) <- rownames(barcodes.proportional)
  barcodes.proportional$Color <- colors

  # maintain Rank information of barcodes. Plot in ascending rank order from original barcode library
  barcodes.proportional$Position <- as.factor(seq(1,length(rownames(counts.obj))))
  barcodes.proportional$Barcode <- rownames(barcodes.proportional)

  # melt data frame and rename columns correctly
  barcodes.proportional.melted <- reshape2::melt(barcodes.proportional, id.vars = c("Color", "Position", "Barcode"))
  colnames(barcodes.proportional.melted) <- c("Color", "Position", "Barcode", "Sample", "Proportion")
  # convert variables to correct form
  barcodes.proportional.melted$Barcode <- as.factor(barcodes.proportional.melted$Barcode)
  barcodes.proportional.melted$Sample <- as.factor(barcodes.proportional.melted$Sample)
  barcodes.proportional.melted$Proportion <- as.numeric(barcodes.proportional.melted$Proportion)

  if(isTRUE(labels)){
    # identify high proportion barcodes for labelling
    Highbarcodes <- barcodes.proportional.melted[barcodes.proportional.melted$Proportion > proportion.cutoff,]
    # only take unique barcode labels
    HighbarcodeOrdered <- unique(Highbarcodes[order(Highbarcodes$Barcode),1:3])

    # generate bubbleplot
    bubble.plot <- ggplot(barcodes.proportional.melted, aes(x=Position, y=Sample, size = Proportion, color=Color)) +
      geom_point(stat = "identity", alpha = 0.6, shape = 16) +
      scale_color_identity() +
      labs(y = "Condition", x = "Barcode", title = name) +
      scale_size_continuous(range = c(0.1, 10)) +
      scale_x_discrete(breaks = HighbarcodeOrdered$Position, labels = HighbarcodeOrdered$Barcode) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle=90, vjust=0.5, colour = HighbarcodeOrdered$Color, size = 5),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(size=8),
            axis.title = element_text(size= 6))
  } else {
    # generate bubbleplot
    bubble.plot <- ggplot(barcodes.proportional.melted, aes(x=Position, y=Sample, size = Proportion, color=Color)) +
      geom_point(stat = "identity", alpha = 0.6, shape = 16) +
      scale_color_identity() +
      labs(y = "Condition", x = "Barcode", title = name) +
      scale_size_continuous(range = c(0.1, 10)) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(size=8),
            axis.title = element_text(size= 6))
  }
  print(bubble.plot)
}
