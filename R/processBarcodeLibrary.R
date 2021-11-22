#' processBarcodeLibrary
#'
#' process barcode reference file from raw sequencing datasets
#'
#' @param filepath path to starcode output for barcode library
#' @return returns data frame containing barcode and read count in reference
#' @export
#' @examples
#' process_barcode_reference(filepath = NULL)

process_barcode_reference <- function(filepath = NULL){
  # imports barcode reference file from starcode or similar data frame.
  # Barcode reference should have two columns as below
  # Barcode_ID Count
  # Barcode_1   100

  require(data.table)
  barcodefile <- fread(file = filepath, header = F)
  barcodefile <- tail(barcodefile, -1) # get rid of the header line
  colnames(barcodefile) <- c("ID", "Count")

  # Subset for barcodes over 2 counts
  barcodefile <- barcodefile[barcodefile$Count >= 2,]
  head(barcodefile)
  sum(barcodefile$Count)

  # order by decending barcode count
  ordered <- order(barcodefile$Count, decreasing = T)
  barcodefile <- barcodefile[ordered,]

  # generate Rank and proportion information
  barcodefile.sumreads <- sum(barcodefile$Count)
  barcodefile$Rank <- paste("Barcode_rank_",seq(1,length(barcodefile$ID)), sep = "")

  barcodefile$Proportion <- barcodefile$Count / sum(barcodefile$Count) * 100
  sum(barcodefile$Proportion) # should total 100
  head(barcodefile)

  return(barcodefile)
}


# get object name as string
getname <- function(v1) {
  deparse(substitute(v1))
}

# fasta output from dataframe function
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,1], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,2]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

