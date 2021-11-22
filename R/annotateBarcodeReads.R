#' Annotate Barcode Reads
#'
#' Matches cells in a single cell experiment to detected DNA barcodes.
#' Dataframe returned will have all cells matched to a barcode
#' If there is no barcode matchable to a cell "not.detected" is returned
#' For cells that have multiple detected barcodes each barcode is returned separated by ';'
#'
#' @param all.cells list of all cells in the experiment to be annotated
#' @param barcoded.cells dataframe a dataframe with cell id and barcode id as columns
#' @return Returns a data-frame containing the the 10X cell ID and the lintraceR DNA barcode ID.
#' @export
#' @examples
#' annotateBarcodeReads(all.cells = cell.list, barcoded.cells = cell.barcode.df)

annotateBarcodeReads <- function(all.cells, cells.w.barcode.df){
  require(dplyr)

  message("Annotating 10X cell ID with DNA barcode ID")
  # make sure this takes the cell id and barcode reference id columns
  cells.w.barcode <- cells.w.barcode.df %>% dplyr::select(Cell.10X.Barcode, referenceID)
  anno.merge <- dplyr::left_join(all.cells,cells.w.barcode,by=c("V1"="Cell.10X.Barcode"))

  # order by cell id and barcode number within cell and get unique entries
  vals <- as.numeric(gsub("[A-Z]*_Barcode_","", anno.merge$referenceID, ignore.case = T))
  anno.order <- anno.merge[order(anno.merge$V1,vals),]
  anno.uniq <- anno.order[!duplicated(anno.order),]

  # collapse multiple barcodes per cell
  anno.collapse <- aggregate(anno.uniq$referenceID, list(anno.uniq$V1), paste, collapse=";")

  # reformat
  anno.collapse$x <- gsub("NA","not.detected",anno.collapse$x)
  colnames(anno.collapse) <- c("cell.id", "barcode")
  anno.final <- data.frame(barcode = anno.collapse$barcode, row.names = gsub('-1','',anno.collapse$cell.id))

  # Annotation summary stats
  print("Total cells annotated:")
  print(nrow(anno.final))

  print("Total cells with barcode:")
  print(nrow(anno.final %>% filter(barcode != "not.detected")))

  print("Total cells without barcode:")
  print(nrow(anno.final %>% filter(barcode == "not.detected")))

  print("Percentage of cells with annotated barcode:")
  print(nrow(anno.final %>% filter(barcode != "not.detected"))/nrow(anno.final)*100)

  print("Number of unique barcodes:")
  print(length(unique(anno.final$barcode)))

  return(anno.final)
}


#-------------------------------------------------------------------------------------



#' Annotate 10X Barcode
#'
#' This function identifies reads containing a lintraceR DNA barcode from the bam alignment file
#' and matches these to the corresponding cell ID in a 10X scRNA-seq experiment. Returns a
#' data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA barcode sequence and
#' lintraceR DNA barcode ID. Inspired by clonehunter.
#'
#' @param bamfile bam file containing barcoded reads
#' @param barcode_pattern pattern to search for in the reads
#' @param reference_library fasta file of barcode reference library to match against
#' @param constant the constant region flanking the barcode, two constant regions can be input as a
#' character vector
#' @param yieldSize how many lines of BAM file to take at each pass of the loop
#' @param threads Number of CPU threads
#' @param bowtie_mismatches number of allowed mismatches in bowtie
#' @param draw_plots logical. do you want plots to be generated?
#' @param cleanup remove intermediate files produced from each BAM chunk
#' @param outfile optional. save output to file
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence and lintraceR DNA barcode ID.
#' @export
#' @examples
#' annotate_10X_barcode(bam = bamfile, pattern = "ACTG", reference_library = ref)


annotate_10X_barcode <- function(bamfile, barcode_pattern, reference_library, constant = "",
                                 yieldSize = 1000000, threads = 4, bowtie_mismatches = 1,
                                 draw_plots = FALSE, cleanup = TRUE, outfile = NULL) {

  # confirm input is bam file and open for reading
  if (endsWith(bamfile, "bam")) {
    bam <- bamfile
    bamFile <- Rsamtools::BamFile(bam)
    Rsamtools::yieldSize(bamFile) <- 1000000
    open(bamFile)
    bam.fields <- c("rname", "qname", "seq", "qual", "groupid")
    flags <- Rsamtools::scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
    parameters <- Rsamtools::ScanBamParam(what = bam.fields, tag = c("CB", "UB"), flag = flags)
    bam.parsed.df <- data.frame()
  } else {
    stop("ERROR: input file must be a BAM file")
  }

  # if input is fasta exit. TODO
  if (endsWith(reference_library, "fasta")) {
    file.exists(reference_library)
  } else {
    stop("ERROR: reference library file must be a fasta file")
  }

  # setup counter to report iteration
  start = 0
  end = Rsamtools::yieldSize(bamFile)
  count = 0

  # setup temp out
  temp.dir <- file.path(getwd(), "bowtie")
  dir.create(temp.dir, showWarnings = F)

  # check if bowtie index exists
  index.dir <- file.path(temp.dir, 'bowtie_index')
  dir.create(path = index.dir, showWarnings = FALSE)

  if (file.exists(file.path(index.dir,"/index.1.ebwt"))) {
    message("bowtie index already exists")
  } else {
    # build bowtie index for alignment
    message("Generating bowtie index")
    Rbowtie::bowtie_build(reference_library, outdir = index.dir, force = TRUE, prefix = 'index')
  }

  # iterate over BAM file by yieldSize
  while(TRUE){
    current.reads <- Rsamtools::scanBam(bamFile, param = parameters)[[1]]
    message("")
    update <- paste("Parsing reads", start, "to", end)
    message(update)
    if (length(current.reads$qname) <= 0){
      message("No reads contained a 10X cell ID. Skip")
      break
    } else {
      # select reads
      current.seqs <- as.character(current.reads$seq)

      # identify matches
      matches <- grep(barcode_pattern, current.seqs, perl = TRUE)

      # keep the full sequences to check
      full.seqs <- current.seqs[matches]
      seqs <- stringr::str_extract(current.seqs, barcode_pattern)
      seqs <- seqs[!is.na(seqs)]

      # check lengths of index and seq are correct
      length(matches) == length(seqs)

      # report counts
      print(paste(length(matches), "reads match barcode pattern"))

      if (length(matches) > 0) {

        # get read ID, UMI and cell ID of match
        match.qname <- current.reads$qname[matches]
        match.quals <- current.reads$qual[matches]
        match.UB <- current.reads$tag$UB[matches]
        match.CB <- current.reads$tag$CB[matches]

        if (!(is.null(match.UB) | is.null(match.CB))) {
          match.seqs <- current.seqs[matches]

          # extract the barcode sequence between constant regions for matches
          DNA.barcodes <- str_extract(match.seqs, barcode_pattern)

          # initialize the current data table
          current.df <- data.table::data.table(Read.ID = match.qname,
                                               Read.Qual = as.character(match.quals),
                                               Cell.10X.Barcode = match.CB,
                                               Read.10X.UMI = match.UB,
                                               DNA.Barcode = DNA.barcodes)

          # filter rows that contain NAs
          NA.reads <- length(which(is.na(current.df$Cell.10X.Barcode)))
          total.reads <- length(current.df$Cell.10X.Barcode)
          filtered.reads <- length(current.df$Cell.10X.Barcode)
          print(paste(filtered.reads, " of ", total.reads,
                      "reads (", floor(100*(filtered.reads/total.reads)),
                      "% ) have a 10X cell ID match"), sep = "")

          current.df <- na.omit(current.df)

          # extract only the barcode region
          constant <- as.vector(constant)
          if (length(constant) == 1){
            current.df$DNA.Barcode <- gsub(pattern = constant[1], replacement = "",
                                           x = current.df$DNA.Barcode)
          }
          if (length(constant) > 1){
            current.df$DNA.Barcode <- gsub(pattern = constant[1], replacement = "",
                                           x = current.df$DNA.Barcode)
            current.df$DNA.Barcode <- gsub(pattern = constant[2], replacement = "",
                                           x = current.df$DNA.Barcode)
          }
          if (length(constant) > 2){
            message("ERROR: constant can have up to two elements only")
            break
          }

          # match cell DNA barcode to reference library
          seqinr::write.fasta(sequences = as.list(current.df$DNA.Barcode),
                              names = paste(current.df$Read.ID, current.df$Cell.10X.Barcode, sep = "_"),
                              file.out = file.path(temp.dir, paste("tmp_reads_",count,".fasta", sep = "")))
          reads <- file.path(temp.dir,paste("tmp_reads_",count,".fasta", sep = ""))

          # run bowtie end-to-end mode (complete alignment no soft-clipping)
          # allow v mismatches in read
          # allowed mismatches and threads can be given as input to function
          Rbowtie::bowtie(sequences = reads,
                          index = paste(index.dir,'/index', sep = ''),
                          outfile = file.path(temp.dir, paste("/tmp_",count,'.sam', sep = "")),
                          force = TRUE, v = bowtie_mismatches, f = TRUE,
                          norc = TRUE, t = TRUE, p = threads, sam = TRUE)

          # extract read group mapped for each read
          bam2 <- file.path(temp.dir, paste("/tmp_", as.character(count), '.sam', sep = ""))
          bamFile2 <- Rsamtools::BamFile(asBam(bam2))
          open(bamFile2)
          bam.fields2 <- c("rname", "qname", "seq", "qual", "groupid")
          flags2 <- Rsamtools::scanBamFlag(isDuplicate=FALSE,
                                           isNotPassingQualityControls=FALSE,
                                           #isUnmappedQuery = FALSE,
                                           isSecondaryAlignment = FALSE)
          parameters2 <- Rsamtools::ScanBamParam(what = bam.fields2, flag = flags2)
          bam.parsed <- Rsamtools::scanBam(bamFile2, param = parameters2)[[1]]
          barcode.match <- data.frame(Read.ID = ldply(strsplit(bam.parsed$qname,
                                                               split = "_"))[[1]],
                                      #Cell_10X_ID = ldply(strsplit(bam.parsed$qname,
                                      #split = "_"))[[2]],
                                      referenceID = bam.parsed$rname)

          # merge in reference barcode information
          current.df <- dplyr::left_join(x = current.df, y = barcode.match,
                                         by = "Read.ID")
          current.df <- na.omit(current.df)

          # keep only unique entries by 10X Cell ID, 10X UMI and reference barcode ID
          current.df <- dplyr::distinct(current.df, Cell.10X.Barcode, Read.10X.UMI,
                                        referenceID, .keep_all = TRUE)


          if (isTRUE(cleanup)) {
            # cleanup tmp files
            files <- list.files(pattern = "tmp", path = paste(getwd(),'/bowtie/',
                                                              sep = ''))
            files <- paste(getwd(),'/bowtie/', files, sep = "")
            file.remove(files)
          }
          # Merge current set into the final data frame
          if (nrow(bam.parsed.df) <= 0) {
            bam.parsed.df <- current.df
          } else {
            bam.parsed.df <- rbind(bam.parsed.df, current.df)
          }
        } else {
          next
        }


      }


    }
    count <- count + 1
    start <- Rsamtools::yieldSize(bamFile) * count
    end <- start + end + 1


  }


  if (isTRUE(draw_plots)) {
    # cell ID frequency
    dir.create(file.path(getwd(),"plots"))
    plots.dir <- file.path(getwd(),"plots")

    png(paste(plots.dir, "/barcodes_per_10X_cellID.png", sep = ""))
    hist(table(bam.parsed.df$Cell.10X.Barcode),
         col = "steelblue",
         breaks = 50,
         main = paste("Distribution of 10X cell ID tags in BAM entries", start, "to", end))
    dev.off()

    png(paste(plots.dir, "/barcode_match_lengths.png", sep = ""))
    lengths <- sapply(X = bam.parsed.df$DNA.Barcode, FUN = nchar, simplify = T)
    hist(lengths,
         col = "steelblue",
         breaks = 20,
         main = paste("Distribution of barcode match lengths in BAM entries", start, "to", end))
    dev.off()

  }

  close(bamFile)
  if (!is.null(outfile)) {
    write.csv(bam.parsed.df, file = paste(getwd(), outfile,'.csv', sep = ""))
  }

  return(bam.parsed.df)
}

