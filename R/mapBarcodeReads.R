#' Map Barcode Reads
#'
#' Map reads containing a lintraceR DNA barcode to a barcode reference library/fasta and return CellID:lintracer barcode ID annotations
#'
#' @param reads Dataframe from extractBarcodeReads containing barcode reads
#' @param bowtie_index Path to bowtie index of barcode reference library. Include index prefix.
#' @param mismatches Number of mismatches allowed in bowtie for read mapping
#' @param threads Number of CPU threads for bowtie to use
#' @param outdir Optional. save output to directory
#' @param prefix Optional. prefix to add to output file
#' @param reference_fasta Optional. Fasta file of reference barcode library for alignment. Use instead of bowtie_index.
#' @param cleanup Logical. Remove intermediate files generated by mapBarcodeReads
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence mapped to the lintraceR DNA barcode ID.
#' @export
#' @examples
#' load(system.file("extdata", "test_map.rda", package = "barista"))
#' test.ref <- system.file("extdata", "barcode_lib_reference_test.fasta", package = "barista")
#' mapBarcodeReads(reads = test.map, reference_fasta = test.ref, outdir = tempdir(), prefix = "map")

mapBarcodeReads <- function(reads, bowtie_index = NULL, mismatches = 2, threads = 4,
                            outdir = getwd(), prefix = "map_barcode_reads", cleanup = TRUE,
                            reference_fasta = NULL){

    # welcome message
    message("mapBarcodeReads")

    # check bowtie index exists
    if (is.null(bowtie_index)) {
        if (is.null(reference_fasta)){
            message("ERROR: no bowtie index or reference fasta supplied. Exiting")
        }
    }

    if (is.null(bowtie_index) & !is.null(reference_fasta)) {
        message("WARNING: No bowtie index given. Generating bowtie index from supplied fasta file.")

        # build bowtie index for alignment
        message(paste("Generating bowtie index from file:", reference_fasta))

        index.dir <- file.path(outdir, 'bowtie_index')
        dir.create(path = index.dir, showWarnings = FALSE)
        Rbowtie::bowtie_build(reference_fasta, outdir = index.dir, force = TRUE, prefix = "index")
        bowtie_index <- paste(index.dir, "/index", sep = '')
    }


    # setup temp out
    if (dir.exists(file.path(getwd(), "tmp"))){
        message("tmp exists, deleting")
        unlink(file.path(getwd(), "tmp"), recursive=TRUE)
    }

    temp.dir <- file.path(getwd(), "tmp")
    dir.create(temp.dir, showWarnings = F)

    # create fasta from input dataframe
    if ("Cell.10X.Barcode" %in% colnames(reads)) {
        seqinr::write.fasta(sequences = as.list(reads$DNA.Barcode),
                        names = paste(reads$Read.ID, reads$Cell.10X.Barcode, reads$Read.10X.UMI, sep = "_"),
                        file.out = file.path(temp.dir, "tmp.fa"))
    } else {
        # remove index data from read name which will prevent matching below
        reads$Read.ID <- reshape2::colsplit(reads$Read.ID, ' ', c("Read.ID", "n"))[,1]
        seqinr::write.fasta(sequences = as.list(reads$DNA.Barcode),
                            names = reads$Read.ID,
                            file.out = file.path(temp.dir, "tmp.fa"))
    }

    # map barcode containing reads to reference barcode library
    # run bowtie end-to-end mode (complete alignment no soft-clipping)
    # allow v mismatches in read
    # allowed mismatches and threads can be given as input to function

    Rbowtie::bowtie(sequences = file.path(temp.dir, "tmp.fa"),
                    index = bowtie_index,
                    outfile = file.path(temp.dir,'tmp.sam'),
                    force = TRUE, v = mismatches, f = TRUE,
                    norc = TRUE, t = TRUE, p = threads, sam = TRUE)

    # extract read group mapped for each read
    bam2 <- file.path(temp.dir,'tmp.sam')
    bamFile2 <- Rsamtools::BamFile(Rsamtools::asBam(bam2))
    open(bamFile2)
    index <- Rsamtools::index(bamFile2)
    bam.fields2 <- c("rname", "qname", "seq", "qual", "groupid")
    flags2 <- Rsamtools::scanBamFlag(isDuplicate=FALSE,
                                 isNotPassingQualityControls=FALSE,
                                 isSecondaryAlignment = FALSE)
    parameters2 <- Rsamtools::ScanBamParam(what = bam.fields2, flag = flags2)
    bam.parsed <- Rsamtools::scanBam(bamFile2, param = parameters2)[[1]]
    barcode.match <- data.frame(Read.ID = plyr::ldply(strsplit(bam.parsed$qname,
                            split = "_"))[[1]],
                            referenceID = bam.parsed$rname)
    barcode.match$referenceID <- as.character(barcode.match$referenceID)
    barcode.match$Read.ID <- as.character(barcode.match$Read.ID)

    # merge in reference barcode information
    current.df <- dplyr::left_join(x = reads, y = barcode.match, by = "Read.ID")
    current.df <- stats::na.omit(current.df)

    # keep only unique entries by 10X Cell ID, 10X UMI and reference barcode ID
    if ("Cell.10X.Barcode" %in% colnames(current.df)) {
        current.df <- dplyr::distinct(current.df, current.df$Cell.10X.Barcode, current.df$Read.10X.UMI,
                                      current.df$referenceID, .keep_all = TRUE)
    }

    # generate raw counts
    counts <- barista::countBarcodeReads(bam = bamFile2, bam_index = index, outdir = outdir, prefix = prefix)

    # remove bamfile
    close(bamFile2)

    #print(dim(current.df))
    if (isTRUE(cleanup)) {
        # cleanup tmp files
        files <- list.files(pattern = "tmp", path = temp.dir)
        files <- paste(temp.dir,"/", files, sep = "")
        file.remove(files)
        file.remove(temp.dir)
    }

    # write outfiles
    if (!is.null(outdir)) {
        utils::write.csv(current.df, file = file.path(outdir, paste(prefix,'.csv', sep = '')), quote = F, row.names = F)
        utils::write.csv(counts, file = file.path(outdir, paste(prefix,'_counts.csv', sep = '')), quote = F, row.names = F)
    }

    return(current.df)
}

#' Count Barcode Reads
#'
#' Count number of occurences of each barcode in the reference library
#'
#' @param bam BAM file from mapBarcodeReads containing mapped barcode reads
#' @param bam_index Path to bowtie index of barcode reference library. Include index prefix.
#' @param outdir Optional. Save output to directory
#' @param prefix Optional. Prefix to add to output file
#'
#' @return Returns a data-frame containing the Barcode ID, and the count per sample.
#' @export

countBarcodeReads <- function(bam, bam_index, outdir = ".", prefix = "map_barcode_reads"){
    counts <- Rsamtools::idxstatsBam(file = bam, index = bam_index)
    counts <- counts[,c(1,3)]
    return(counts)
}



