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
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence mapped to the lintraceR DNA barcode ID.
#' @export
#' @examples
#' mapBarcodeReads(reads = file.bam, reference_library = ref.fasta, mismatches = 0, outdir = "./", prefix = "map_barcode_reads")

mapBarcodeReads <- function(reads, bowtie_index = NULL, mismatches = 2, threads = 4,
                            outdir = getwd(), prefix = "map_barcode_reads", cleanup = TRUE,
                            reference_fasta = NULL){

    # welcome message
    base::message("mapBarcodeReads")

    # check bowtie index exists
    if (base::is.null(bowtie_index)) {
        if (base::is.null(reference_fasta)){
            base::message("ERROR: no bowtie index or reference fasta supplied. Exiting")
        }
    }

    if (base::is.null(bowtie_index) & !base::is.null(reference_fasta)) {
        base::message("WARNING: No bowtie index given. Generating bowtie index from supplied fasta file.")

        # build bowtie index for alignment
        base::message(base::paste("Generating bowtie index from file:", reference_fasta))

        index.dir <- base::file.path(outdir, 'bowtie_index')
        base::dir.create(path = index.dir, showWarnings = FALSE)
        Rbowtie::bowtie_build(reference_fasta, outdir = index.dir, force = TRUE, prefix = "index")
        bowtie_index <- base::paste(index.dir, "/index", sep = '')
    }


    # setup temp out
    if (base::dir.exists(file.path(getwd(), "tmp"))){
        base::message("tmp exists, deleting")
        base::unlink(file.path(getwd(), "tmp"), recursive=TRUE)
    }

    temp.dir <- base::file.path(getwd(), "tmp")
    base::dir.create(temp.dir, showWarnings = F)

    # create fasta from input dataframe
    if ("Cell.10X.Barcode" %in% base::colnames(reads)) {
        seqinr::write.fasta(sequences = base::as.list(reads$DNA.Barcode),
                        names = base::paste(reads$Read.ID, reads$Cell.10X.Barcode, reads$Read.10X.UMI, sep = "_"),
                        file.out = base::file.path(temp.dir, "tmp.fa"))
    } else {
        # remove index data from read name which will prevent matching below
        reads$Read.ID <- reshape2::colsplit(reads$Read.ID, ' ', c("Read.ID", "n"))[,1]
        seqinr::write.fasta(sequences = as.list(reads$DNA.Barcode),
                            names = reads$Read.ID,
                            file.out = base::file.path(temp.dir, "tmp.fa"))
    }

    # map barcode containing reads to reference barcode library
    # run bowtie end-to-end mode (complete alignment no soft-clipping)
    # allow v mismatches in read
    # allowed mismatches and threads can be given as input to function

    Rbowtie::bowtie(sequences = base::file.path(temp.dir, "tmp.fa"),
                    index = bowtie_index,
                    outfile = base::file.path(temp.dir,'tmp.sam'),
                    force = TRUE, v = mismatches, f = TRUE,
                    norc = TRUE, t = TRUE, p = threads, sam = TRUE)

    # extract read group mapped for each read
    bam2 <- base::file.path(temp.dir,'tmp.sam')
    bamFile2 <- Rsamtools::BamFile(Rsamtools::asBam(bam2))
    base::open(bamFile2)
    index <- Rsamtools::index(bamFile2)
    bam.fields2 <- c("rname", "qname", "seq", "qual", "groupid")
    flags2 <- Rsamtools::scanBamFlag(isDuplicate=FALSE,
                                 isNotPassingQualityControls=FALSE,
                                 isSecondaryAlignment = FALSE)
    parameters2 <- Rsamtools::ScanBamParam(what = bam.fields2, flag = flags2)
    bam.parsed <- Rsamtools::scanBam(bamFile2, param = parameters2)[[1]]
    barcode.match <- base::data.frame(Read.ID = plyr::ldply(base::strsplit(bam.parsed$qname,
                            split = "_"))[[1]],
                            referenceID = bam.parsed$rname)
    barcode.match$referenceID <- base::as.character(barcode.match$referenceID)
    barcode.match$Read.ID <- base::as.character(barcode.match$Read.ID)

    # merge in reference barcode information
    current.df <- dplyr::left_join(x = reads, y = barcode.match, by = "Read.ID")
    current.df <- stats::na.omit(current.df)

    # keep only unique entries by 10X Cell ID, 10X UMI and reference barcode ID
    if ("Cell.10X.Barcode" %in% base::colnames(current.df)) {
        current.df <- dplyr::distinct(current.df, Cell.10X.Barcode, Read.10X.UMI,
                                      referenceID, .keep_all = TRUE)
    }

    # generate raw counts
    counts <- barcodeR::countBarcodeReads(bam = bamFile2, bam_index = index, outdir = outdir, prefix = prefix)

    # remove bamfile
    base::close(bamFile2)

    #print(dim(current.df))
    if (base::isTRUE(cleanup)) {
        # cleanup tmp files
        files <- base::list.files(pattern = "tmp", path = temp.dir)
        files <- base::paste(temp.dir,"/", files, sep = "")
        base::file.remove(files)
        base::file.remove(temp.dir)
    }

    # write outfiles
    if (!base::is.null(outdir)) {
        utils::write.csv(current.df, file = base::file.path(outdir, base::paste(prefix,'.csv', sep = '')), quote = F, row.names = F)
        utils::write.csv(counts, file = base::file.path(outdir, base::paste(prefix,'_counts.csv', sep = '')), quote = F, row.names = F)
    }

    return(current.df)
}


#-------------------------------------------------------------------------------------

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
#' @examples
#' countBarcodeReads(bam = bamfile, reference_library, mismatches = 1, outdir = "./", prefix = "map_barcode_reads")

countBarcodeReads <- function(bam, bam_index, outdir = ".", prefix = "map_barcode_reads"){
    counts <- Rsamtools::idxstatsBam(file = bam, index = bam_index)
    counts <- counts[,c(1,3)]
    return(counts)
}


#' Align and Count DNA Barcodes
#'
#' Align lintraceR DNA barcodes to reference barcode library and summarise counts
#'
#' @param reads REQUIRED. Path to fasta file from extract_DNA_barcodes containing barcode reads to align.
#' @param reference REQUIRED. Path to fasta file containing barcode reference library to align against.
#' @param outfile REQUIRED. path of output SAM file. Will be converted to BAM. Counts file will be written as <prefix>.counts.
#' @param bowtie_mismatches Numeric. Allowe number of mismatches between read and reference index. Default = 1.
#' @param overwrite Boolean. Overwrite output files? Default = FALSE.
#' @param threads Numeric. Number of CPU threads to use. Default = 2.
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence and lintraceR DNA barcode ID.
#' @export
#' @examples
#' align_count_DNA_barcodes(reads = file.fa, reference_library = ref, outfile = out.bam, bowtie_mismatches = 1, overwrite = T)

align_count_DNA_barcodes <- function(reads, reference, outfile = NULL, bowtie_mismatches = 1, overwrite = FALSE, threads = 2){

    # confirm input fasta file
    if (endsWith(reads, ".fasta") | endsWith(reads, ".fa")) {
        message(paste("Aligning barcode reads from", reads))
    } else {
        stop("ERROR: reads is not a valid fasta file. Exiting.")
    }

    # confirm reference library is fasta
    if (endsWith(reference, "fasta")) {
        file.exists(reference)
    } else {
        stop("ERROR: reference library file must be a fasta file")
    }

    # check whether outfile already exists
    if (is.null(outfile)){
        stop("ERROR: No output file given.")
    }

    # generate bowtie reference index
    message("Generating Bowtie index")
    index.dir <- file.path(getwd(), 'bowtie_index')
    dir.create(path = index.dir, showWarnings = FALSE)
    if (length(list.files(index.dir, all.files = TRUE, no.. = TRUE, include.dirs = TRUE)) == 0) {
        Rbowtie::bowtie_build(reference, outdir = index.dir, force = TRUE)
    } else {
        message("Bowtie index exists. Skipping")
    }

    # setup checks for output bam file and index
    if (file.exists(outfile) & isTRUE(overwrite)){
        message(paste("WARNING:", outfile, "exists and overwrite is set to TRUE.", outfile, "will be overwritten."))
        file.remove(outfile)
    }

    if(file.exists(outfile) && isFALSE(overwrite)){
        stop("ERROR: Outfile exists and overwrite is set to FALSE. Exiting.")
    }

    if (file.exists(paste(tools::file_path_sans_ext(outfile), ".bai", sep = ""))){
        file.remove(paste(tools::file_path_sans_ext(outfile), ".bai", sep = ""))
    }

    # map barcode containing reads to reference barcode library
    # run bowtie end-to-end mode (complete alignment no soft-clipping)
    # allow v mismatches in read
    # allowed mismatches and threads can be given as input to function
    out.tmp <- paste(tools::file_path_sans_ext(outfile), ".sam", sep = "")
    Rbowtie::bowtie(sequences = reads,
                    index = paste(index.dir,"/index", sep = ""),
                    outfile = out.tmp,
                    force = TRUE, v = bowtie_mismatches, f = TRUE,
                    norc = TRUE, t = TRUE, p = threads, sam = TRUE)
    bam <- BamFile(asBam(out.tmp))
    index <- Rsamtools::index(bam)
    file.remove(out.tmp)

    counts <- Rsamtools::idxstatsBam(file = bam, index = index)
    counts <- counts[,c(1,3)]
    head(counts)
    write.table(counts, file = paste(tools::file_path_sans_ext(outfile), ".counts", sep = ""), quote = F, sep = "\t", append = FALSE, row.names = FALSE)
}




