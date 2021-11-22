#' extractBarcodeReads
#'
#' Identify and extract reads containing DNA barcodes from a FASTQ or BAM alignment file
#'
#' @param infile Input fastq or bam file containing barcoded reads
#' @param barcode_pattern Pattern to search for in the reads
#' @param constant The constant region flanking the barcode, two constant regions can be input as a
#' character vector. This should be set if the constant region is part of the barcode_pattern.
#' Up to two constant regions can be set e.g. c("ACTG","TACG")
#' @param yieldSize How many lines of input file to stream at each pass
#' @param outfile Optional. path to save parsed reads as output file
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence and lintraceR DNA barcode ID.
#' @export
#' @examples
#' extractBarcodeReads(infile = file.fastq, barcode_pattern = "ACTG")
#' extractBarcodeReads(infile = file.bam, barcode_pattern = "ACTG")

extractBarcodeReads <- function(infile,
                                barcode_pattern,
                                constant,
                                yieldSize = 1e6,
                                outfile = NULL,
                                clean = TRUE,
                                verbose = FALSE,
                                overwrite = TRUE){

    # welcome message
    base::message("extractBarcodeReads")

    # confirm if input is fastq or bam file and make call to appropriate function
    if (endsWith(infile, ".bam")) {
        message(paste("Input file: ", infile))
        message("Format is BAM.")
        out <- barcodeR::extractBAM(infile, barcode_pattern = barcode_pattern,
                                    constant = constant, yieldSize = yieldSize, outfile = outfile)
    }

    if (endsWith(infile, ".fastq") | endsWith(infile, ".fq") |
        endsWith(infile, ".fastq.gz") | endsWith(infile, ".fq.gz")) {

        message(paste("Input file: ", infile))
        message("Format is FASTQ.")

        out <- barcodeR::extractFASTQ(fastq = infile, barcode_pattern = barcode_pattern,
                                      constant = constant, outfile = outfile, block_size = yieldSize,
                                      verbose = verbose, overwrite = overwrite)

    }
    return(out)
}


#' extractBAM
#'
#' Barcode extraction function for BAM alignment files
#'
#' @param bam BAM file containing barcoded reads
#' @param barcode_pattern Pattern to search for in the reads
#' @param constant The constant region flanking the barcode, two constant regions can be input as a
#' character vector. This should be set if the constant region is part of the barcode_pattern.
#' Up to two constant regions can be set e.g. c("ACTG","TACG")
#' @param yieldSize How many lines of input file to stream at each pass
#' @param outfile Optional. path to save parsed reads as output file
#'
#' @return Returns a data-frame containing the read ID, the 10X cell ID, 10X umi, lintraceR DNA
#' barcode sequence and lintraceR DNA barcode ID.
#' @export
#' @examples
#' extractBAM(infile = fastq/bam, barcode_pattern = "ACTG")

extractBAM <- function(bam, barcode_pattern, constant = "",
                       yieldSize = 1e6, outfile = NULL) {
    # confirm input is bam file and open for reading
    if (base::endsWith(bam, ".bam")) {
        bamFile <- Rsamtools::BamFile(bam)
        Rsamtools::yieldSize(bamFile) <- yieldSize
        base::open(bamFile)
        bam.fields <- c("rname", "qname", "seq", "qual", "groupid")
        flags <- Rsamtools::scanBamFlag(isDuplicate = FALSE, isNotPassingQualityControls = FALSE)

        # these flags are the cell IDs and UMIs
        parameters <- Rsamtools::ScanBamParam(what = bam.fields, tag = c("CB", "UB"), flag = flags)

    } else {
        stop("ERROR: input file must be a BAM file")
    }

    # setup counter to report iteration
    start = 0
    end = Rsamtools::yieldSize(bamFile)
    count = 1

    # setup container to hold all matching reads
    bam.parsed.df <- base::data.frame()

    # iterate over BAM file by yieldSize
    while(TRUE){
        current.reads <- Rsamtools::scanBam(bamFile, param = parameters)[[1]]
        base::message("")
        update <- base::paste("Parsing reads", as.numeric(start), "to", as.numeric(end))
        message(update)
        if (base::length(current.reads$qname) <= 0){
            base::message("No reads contained a 10X cell ID. Skip")
            break
        } else {
            # select reads
            current.seqs <- base::as.character(current.reads$seq)

            # identify matches
            matches <- base::grep(barcode_pattern, current.seqs, perl = FALSE)
            pc.matches <- 100*(base::length(matches) / base::length(current.seqs))
            update <- base::paste(pc.matches, "% of reads contain the barcode pattern", sep = '')
            base::message(update)

            # keep the full sequences to check
            full.seqs <- current.seqs[matches]
            seqs <- stringr::str_extract(current.seqs, barcode_pattern)
            seqs <- seqs[!is.na(seqs)]

            # check lengths of index and seq are correct
            base::length(matches) == base::length(seqs)

            # process only reads that match the barcode pattern
            if (base::length(matches) > 0) {

                # get read ID, UMI and cell ID of match
                match.qname <- current.reads$qname[matches]
                match.quals <- current.reads$qual[matches]
                match.UB <- current.reads$tag$UB[matches]
                match.CB <- current.reads$tag$CB[matches]

                # skip if there are no cell barcodes or UMIs detected
                if (!(base::is.null(match.UB) | base::is.null(match.CB))) {
                    match.seqs <- current.seqs[matches]

                    # extract the barcode sequence between constant regions for matches
                    DNA.barcodes <- stringr::str_extract(match.seqs, barcode_pattern)

                    # initialize the current data table
                    current.df <- data.table::data.table(Read.ID = match.qname,
                                                         Read.Qual = as.character(match.quals),
                                                         Cell.10X.Barcode = match.CB,
                                                         Read.10X.UMI = match.UB,
                                                         DNA.Barcode = DNA.barcodes)

                    # filter rows that contain NAs
                    total.reads <- base::length(current.df$Cell.10X.Barcode)
                    current.df <- stats::na.omit(current.df)
                    filtered.reads <- base::length(current.df$Cell.10X.Barcode)

                    update <- base::paste(base:floor(100 * (filtered.reads/total.reads)),
                                          "% of reads have a 10X cell ID match", sep = "")
                    base::message(update)

                    # barcode currently also contains the constant region, extract only the variable barcode
                    constant <- base::as.vector(constant)
                    if (base::length(constant) == 1){
                        current.df$DNA.Barcode <- base::gsub(pattern = constant[1],
                                                             replacement = "",
                                                             x = current.df$DNA.Barcode)
                    }
                    if (base::length(constant) > 1){
                        current.df$DNA.Barcode <- base::gsub(pattern = constant[1],
                                                       replacement = "",
                                                       x = current.df$DNA.Barcode)
                        current.df$DNA.Barcode <- base::gsub(pattern = constant[2],
                                                       replacement = "",
                                                       x = current.df$DNA.Barcode)
                    }
                    if (base::length(constant) > 2){
                        base::message("ERROR: constant can have up to two elements only")
                        break
                    }

                    # add good reads to output data.frame
                    bam.parsed.df <- base::rbind(current.df, bam.parsed.df)
                }
            }

            # update counters
            count <- count + 1
            start <- start + yieldSize
            end <- start + yieldSize

        }
    }
    # report number of barcode reads
    base::message(base::paste("Identified", base::nrow(bam.parsed.df), "reads that contained the barcode pattern"))

    # close Bamfile and return results
    base::close(bamFile)

    # if outdir write resulting barcode reads to file
    if(!base::is.null(outfile)){
        if(base::endsWith(outfile, ".fa") | base::endsWith(outfile, ".fasta")){
            seqinr::write.fasta(sequences = base::as.list(bam.parsed.df$DNA.Barcode),
                                names = base::paste(bam.parsed.df$Read.ID, bam.parsed.df$Cell.10X.Barcode, sep = "_"),
                                file.out = base::file.path(outfile))
        } else {
            seqinr::write.fasta(sequences = base::as.list(bam.parsed.df$DNA.Barcode),
                                names = base::paste(bam.parsed.df$Read.ID, bam.parsed.df$Cell.10X.Barcode, sep = "_"),
                                file.out = base::file.path(paste(outfile, '.fa', sep = '')))

        }

    }
    return(bam.parsed.df)
}


#' extractFASTQ
#'
#' Barcode extraction function for raw FASTQ files
#'
#' @param fastq REQUIRED. Path to input fastq file
#' @param barcode_pattern REQUIRED. Pattern to search for in the reads.
#' @param constant The constant region(s) flanking the barcode. This should be set if the
#' constant region is part of the barcode_pattern.
#' Up to two constant regions can be set e.g. constant = c("ACTG","TACG")
#' @param outfile OPTIONAL. path to save filtered reads as output fasta file with read name as sequence ID
#' @param block_size Number of reads to load into memory per chunk. Default = 5e5
#' @param verbose Boolean. Show verbose output from ShortRead? Default = FALSE.
#' @param overwrite Boolean. Overwrite output files? Default = FALSE.
#'
#' @return Returns a data-frame containing the read ID and the extracted DNA
#' barcode sequence
#' @export
#' @examples
#' extractFASTQ(fastq = file.fq.gz, barcode_pattern = "ACTG", constant = NULL, outfile = NULL, block_size = 1e6, verbose = FALSE, overwrite = FALSE)

extractFASTQ <- function(fastq, barcode_pattern, constant = NULL, outfile = NULL, block_size = 1e6, verbose = FALSE, overwrite = FALSE){
    # confirm input fastq file
    #print(fastq)
    if (base::endsWith(fastq, ".fastq") | base::endsWith(fastq, ".fq") |
        base::endsWith(fastq, ".fastq.gz") | base::endsWith(fastq, ".fq.gz")) {
        base::message("Extracting barcode reads.")

        # open input stream
        fq.stream <- ShortRead::FastqStreamer(fastq, verbose = verbose, n = block_size)
        base::on.exit(base::close(fq.stream))

    } else {
        stop("ERROR: input file is not a valid fastq file. Exiting.")
    }

    # confirm barcode pattern is set
    if(base::is.null(barcode_pattern)){
        stop("ERROR: Barcode pattern has not been set. Exiting.")
    }

    # check whether outfile already exists
    if (!base::is.null(outfile)){
        if(base::file.exists(outfile) && base::isFALSE(overwrite)){
            stop("ERROR: Outfile exists and overwrite is set to FALSE. Exiting.")
        }
        if (base::file.exists(outfile) && base::isTRUE(overwrite)){
            base::message(paste("WARNING:", outfile, "exists and overwrite is set to TRUE.", outfile, "will be overwritten."))
            base::file.remove(outfile)
        }
    }
    # setup container to hold all matching reads
    fastq.parsed.df <- base::data.frame()

    # setup counters
    total.reads <- 0
    total.reads.barcode <- 0

    # setup output list
    output <- base::vector(mode = "list", length = 100)
    iter <- 1

    # stream fastq file into memory
    repeat {
        fq <- ShortRead::yield(fq.stream)
        if (base::length(fq) == 0)
            break

        # increment
        iter <- iter + 1

        # determine read length
        avg.read.len <- base::mean(BiocGenerics::width(fq@sread))
        reads.w.barcode <- base::grep(fq@sread, pattern = barcode_pattern)

        # count total reads matching criteria and update global counts
        total.reads <- total.reads + base::length(fq@sread)
        total.reads.barcode <- total.reads.barcode + base::length(reads.w.barcode)

        # subset barcode containing reads
        barcode.reads <- fq[grep(fq@sread, pattern = barcode_pattern)]

        #setup output dataframe
        df <- base::data.frame(id = base::as.data.frame(barcode.reads@id)$x,
                         seq = base::lapply(base::as.data.frame(barcode.reads@sread),
                                      FUN = stringr::str_extract, pattern = barcode_pattern))
        # bind into list
        output[[iter]] <- df

    }

    # calculate performance metrics
    percent.barcode <- 100*(total.reads.barcode / total.reads)

    # combine rows into single data frame
    fastq.parsed.df <- dplyr::bind_rows(output)

    # fix colnames
    base::colnames(fastq.parsed.df) <- c("Read.ID", "DNA.Barcode")

    # remove constant regions
    if (!base::is.null(constant))
        if(base::length(constant) > 2){
            stop("ERROR: A maximum of two constant regions are allowed. Exiting.")
        }
    if(base::length(constant) == 1){
        fastq.parsed.df$DNA.Barcode <- base::gsub(pattern = constant, replacement = '', fastq.parsed.df$DNA.Barcode)
    }
    if(base::length(constant) == 2){
        tmp <- base::gsub(pattern = constant[1], replacement = '', fastq.parsed.df$DNA.Barcode)
        fastq.parsed.df$DNA.Barcode <- base::gsub(pattern = constant[2], replacement = '', tmp)
    }

    # write entire data frame to fasta file.
    if (!base::is.null(outfile)){
        seqinr::write.fasta(sequences = base::as.list(fastq.parsed.df$DNA.Barcode), names = fastq.parsed.df$Read.ID,
                            file.out = base::file.path(outfile), as.string = TRUE)
        base::message(paste("Percentage of reads with barcode: ", percent.barcode))
        base::message(paste("Parsing complete for:", fastq))
    } else {
        base::message(paste("Percentage of reads with barcode: ", percent.barcode))
        base::message(paste("Parsing complete for:", fastq))

    }
    return(fastq.parsed.df)
}






#-------------------------------------------------------------
#' Extract DNA Barcode Reads
#'
#' Extract lintraceR DNA barcodes from raw fastq files
#'
#' @param read1 REQUIRED. Path to fastq file read1
#' @param barcode_pattern REQUIRED. Pattern to search for in the reads.
#' @param constant The constant region(s) flanking the barcode. This should be set if the
#' constant region is part of the barcode_pattern.
#' Up to two constant regions can be set e.g. constant = c("ACTG","TACG")
#' @param outfile OPTIONAL. path to save filtered reads as output fasta file with read name as sequence ID
#' @param block_size Number of reads to load into memory per chunk. Default = 5e5
#' @param verbose Boolean. Show verbose output from ShortRead? Default = FALSE.
#' @param overwrite Boolean. Overwrite output files? Default = FALSE.
#'
#' @return Returns a data-frame containing the read ID and the extracted DNA
#' barcode sequence
#' @export
#' @examples
#' extract_DNA_barcodes(read1 = file.fq.gz, barcode_pattern = "ACTG", constant = NULL, outfile = NULL, block_size = 5e5, verbose = FALSE, overwrite = FALSE)

extract_DNA_barcodes <- function(read1, barcode_pattern, constant = NULL, outfile = NULL, block_size = 2e6, verbose = FALSE, overwrite = FALSE){
    # confirm input fastq file
    if (endsWith(read1, ".fastq") | endsWith(read1, ".fq") | endsWith(read1, ".fastq.gz") | endsWith(read1, ".fq.gz")) {
        message(paste("Extracting barcode reads from", read1))
        # open input stream
        fq.stream <- FastqStreamer(read1, verbose = verbose, n = block_size)
        on.exit(close(fq.stream))

    } else {
        stop("ERROR: read1 is not a valid fastq file. Exiting.")
    }

    # confirm barcode pattern is set
    if(is.null(barcode_pattern)){
        stop("ERROR: Barcode pattern has not been set. Exiting.")
    }

    # check whether outfile already exists
    if (!is.null(outfile)){
        if(file.exists(outfile) && isFALSE(overwrite)){
            stop("ERROR: Outfile exists and overwrite is set to FALSE. Exiting.")
        }
        if (file.exists(outfile) && isTRUE(overwrite)){
            message(paste("WARNING:", outfile, "exists and overwrite is set to TRUE.", outfile, "will be overwritten."))
            file.remove(outfile)
        }
    }
    # setup container to hold all matching reads
    fastq.parsed.df <- data.frame()

    # setup counters
    total.reads <- 0
    total.reads.barcode <- 0

    # setup output list
    output <- vector(mode = "list", length = 100)
    iter <- 1
    # stream fastq file into memory
    repeat {
        fq <- yield(fq.stream)
        if (length(fq) == 0)
            break

        # increment
        iter <- iter + 1

        # determine read length
        avg.read.len <- mean(width(fq@sread))
        reads.w.barcode <- grep(fq@sread, pattern = barcode_pattern)

        # count total reads matching criteria and update global counts
        total.reads <- total.reads + length(fq@sread)
        total.reads.barcode <- total.reads.barcode + length(reads.w.barcode)

        # subset barcode containing reads
        barcode.reads <- fq[grep(fq@sread, pattern = barcode_pattern)]

        #setup output dataframe
        df <- data.frame(id = as.data.frame(barcode.reads@id)$x,
                         seq = lapply(as.data.frame(barcode.reads@sread),
                                      FUN = str_extract, pattern = barcode_pattern))
        # bind into list
        output[[iter]] <- df

    }

    # calculate performance metrics
    percent.barcode <- 100*(total.reads.barcode / total.reads)

    # combine rows into single data frame
    fastq.parsed.df <- dplyr::bind_rows(output)

    # fix colnames
    colnames(fastq.parsed.df) <- c("id", "seq")

    # remove constant regions
    if (!is.null(constant))
        if(length(constant) > 2){
            stop("ERROR: A maximum of two constant regions are allowed. Exiting.")
        }
    if(length(constant) == 1){
        fastq.parsed.df$seq <- gsub(pattern = constant, replacement = '', fastq.parsed.df$seq)
    }
    if(length(constant) == 2){
        tmp <- gsub(pattern = constant[1], replacement = '', fastq.parsed.df$seq)
        fastq.parsed.df$seq <- gsub(pattern = constant[2], replacement = '', tmp)
    }

    # write entire data frame to file.
    if (!is.null(outfile)){
        seqinr::write.fasta(sequences = as.list(fastq.parsed.df$seq), names = fastq.parsed.df$id, file.out = file.path(outfile), as.string = TRUE)
        message(paste("Percentage of reads with barcode: ", percent.barcode))
        message(paste("Parsing complete for:", read1))
    } else {
        message(paste("Percentage of reads with barcode: ", percent.barcode))
        message(paste("Parsing complete for:", read1))
        return(fastq.parsed.df)
    }
}

