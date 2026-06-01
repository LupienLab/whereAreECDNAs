#' Perform peak calling with MACS2
#'
#' @description
#' Performs peak calling with MACS2 on a three-column BED file derived from
#' pseudobulk fragments. The resulting `.narrowPeak` file is imported as a
#' `GRanges` object and stored in the `narrowPeak` slot.
#'
#' @details
#' The default MACS2 settings are chosen for pseudobulk single-cell ATAC-seq
#' fragments represented as BED-like intervals. `--nomodel`, `--shift -100`,
#' and `--extsize 200` are commonly used to center and extend Tn5 insertion
#' signal when calling narrow ATAC-seq peaks from BED input. For more stringent
#' peak sets, consider reducing `qval` from `0.05` to `0.01`.
#'
#' @param object A `whereAreECDNAs` object.
#' @param macs2Path Character(1). Path to `macs2` executable.
#' @param outputDir Character(1). Base output directory. If `NULL`, uses the
#'   `outputDirectory` stored in `object`, or `paste0(sampleName, ".outputs")`
#'   when no stored directory is available.
#' @param format Character(1). MACS2 input format. The current implementation
#'   writes three-column BED input and therefore supports `"BED"` (default).
#' @param genome Character(1), numeric(1), or `NULL`. MACS2 `-g` argument. If
#'   `NULL`, the value is inferred from `object@genomeAnnotation$genome` when
#'   possible.
#' @param shift Integer(1). MACS2 `--shift`.
#' @param extsize Integer(1). MACS2 `--extsize`.
#' @param qval Numeric(1). MACS2 `-q` in (0, 1]. Default is 0.05.
#' @param nThreads Integer(1) >= 1. Threads for `data.table::fwrite()` (default: 1).
#' @param extraArgs Character vector of additional MACS2 arguments (default: `NULL`).
#' @param keepStandard Logical(1). Whether to retain only standard chromosomes.
#' @param verbose Logical(1). Whether to print progress messages.
#' @return The input object with narrow peak information stored in `object@narrowPeak`.
#'
#' @importFrom methods is
#' @importFrom S4Vectors SimpleList
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom data.table data.table fwrite
#' @importFrom IRanges IRanges
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- callPeak(
#'     object = object,
#'     macs2Path = "/path/to/macs2",
#'     nThreads = 10
#'   )
#' }
#'
callPeak <- function(object = NULL,
                     macs2Path = NULL,
                     outputDir = NULL,
                     format = "BED",
                     genome = NULL,
                     shift = -100L,
                     extsize = 200L,
                     qval = 0.05,
                     nThreads = 1,
                     extraArgs = NULL,
                     keepStandard = TRUE,
                     verbose = TRUE){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  if(missing(object) || !methods::is(object, "whereAreECDNAs")){
    stop("'object' is required and it must be a 'whereAreECDNAs' object.", call. = FALSE)
  }

  sampleName = object@sampleMetadata$sample
  fragment <- object@fragment
  format <- toupper(format)

  if(!identical(format, "BED")){
    stop("callPeak() currently writes three-column BED input; please use format = 'BED'.", call. = FALSE)
  }

  if(!is.numeric(qval) || length(qval) != 1 || is.na(qval) || qval <= 0 || qval > 1){
    stop("'qval' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }

  ## Resolve the MACS2 executable path.
  if(is.null(macs2Path)){
    macs2Path <- Sys.which("macs2")
    if(macs2Path == ""){
      stop("Could not find 'macs2' on PATH. Please provide 'macs2Path'.")
    }
  }else{
    if(!file.exists(macs2Path)){
      stop("'macs2Path' does not exist: ", macs2Path)
    }
  }

  ## Resolve and create output directories.
  if(is.null(outputDir)){
    # Prefer the object-level output directory when it is available.
    od <- tryCatch(object@objectMetadata$outputDirectory, error = function(e) NULL)
    outputDir <- if(is.character(od) && length(od) == 1L && !is.na(od) && nzchar(od)) od else paste0(sampleName, ".outputs")
  }
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  outdir <- file.path(outputDir, "macs2")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  msg("[callPeak] Set MACS2 output directory: %s ", outputDir)

  ## Write pseudobulk fragments to BED format.
  msg("[callPeak] Writing fragments to a three-column BED file for MACS2 input...")
  ## BED3 columns: chromosome, start, and end.
  # Convert GRanges coordinates to a three-column table for MACS2.
  chrom <- as.character(GenomicRanges::seqnames(fragment))
  start0 <- as.integer(GenomicRanges::start(fragment))  # start coordinate
  end1 <- as.integer(GenomicRanges::end(fragment)) # end coordinate

  bedDt <- data.table::data.table(chrom, start0, end1)
  bedPath <- file.path(outputDir, paste0(sampleName, "_fragments.bed"))
  data.table::fwrite(bedDt, bedPath, sep = "\t", col.names = FALSE, nThread = nThreads)

  ## Build and execute the MACS2 command.
  msg("[callPeak] Running MACS2 to call peaks...")

  if(is.null(genome)){
    genome <- .inferMacs2Genome(object)
  }

  gArg <- if(is.numeric(genome)){
    format(genome, scientific = FALSE)
  }else{
    genome
  }

  args <- c(
    "callpeak",
    "-t", bedPath,
    "-n", sampleName,
    "--outdir", outdir,
    "-f", format,
    "-g", gArg,
    "--nomodel",
    "--shift", as.character(shift),
    "--extsize", as.character(extsize),
    "-q", as.character(qval)
  )

  if (!is.null(extraArgs) && length(extraArgs) > 0L) {
    args <- c(args, extraArgs)
  }

  res <- system2(macs2Path, args = args, wait = TRUE)
  exitCode <- attr(res, "status")
  if (is.null(exitCode)) exitCode <- 0L

  msg("[callPeak] MACS2 peak calling completed.")

  if (!identical(exitCode, 0L)) {
    warning("MACS2 returned a non-zero exit code (", exitCode, "). ",
            "Please inspect the MACS2 log in '", outdir, "'.")
  }

  ## Locate the narrowPeak output file.
  expected_narrowPeak <- file.path(outdir, paste0(sampleName, "_peaks.narrowPeak"))

  if(file.exists(expected_narrowPeak)){
    narrowPeakFile <- expected_narrowPeak
  } else {
    files <- list.files(path = outdir, pattern = "\\.narrowPeak$", full.names = TRUE)
    if(length(files) == 0L){
      stop("No '.narrowPeak' file was found in output directory: ", outdir)
    }
    if(length(files) > 1L){
      warning("Multiple '.narrowPeak' files found. Using the first one: ", basename(files[1L]))
    }
    narrowPeakFile <- files[1L]
  }

  ## Read MACS2 peaks and convert them to a GRanges object.
  peaksColnames <- c(
    "chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"
  )

  peakDf <- utils::read.table(
    file = narrowPeakFile,
    sep = "\t",
    col.names = peaksColnames,
    stringsAsFactors = FALSE
  )

  peakRanges <- IRanges::IRanges(
    start = peakDf$start,
    end = peakDf$end
  )

  peakGRanges <- GenomicRanges::GRanges(
    seqnames = peakDf$chrom,
    ranges = peakRanges,
    strand = rep("*", nrow(peakDf)),
    symbol = paste0(peakDf$chrom, ":", peakDf$start, "-", peakDf$end)
  )

  ## Load and validate genome annotation.
  genomeAnnotation <- setGenomeAnnotation()
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)

  if(isTRUE(keepStandard)){
    msg("[callPeak] Keeping standard chromosomes...")
    peakGRanges <- GenomeInfoDb::keepStandardChromosomes(peakGRanges, pruning.mode = "coarse")
  }

  attr(peakGRanges, "outdir") <- outdir
  object@narrowPeak <- peakGRanges

  msg("[callPeak] Peak calling completed. A total of %d narrow peaks have been added to 'object@narrowPeak'.", length(peakGRanges))
  return(object)
}
