#' Read cell barcodes
#'
#' @description
#' Reads a cell barcode file and returns the barcodes as a character vector.
#'
#' @param barcodePath Character(1). Path to a barcode file, such as a `.txt` file.
#' @param verbose Logical(1). Whether to print progress messages.
#' @return A character vector of cell barcodes.
#'
#' @importFrom  data.table  fread
#' @export
#'
#' @examples
#' \dontrun{
#'   barcode <- getBarcode(
#'     barcodePath = "/path/to/barcode"
#'   )
#'   class(barcode)  # "character"
#' }
#'
getBarcode <- function(barcodePath, verbose = TRUE){

  stopifnot(file.exists(barcodePath))
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  return(data.table::fread(barcodePath, header = FALSE)[[1]])
}

#' Read ecDNA amplicon intervals
#'
#' @description
#' Reads an ecDNA amplicon interval file and returns a `GRanges` object.
#'
#' @param ampliconPath Character(1). Path to a BED-like or TSV amplicon interval file.
#' @param verbose Logical(1). Whether to print progress messages.
#' @return A `GRanges` object with ranges defined by chromosome, start, and end columns.
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
#' @examples
#' \dontrun{
#'   amplicon <- getAmplicon(
#'     ampliconPath = "/path/to/amplicon"
#'   )
#'   class(amplicon)  # "GRanges"
#' }
#'
getAmplicon <- function(ampliconPath, verbose = TRUE){

  stopifnot(file.exists(ampliconPath))
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  return(fread(ampliconPath, header = FALSE) %>%
           as.data.frame() %>%
           {colnames(.) = c("chr", "start", "end");.} %>%
           GenomicRanges::makeGRangesFromDataFrame())
}

#' Read fragment intervals
#'
#' @description
#' Reads a single-cell ATAC-seq fragment file and returns a `GRanges` object.
#'
#' @param fragmentPath Character(1). Path to a fragment file (`.tsv` or `.tsv.gz`).
#'   The expected columns are chromosome, start, end, barcode, and optional count.
#' @param nThreads Integer(1). Number of threads for `data.table::fread()` (default: `1`).
#' @param keepStandard Logical(1). Whether to retain only standard chromosomes.
#' @param verbose Logical(1). Whether to print progress messages.
#' @return A `GRanges` object containing fragment intervals, with `barcode` and
#'   `count` metadata columns.
#'
#' @importFrom data.table fread
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges sort
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
#'
#' @examples
#' \dontrun{
#'   fragment <- getFragment(
#'     fragmentPath = "/path/to/fragment",
#'     nThreads = 1
#'   )
#'   class(fragment)  # "GRanges"
#' }
#'

getFragment <- function(fragmentPath,
                        nThreads = 1,
                        keepStandard = TRUE,
                        verbose = TRUE){

  stopifnot(is.character(fragmentPath), length(fragmentPath) == 1L)
  if (!file.exists(fragmentPath)) stop("File not found: ", fragmentPath)

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  ## Read fragment files with either four or five leading columns.
  ## Four columns: chromosome, start, end, barcode.
  ## Five columns: chromosome, start, end, barcode, count.
  ## Additional columns are ignored.
  colnames5 <- c("chr","start","end","barcode","count")
  colnames4 <- c("chr","start","end","barcode")

  ## Inspect the file header to determine the number of available columns.
  peek <- data.table::fread(fragmentPath, header = FALSE, nrows = 0, showProgress = isTRUE(verbose))
  ncols <- ncol(peek)
  if (ncols < 4L) stop("Fragment file must have at least 4 columns: chr, start, end, barcode")
  use_cols <- min(5L, ncols)
  use_names <- if (use_cols >= 5L) colnames5[seq_len(use_cols)] else colnames4[seq_len(use_cols)]
  colnames(peek) <- use_names

  msg("[getFragment] Loading the fragment file (columns detected: %d).", use_cols)
  fr_dt <- data.table::fread(
    file = fragmentPath,
    header = FALSE,
    select = seq_len(use_cols),
    col.names = use_names,
    nThread = nThreads,
    showProgress = isTRUE(verbose)
  )

  ## Ensure fragment coordinates are stored as integers.
  fr_dt[["start"]] <- as.integer(fr_dt[["start"]])
  fr_dt[["end"]] <- as.integer(fr_dt[["end"]])

  ## Remove rows with missing or invalid coordinates.
  fr_dt <- fr_dt[!is.na(chr) & !is.na(start) & !is.na(end) & end >= start]

  if(nrow(fr_dt) == 0L){
    msg("[getFragment] No valid fragments after coordinate QC; returning empty GRanges.")
    return(GenomicRanges::GRanges())
  }

  ## Convert the fragment table to a GRanges object.
  gr <- GenomicRanges::GRanges(
    seqnames = fr_dt[["chr"]],
    ranges   = IRanges::IRanges(start = fr_dt[["start"]], end = fr_dt[["end"]]),
    strand   = "*"
  )
  ## Attach barcode and count metadata.
  mcols(gr)$barcode <- fr_dt[["barcode"]]
  mcols(gr)$count <- as.integer(fr_dt[["count"]])

  ## Sort fragments to provide consistent downstream behavior.
  gr <- GenomicRanges::sort(gr, ignore.strand = TRUE)

  if(isTRUE(keepStandard)){
    msg("[getFragment] Keeping standard chromosomes...")
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  }

  msg("[getFragment] Fragment loading completed. A total of %d raw fragments were loaded. ", length(gr))
  gr
}
