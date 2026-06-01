#' Quantify ecDNA content using the Mann-Whitney U statistic
#'
#' @description
#' Quantifies cell-level ecDNA content using an ecDNA-restricted cell-by-peak
#' accessibility matrix. Peaks overlapping ecDNA amplicons are first selected,
#' and each cell is assigned an ecDNA content score based on a
#' Mann-Whitney U statistic.
#'
#' @param object A `whereAreECDNAs` object containing ecDNA amplicons,
#'   a cell-by-peak accessibility matrix, and cell metadata.
#' @param tiesMethod Character(1). Ties handling passed to `base::rank()` when
#' ranking accessibility values. Options include `"average"`, `"first"`, `"last"`,
#' `"random"`, `"max"`, and `"min"`. Default is `"min"`.
#' @param verbose Logical(1). Whether to print progress messages (default: `TRUE`).
#'
#' @return A `whereAreECDNAs` object with ecDNA content scores added to
#'   `object@cellColData$content`. The ecDNA-restricted cell-by-peak matrix is
#'   stored in `object@subcellPeakMatrix`.
#'
#' @importFrom methods is slotNames
#' @importFrom GenomicRanges GRanges findOverlaps seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom stats median
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- quantifyECDNAContent(
#'     object = object,
#'     tiesMethod = "min"
#'   )
#' }


quantifyECDNAContent <- function(
    object = NULL,
    tiesMethod = "min",
    verbose = TRUE){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  ## Validate the input object.
  if(missing(object) || !methods::is(object, "whereAreECDNAs")){
    stop("'object' is required and must be a 'whereAreECDNAs' object.", call. = FALSE)
  }

  validTieMethods <- c("average", "first", "last", "random", "max", "min")
  if(!tiesMethod %in% validTieMethods){
    stop("'tiesMethod' must be one of: ", paste(validTieMethods, collapse = ", "), call. = FALSE)
  }

  ## Extract sample metadata used for output files.
  sampleName <- object@sampleMetadata$sample
  nth <- object@sampleMetadata$nth

  if(is.null(sampleName) || !nzchar(sampleName)){
    stop("Sample name is missing from 'object@sampleMetadata$sample'.", call. = FALSE)
  }

  if(is.null(nth) || !nzchar(as.character(nth))){
    stop("ecDNA species index is missing from 'object@sampleMetadata$nth'.", call. = FALSE)
  }

  ## Validate ecDNA amplicon intervals.
  amplicon <- object@amplicon

  if(is.null(amplicon) || !methods::is(amplicon, "GRanges") || length(amplicon) == 0L){
    stop("object@amplicon must be a non-empty GRanges of ecDNA amplicons.", call. = FALSE)
  }

  ## Validate cell metadata and barcode names.
  cellBarcode <- rownames(object@cellColData)

  if(is.null(cellBarcode) || length(cellBarcode) == 0L){
    stop("`rownames(object@cellColData)` are missing/empty; barcodes are required.", call. = FALSE)
  }

  ## Validate the cell-by-peak matrix.
  if(!"cellPeakMatrix" %in% methods::slotNames(object)){
    stop("Slot 'cellPeakMatrix' not found in the object. Store your peak-by-cell matrix there.", call. = FALSE)
  }

  cellPeakMatrix <- object@cellPeakMatrix

  if(is.null(cellPeakMatrix)) {
    stop("object@cellPeakMatrix must be a non-null matrix (cells x peaks).", call. = FALSE)
  }

  if(is.null(colnames(cellPeakMatrix))){
    stop("object@cellPeakMatrix must have colnames like 'chr:start-end'.", call. = FALSE)
  }

  ## Resolve and create the output directory.
  outputDir <- tryCatch(
    object@objectMetadata$outputDirectory,
    error = function(e) NULL
  )

  if(is.null(outputDir) || !nzchar(outputDir)){
    outputDir <- paste0(sampleName, ".outputs")
  }

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  ## Convert peak names into a GRanges object.
  peakNames <- colnames(cellPeakMatrix)

  isValidPeakName <- grepl("^[^:]+:[0-9]+-[0-9]+$", peakNames)

  if(!all(isValidPeakName)){
    stop("All peak names in 'object@cellPeakMatrix' must follow the format 'chr:start-end'.", call. = FALSE)
  }

  peakSeqnames <- sub(":.*$", "", peakNames)
  peakStart <- as.integer(sub("^[^:]+:([0-9]+)-[0-9]+$", "\\1", peakNames))
  peakEnd <- as.integer(sub("^[^:]+:[0-9]+-([0-9]+)$", "\\1", peakNames))

  peak <- GenomicRanges::GRanges(
    seqnames = peakSeqnames,
    ranges = IRanges::IRanges(
      start = peakStart,
      end = peakEnd
    )
  )

  ## Identify peaks that overlap ecDNA amplicons.
  hits <- S4Vectors::queryHits(GenomicRanges::findOverlaps(query = peak, subject = amplicon))

  if(length(hits) == 0L){
    stop("[quantifyECDNAContent] No peaks overlap ecDNA amplicons. ecDNA content cannot be quantified.", call. = FALSE)
  }

  msg("[quantifyECDNAContent] Identified %d peaks overlapping ecDNA amplicons.",length(hits))

  subcellPeakMatrix <- cellPeakMatrix[, hits, drop = FALSE]

  msg("[quantifyECDNAContent] Constructed an ecDNA-restricted cell-by-peak accessibility matrix and added it to 'object@subcellPeakMatrix'.")
  object@subcellPeakMatrix <- subcellPeakMatrix

  ## Quantify ecDNA content with a rank-based statistic.
  msg("[quantifyECDNAContent] Quantifying the ecDNA content per cell using '%s' for tie handling.", tiesMethod)

  n_cells <- nrow(subcellPeakMatrix)
  n_peaks <- ncol(subcellPeakMatrix)

  ## Convert the sparse matrix to a numeric vector for global ranking.
  ## Higher accessibility values receive smaller ranks because of the negative sign.
  rankedValues <- rank(
    -as.vector(as.matrix(subcellPeakMatrix)),
    ties.method = tiesMethod
  )

  uMatrixRank <- matrix(
    rankedValues,
    nrow = n_cells,
    ncol = n_peaks
  )

  rownames(uMatrixRank) <- rownames(subcellPeakMatrix)
  colnames(uMatrixRank) <- colnames(subcellPeakMatrix)

  ## Compute the Mann-Whitney U statistic.
  uStatistic <- rowSums(uMatrixRank) - (n_peaks * (n_peaks + 1L)) / 2

  uScore <- 1 - uStatistic / (n_peaks * max(uMatrixRank))

  names(uScore) <- rownames(subcellPeakMatrix)

  ## Store ecDNA content scores in cell metadata.
  object@cellColData$content <- as.numeric(uScore[cellBarcode])

  msg("[quantifyECDNAContent] Quantification completed. Cell-level ecDNA content scores have been added to 'object@cellColData$content'.")

  return(object)
}
