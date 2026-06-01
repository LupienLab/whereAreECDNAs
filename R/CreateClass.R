setClassUnion("characterOrNull", c("character", "NULL"))
setClassUnion("GRangesOrNull", c("GRanges", "NULL"))
setClassUnion("dgCMatrixOrNull", c("dgCMatrix", "NULL"))
setClassUnion("dgCMatrixOrMatrixOrNull", c("dgCMatrix", "matrix", "NULL"))

setClass("whereAreECDNAs",
         representation(
           objectMetadata = "SimpleList",
           sampleMetadata = "SimpleList",
           fragment = "GRangesOrNull",
           cellColData = "DataFrame",
           amplicon = "GRangesOrNull",
           genomeAnnotation = "SimpleList",
           narrowPeak = "GRangesOrNull",
           cellPeakMatrix = "dgCMatrixOrNull",
           subcellPeakMatrix = "dgCMatrixOrNull",
           anyVector = "characterOrNull",
           anyMatrix = "dgCMatrixOrMatrixOrNull"
         )
)

setMethod(
  "show", "whereAreECDNAs",
  function(object){

    ## Define a helper for concise object summaries.
    scat <- function(fmt, vals = character(), exdent = 2, n = 5, ...){

      vals <- as.character(vals)
      if(length(vals) == 0L){lbls <- "none"}
      else{
        lbls <- paste(Biobase::selectSome(vals, maxToShow = n), collapse = " ")
      }
      txt <- sprintf(fmt, length(vals), lbls)
      cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
    }

    ## Confirm that sample metadata is available before printing the summary.
    tryCatch({
      object@cellColData$sample
    }, error = function(e){
              stop(paste0("\nError accessing sample info from a 'whereAreECDNAs' object.",
                          "\nPlease report this issue on GitHub: https://github.com/LupienLab/whereAreECDNAs/issues"
              ),
            call. = FALSE
          )
      })

    cat("class:", class(object), "\n")
    cat("sample name:", object@sampleMetadata$sample, "\n")
    cat("ecDNA species:", paste0(object@sampleMetadata$sample, "-", object@sampleMetadata$nth), "\n")
    cat("number of cells:", nrow(object@cellColData), "\n")
    scat("column of cellColData (%d): %s\n", names(object@cellColData))
    cat("output directory:", object@objectMetadata$outputDirectory, "\n")
  }
)

#' Create a `whereAreECDNAs` object
#'
#' @description
#' Initializes a `whereAreECDNAs` S4 object containing sample metadata,
#' high-quality cell barcodes, fragments, ecDNA amplicons, genome annotation,
#' and output-directory information. The object is the central container for
#' downstream analysis steps, including pseudobulk fragment generation, peak
#' calling, weighted cell-by-peak matrix construction, ecDNA content
#' quantification, and ecDNA-positive cell prediction.
#'
#' @param sample Character(1). Sample name.
#' @param nth Character(1). Identifier for a specific ecDNA amplicon or species.
#' @param fragment A `GRanges` storing fragments with at least a `barcode` metadata column (and optionally `count`).
#' @param barcode Character vector of high-quality cell barcodes for this sample.
#' @param amplicon A `GRanges` representing ecDNA amplicon intervals.
#' @param nThreads Integer(1) >= 1. Number of threads reserved for downstream steps.
#' @param outputDirectory Character(1). Output directory. If `NULL`, defaults to `paste0(sample, ".outputs")`.
#' @param genomeAnnotation Genome annotation object returned by `setGenomeAnnotation()`.
#'   The value is validated with `.validGenomeAnnotation()`.
#' @param verbose Logical(1). Whether to print progress messages.
#' @return A `whereAreECDNAs` S4 object.
#'
#' @importFrom methods new is
#' @importFrom S4Vectors SimpleList DataFrame Rle
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix dgCMatrix
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- initwhereAreObject(
#'     sample = "sample",
#'     nth = "nth",
#'     fragment = fragment,
#'     barcode = barcode,
#'     amplicon = amplicon
#'   )
#' }
initwhereAreObject <- function(
    sample = NULL,
    nth = NULL,
    fragment = NULL,
    barcode = NULL,
    amplicon = NULL,
    nThreads = 1L,
    outputDirectory = NULL,
    genomeAnnotation = setGenomeAnnotation(),
    verbose = TRUE){

  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))
  stopifnot(is.character(sample), is.character(barcode))

  if(is.null(outputDirectory)){
    outputDirectory <- paste0(sample, ".outputs")
  }

  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)

  if(grepl(" ", outputDirectory)){
    stop("outputDirectory cannot have a space in the path! Path : ", outputDirectory)
  }
  dir.create(outputDirectory, showWarnings = FALSE)
  if (grepl(" ", normalizePath(outputDirectory))) {
    stop("outputDirectory cannot have a space in the full path! Full path : ", normalizePath(outputDirectory))
  }

  ## Assemble sample-level and cell-level metadata.
  sampleMetadata <- S4Vectors::SimpleList()
  sampleMetadata[["sample"]] <- sample
  sampleMetadata[["nth"]] <- nth

  metaData <- DataFrame(barcode)
  colnames(metaData) <- "barcode"
  metaData$sample <- Rle(sample, nrow(metaData))
  metaData$index <- 1:nrow(metaData)
  rownames(metaData) <- metaData$barcode
  metaData <- metaData[, -which(colnames(metaData) == "barcode")]
  metadataList <- metaData

  allCols <- unique(c("sample", rev(sort(colnames(metadataList)))))
  cellColData <- metadataList

  msg("[initwhereAreObject] Initializing a whereAreECDNAs object.")
  object <- methods::new(
    "whereAreECDNAs",
    objectMetadata = S4Vectors::SimpleList(outputDirectory = normalizePath(outputDirectory)),
    sampleMetadata = sampleMetadata,
    fragment = fragment,
    cellColData = cellColData,
    amplicon = amplicon,
    genomeAnnotation = genomeAnnotation,
    anyVector = NULL,
    anyMatrix = NULL
  )

  msg("[initwhereAreObject] A new whereAreECDNAs object has been initialized.")

  return(object)
}

#' @importFrom utils .DollarNames
#' @export
".DollarNames.whereAreECDNAs" <- function(x, pattern = ''){
  cpan <- as.list(c("barcode",colnames(x@cellColData)))
  names(cpan) <- c("barcode",colnames(x@cellColData))
  return(.DollarNames(x = cpan, pattern = pattern))
}

#' @export
"$.whereAreECDNAs" <- function(x, i){
  if(i=="barcode"){
    return(rownames(x@cellColData))
  } else {
    val <- x@cellColData[[i, drop = TRUE]]
    return(as.vector(val))
  }
}

#' @export
"$<-.whereAreECDNAs" <- function(x, i, value) {
  if(i == "sample"){
    stop("sample is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(i == "barcode"){
    stop("barcode is the column name in cellColData. Please do not try to overwrite this!")
  }
  if(object.size(Rle(value)) < 2 * object.size(value)){
    value <- Rle(value)
  }
  if(!is.null(value)){
    if(length(value) == 1) {
      value <- Rle(value, lengths = nrow(x@cellColData))
    }
  }
  x@cellColData[[i]] <- value
  return(x)
}
