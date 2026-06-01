#' Build a weighted cell-by-peak accessibility matrix
#'
#' @description
#' Builds a weighted cell-by-peak accessibility matrix from single-cell
#' ATAC-seq fragments and peak regions. Fragment start/end positions are
#' discretized into genomic tiles, cell-by-tile coverage is computed, and
#' cell-by-peak accessibility is calculated using peak-specific weights.
#'
#' @details
#' The matrix is constructed in four steps. First, fragment start and end
#' positions are binned into fixed-width genomic tiles to represent local
#' insertion signal for each cell. Second, tiles overlapping each peak are
#' linked to that peak through a tile-by-peak mapping matrix and assigned a
#' peak-level weight. The default weighting scheme gives narrower peaks higher
#' weights, because focal accessibility events are more likely to represent
#' specific regulatory signal than broad regions. Weights are linearly scaled
#' from `1` for the widest peak to `peakScaleFactor` for the narrowest peak.
#' Third, blacklist regions are removed from the tile-by-peak mapping when
#' available. Finally, chromosome
#' matrices are merged and normalized within each cell by the total weighted
#' accessibility signal across candidate peaks, followed by scaling to
#' `scaleTo`. This produces a sparse cells-by-peaks accessibility matrix
#' suitable for downstream ecDNA content quantification.
#'
#' @param object A `whereAreECDNAs` object containing `fragment`, `narrowPeak`,
#'   and cell barcodes accessible as `rownames(object@cellColData)`.
#' @param tileSize Integer(1). Tile size, in base pairs, for discretizing fragment positions (default: `500`).
#' @param peakScaleFactor Numeric(1). Scaling factor for peak-width-based weighting (default: `5`).
#' @param scaleTo Numeric(1). Target scale factor applied after per-cell
#'   library-size normalization by total weighted peak accessibility
#'   (default: `10000`).
#' @param ceiling Numeric(1). Maximum value for count clipping (default: `8`).
#' @param excludeChr Character vector of chromosome names to exclude (default: `c("chrM")`).
#' @param blacklist Optional `GRanges` object containing blacklist regions to exclude.
#' @param tmpFile Character(1). Prefix/path for temporary `.rds` files (default: `NULL`, auto-generated).
#' @param nThreads Integer(1) >= 1. Number of threads (default: `1`).
#' @param verbose Logical(1). Whether to print progress messages.
#' @return The input `whereAreECDNAs` object with `object@cellPeakMatrix`
#'   populated as a sparse cells-by-peaks matrix.
#'
#' @importFrom methods is
#' @importFrom GenomicRanges GRanges GRangesList seqnames start end width resize sort ranges findOverlaps
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom S4Vectors mcols mcols<- match queryHits subjectHits
#' @importFrom BiocGenerics which
#' @importFrom GenomeInfoDb seqlevels seqlevels<- sortSeqlevels
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- buildCellPeakMatrix(
#'     object = object,
#'     nThreads = 10
#'   )
#' }

buildCellPeakMatrix <- function(
    object = NULL,
    tileSize = 500,
    peakScaleFactor = 5,
    scaleTo = 10000,
    ceiling = 8,
    excludeChr = c("chrM"),
    blacklist = NULL,
    tmpFile = NULL,
    nThreads = 1,
    verbose = TRUE
){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  if(missing(object) || !methods::is(object, "whereAreECDNAs")){
    stop("'object' is required and it must be a 'whereAreECDNAs' object.", call. = FALSE)
  }

  if(is.null(tmpFile)){
    tmpFile <- tempfile("cellPeakMat")
  }

  if(!is.numeric(tileSize) || length(tileSize) != 1 || tileSize <= 0){
    stop("'tileSize' must be a positive numeric scalar.", call. = FALSE)
  }

  if(!is.numeric(peakScaleFactor) || length(peakScaleFactor) != 1 || peakScaleFactor <= 0){
    stop("'peakScaleFactor' must be a positive numeric scalar.", call. = FALSE)
  }

  if(!is.numeric(scaleTo) || length(scaleTo) != 1 || scaleTo <= 0){
    stop("'scaleTo' must be a positive numeric scalar.", call. = FALSE)
  }

  if(!is.null(ceiling) && (!is.numeric(ceiling) || length(ceiling) != 1 || ceiling <= 0)){
    stop("'ceiling' must be NULL or a positive numeric scalar.", call. = FALSE)
  }

  if(!is.numeric(nThreads) || length(nThreads) != 1 || nThreads < 1){
    stop("'nThreads' must be a positive integer.", call. = FALSE)
  }

  ## Extract cell barcodes from the object metadata.
  if (!is.null(object@cellColData) && nrow(object@cellColData) > 0){
    cellBarcode <- rownames(object@cellColData)
  } else {
    stop("Cell barcodes could not be extracted from 'object@cellColData'.", call. = FALSE)
  }

  narrowPeak <- object@narrowPeak
  fragment <- object@fragment

  if(length(narrowPeak) == 0){
    stop("'object@narrowPeak' is empty.", call. = FALSE)
  }
  if(length(fragment) == 0){
    stop("'object@fragment' is empty.", call. = FALSE)
  }

  chrs <- unique(as.character(GenomicRanges::seqnames(fragment)))

  ## Filter narrow peaks to retained chromosomes with valid peak labels.
  narrowPeak <- narrowPeak[BiocGenerics::which(seqnames(narrowPeak) %bcni% excludeChr)]
  narrowPeak <- narrowPeak[BiocGenerics::which(seqnames(narrowPeak) %bcin% chrs)]
  seqlevels(narrowPeak) <- as.character(unique(seqnames(narrowPeak)))
  narrowPeak <- narrowPeak[!is.na(mcols(narrowPeak)$symbol) ]
  narrowPeak <- sort(sortSeqlevels(narrowPeak), ignore.strand = TRUE)

  if(length(narrowPeak) == 0){
    stop("No narrow peaks remain after chromosome and metadata filtering.", call. = FALSE)
  }

  ## Load and validate the genome annotation stored in the object.
  genomeAnnotation <- setGenomeAnnotation(object)
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)

  ## Store peak boundary coordinates.
  narrowPeak$peakStart <- start(GenomicRanges::resize(narrowPeak, 1, "start"))
  narrowPeak$peakEnd <- start(GenomicRanges::resize(narrowPeak, 1, "end"))

  ## Compute inverse-width peak weights scaled from 1 to peakScaleFactor.
  inverseWidth <- 1 / width(narrowPeak)
  inverseWidthRange <- range(inverseWidth)
  if(diff(inverseWidthRange) == 0 || peakScaleFactor == 1){
    narrowPeak$peakWeight <- rep(1, length(narrowPeak))
  }else{
    narrowPeak$peakWeight <- 1 + (inverseWidth - inverseWidthRange[1]) *
      (peakScaleFactor - 1) / diff(inverseWidthRange)
  }
  narrowPeak <- sort(sortSeqlevels(narrowPeak), ignore.strand = TRUE)

  ## Split peaks by chromosome for chromosome-level matrix construction.
  narrowPeak <- split(narrowPeak, seqnames(narrowPeak))
  narrowPeak <- lapply(narrowPeak, function(x){
    mcols(x)$idx <- seq_along(x)
    x
  })

  ## Load blacklist regions and split them by chromosome when available.
  if(is.null(blacklist) && !is.null(genomeAnnotation$blacklist)){
    blacklist <- genomeAnnotation$blacklist
  }

  if(!is.null(blacklist) && length(blacklist) > 0){
    blacklist <- split(blacklist, seqnames(blacklist))
  }else{
    blacklist <- NULL
  }

  cellLibrarySize <- .safelapply(seq_along(narrowPeak), function(z){

    chr <- names(narrowPeak)[z]
    msg("[buildCellPeakMatrix] Computing cell-by-peak accessibility submatrix for %s ...", chr)

    cellLibrarySizez <- tryCatch({

      narrowPeakz <- narrowPeak[[z]]
      narrowPeakz <- narrowPeakz[order(narrowPeakz$idx)]
      chrz <- paste0(unique(seqnames(narrowPeakz)))

      fragGRanges <- fragment[seqnames(fragment) == chrz]
      frag <- IRanges(start = start(fragGRanges),
                      end = end(fragGRanges))
      mcols(frag)$barcode <- mcols(fragGRanges)$barcode
      barcodeIdx <- S4Vectors::match(mcols(frag)$barcode, cellBarcode)
      frag <- frag[!is.na(barcodeIdx)]
      barcodeIdx <- barcodeIdx[!is.na(barcodeIdx)]
      fragmentStartTile <- trunc(start(frag)/tileSize) * tileSize
      fragmentEndTile <- trunc(end(frag)/tileSize) * tileSize
      fragmentCellIndex <- rep(barcodeIdx, 2)
      rm(frag); gc()

      uniqueTileStarts <- sort(unique(c(fragmentStartTile, fragmentEndTile)))
      cellbytileMat <- Matrix::sparseMatrix(
        i = fragmentCellIndex,
        j = match(c(fragmentStartTile, fragmentEndTile), uniqueTileStarts),
        x = rep(1, 2 * length(fragmentStartTile)),
        dims = c(length(cellBarcode), length(uniqueTileStarts))
      )
      if(!is.null(ceiling)){
        cellbytileMat@x[cellbytileMat@x > ceiling] <- ceiling
      }
      uniqueTiles <- IRanges(start = uniqueTileStarts, width = tileSize)
      rm(fragmentStartTile, fragmentEndTile, fragmentCellIndex, uniqueTileStarts); gc()

      peakRanges <- ranges(narrowPeakz)

      ## Identify overlaps between peaks and genomic tiles.
      tilePeakHits <- findOverlaps(uniqueTiles, peakRanges)
      tilePeakWeights <- mcols(narrowPeakz)$peakWeight[subjectHits(tilePeakHits)]

      ## Remove contributions from blacklisted tiles.
      blacklistChr <- if(!is.null(blacklist) && chrz %in% names(blacklist)){
        blacklist[[chrz]]
      }else{
        NULL
      }

      if(!is.null(blacklistChr) && length(blacklistChr) > 0){
        tileBlacklistIndicator <- 1 * (!overlapsAny(uniqueTiles, ranges(blacklistChr)))
        tilePeakWeights <- tilePeakWeights * tileBlacklistIndicator[queryHits(tilePeakHits)]
      }

      tilebypeakMat <- Matrix::sparseMatrix(
        i    = queryHits(tilePeakHits),
        j    = subjectHits(tilePeakHits),
        x    = tilePeakWeights,
        dims = c(ncol(cellbytileMat), length(narrowPeakz))
      )

      ## Build a cell-by-peak matrix.
      cellbypeakMat <- Matrix::drop0(cellbytileMat %*% tilebypeakMat)
      rownames(cellbypeakMat) <- cellBarcode

      ## Compute per-cell weighted peak signal for library-size normalization.
      cellLibrarySizez <- Matrix::rowSums(cellbypeakMat)

      ## Save the chromosome-specific matrix for subsequent normalization.
      .safeSaveRDS(cellbypeakMat, file = paste0(tmpFile, "-", chrz, ".rds"), compress = FALSE)
      rm(tilebypeakMat, cellbytileMat, cellbypeakMat, tilePeakWeights, tilePeakHits, uniqueTiles)
      gc()

      cellLibrarySizez

    }, error = function(e){
      ## Report chromosome-specific matrix construction errors.
      message("Error while building the cell-by-peak submatrix for chromosome: ", z)
      message(conditionMessage(e))
      NULL
    })
    cellLibrarySizez
  }, threads = nThreads) %>% Reduce("+", .)

  if(any(!is.finite(cellLibrarySize))){
    stop("Per-cell weighted peak library sizes contain non-finite values.", call. = FALSE)
  }

  zeroLibraryCells <- which(cellLibrarySize <= 0)
  if(length(zeroLibraryCells) > 0){
    msg("[buildCellPeakMatrix] %d cells have zero weighted peak signal and will remain zero after normalization.", length(zeroLibraryCells))
  }

  msg("[buildCellPeakMatrix] Merging chromosome-level submatrices and applying per-cell library-size normalization to scaleTo = %s.", scaleTo)

  cellbypeakMatList <- list()
  allPeakGR <- GRangesList()

  for(z in seq_along(narrowPeak)){
    chrz <- paste0(unique(seqnames(narrowPeak[[z]])))
    narrowPeakz <- narrowPeak[[z]]
    narrowPeakz <- narrowPeakz[order(narrowPeakz$idx)]
    allPeakGR[[chrz]] <- narrowPeakz

    cellbypeakMat <- readRDS(paste0(tmpFile, "-", chrz, ".rds"))
    file.remove(paste0(tmpFile, "-", chrz, ".rds"))

    ## Normalize each cell by total weighted peak signal, then scale to scaleTo.
    normalizationFactor <- rep(0, length(cellLibrarySize))
    positiveCells <- cellLibrarySize > 0
    normalizationFactor[positiveCells] <- scaleTo / cellLibrarySize[positiveCells]
    cellbypeakMat <- Matrix::Diagonal(x = normalizationFactor) %*% cellbypeakMat
    cellbypeakMat@x <- round(cellbypeakMat@x, 3)
    cellbypeakMat <- Matrix::drop0(cellbypeakMat)

    cellbypeakMatList[[chrz]] <- cellbypeakMat
    rm(cellbypeakMat); gc()
  }

  ## Merge chromosome-specific cell-by-peak matrices.
  cellPeakMat <- do.call(cbind, cellbypeakMatList)
  ## Combine chromosome-specific peaks into one GRanges object.
  allPeakGR <- unlist(allPeakGR)

  rownames(cellPeakMat) <- cellBarcode
  colnames(cellPeakMat) <- paste0(seqnames(allPeakGR), ":", start(allPeakGR), "-", end(allPeakGR))

  object@cellPeakMatrix <- cellPeakMat

  msg("[buildCellPeakMatrix] Matrix construction completed. A weighted cell-by-peak accessibility matrix has been added to 'object@cellPeakMatrix'.")

  return(object)
}
