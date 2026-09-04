#' Generate pseudobulk fragments by filtering for valid cell barcodes
#'
#' @description
#' This function subsets the fragment `GRanges` object to retain only entries whose
#' `barcode` metadata matches barcodes present in `rownames(object@cellColData)`.
#' The resulting fragments represent a pseudobulk aggregation of high-quality cells
#' and are typically used for downstream analyses such as MACS2 peak calling.
#'
#' @param object A `whereAreECDNAs` object.
#' @param verbose Logical(1). Whether to print progress messages.
#'
#' @return The input `whereAreECDNAs` object with `@fragment` filtered to
#'   high-quality cell barcodes.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#'
#' @examples
#' \dontrun{
#'   object <- generatePseudobulk(object = object)
#' }
generatePseudobulk <- function(object, verbose = TRUE){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  ## Validate the input object and required fragment metadata.
  if(missing(object) || !methods::is(object, "whereAreECDNAs")){
    stop("'object' is required and must be a 'whereAreECDNAs' object.", call. = FALSE)
  }

  frag <- object@fragment
  if(is.null(frag)){
    stop("object@fragment is NULL. Cannot generate pseudobulk fragments without input fragments.", call. = FALSE)
  }
  if(!methods::is(frag, "GRanges")){
    stop("object@fragment must be a GRanges.", call. = FALSE)
  }

  if (!"barcode" %in% names(S4Vectors::mcols(frag))) {
    stop("object@fragment must contain a metadata column named 'barcode'.", call. = FALSE)
  }

  barcodes <- rownames(object@cellColData)
  if (is.null(barcodes) || length(barcodes) == 0L) {
    stop("rownames(object@cellColData) are empty - no barcodes to filter by.", call. = FALSE)
  }

  ## Filter fragments to high-quality cell barcodes.
  msg("[generatePseudobulk] Filtering fragments by %d high-quality barcodes...", length(barcodes))

  keep_idx <- frag$barcode %in% barcodes
  n_before <- length(frag)
  n_after <- sum(keep_idx, na.rm = TRUE)

  object@fragment <- frag[keep_idx]

  msg("[generatePseudobulk] Fragment filtering completed. Kept %d fragments in 'object@fragment'.", n_after)

  return(object)
}



