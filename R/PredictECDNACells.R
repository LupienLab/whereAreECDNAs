#' Predict ecDNA-positive cells using TabPFN
#'
#' @description
#' Predicts ecDNA-positive and ecDNA-negative cells from an ecDNA-restricted
#' cell-by-peak accessibility matrix using a TabPFN classifier. Putative
#' ecDNA-positive and ecDNA-negative training cells are first defined based on
#' density-derived ecDNA content thresholds. The trained classifier is then used
#' to estimate ecDNA-positivity probabilities for all cells.
#'
#' @details
#' Putative training labels are derived from the distribution of cell-level
#' ecDNA content. When the density shows two well-separated modes, the function
#' identifies the dominant low-content and high-content summits and the lowest
#' intervening valley. Cells at or below the low-content summit are treated as
#' putative ecDNA-negative training cells, and cells above the high-content
#' summit are treated as putative ecDNA-positive training cells. This provides
#' conservative training labels while leaving the intermediate region
#' unlabelled. If a clear bimodal structure is not detected, or if
#' `manualForce = TRUE`, the function uses the lower and upper quantiles
#' specified by `proportion`.
#'
#' @param object A `whereAreECDNAs` object.
#' @param tabpfn Python module for `tabpfn`, imported through `reticulate`.
#' @param numpy Python module for `numpy`, imported through `reticulate`.
#' @param torch Python module for `torch`, imported through `reticulate`.
#' @param manualForce Logical(1). Whether to use manually defined training thresholds
#'   (default: `FALSE`).
#' @param proportion Numeric vector of length 2. The lower and upper tail
#'   proportions used to define putative negative and positive training cells
#'   when quantile-based thresholds are used. Default is `c(0.1, 0.1)`.
#' @param batchSize Integer(1). Number of cells processed per inference batch (default: `500`).
#' @param cpuForce Logical(1). Whether to force TabPFN to run on CPU (default: `FALSE`).
#' @param proThreshold Numeric(1). Probability threshold used to classify cells as ecDNA-positive (default: `0.5`).
#' @param seed Integer(1). Random seed passed to TabPFN (default: `1`).
#' @param verbose Logical(1). Whether to print progress messages (default: `TRUE`).
#'
#' @return An updated `whereAreECDNAs` object with two additional columns in
#'   `object@cellColData`: `probability` and `status`.
#'
#' @importFrom methods is
#' @importFrom utils write.table
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter select pull
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- predictECDNACells(
#'     object = object,
#'     tabpfn = tabpfn,
#'     numpy = numpy,
#'     torch = torch
#'   )
#' }

predictECDNACells <- function(object = NULL,
                              tabpfn = NULL,
                              numpy = NULL,
                              torch = NULL,
                              manualForce = FALSE,
                              proportion = c(0.1, 0.1),
                              batchSize = 500,
                              cpuForce = FALSE,
                              proThreshold = 0.5,
                              seed = 2026,
                              verbose = TRUE
                            ){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  ## Validate the input object.
  if(missing(object) || !methods::is(object, "whereAreECDNAs")){
    stop("'object' is required and must be a 'whereAreECDNAs' object.", call. = FALSE)
  }

  ## Validate required Python modules.
  if(missing(tabpfn) || is.null(tabpfn)){
    stop("tabpfn is required. Please import it first:\n",
         "library(reticulate)\n",
         "tabpfn <- reticulate::import('tabpfn')")
  }

  if(missing(numpy) || is.null(numpy)){
    stop("numpy is required. Please import it first:\n",
         "numpy <- reticulate::import('numpy')")
  }

  if(missing(torch) || is.null(torch)){
    stop("torch is required. Please import it first:\n",
         "torch <- reticulate::import('torch')")
  }

  ## Validate model and threshold parameters.
  if(!is.numeric(proportion) || length(proportion) != 2 || any(proportion <= 0) || any(proportion >= 1)){
    stop("'proportion' must be a numeric vector of length 2 with values between 0 and 1.", call. = FALSE)
  }

  if(!is.numeric(batchSize) || length(batchSize) != 1 || batchSize < 1){
    stop("'batchSize' must be a positive integer.", call. = FALSE)
  }

  if(!is.numeric(proThreshold) || length(proThreshold) != 1 || proThreshold < 0 || proThreshold > 1){
    stop("'proThreshold' must be a numeric value between 0 and 1.", call. = FALSE)
  }

  batchSize <- as.integer(batchSize)
  seed <- as.integer(seed)

  ## Extract sample metadata used for output files.
  sampleName <- object@sampleMetadata$sample
  nth <- object@sampleMetadata$nth

  if(is.null(sampleName) || !nzchar(sampleName)){
    stop("Sample name is missing from 'object@sampleMetadata$sample'.", call. = FALSE)
  }

  if(is.null(nth) || !nzchar(as.character(nth))){
    stop("ecDNA species index is missing from 'object@sampleMetadata$nth'.", call. = FALSE)
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

  ## Define training thresholds from the ecDNA content distribution.
  extrema <- .findDensityExtrema(object = object, manualForce = manualForce, proportion = proportion, verbose = verbose)

  ## Prepare the input matrix and cell metadata.
  msg("[predictECDNACells] Preparing an ecDNA-restricted cell-by-peak matrix: %d cells x %d peaks.", dim(object@subcellPeakMatrix)[1], dim(object@subcellPeakMatrix)[2])

  X_raw <- as.matrix(object@subcellPeakMatrix)
  meta <- as.data.frame(object@cellColData) %>% tibble::rownames_to_column(var = "barcode")

  if(!all(rownames(X_raw) %in% meta$barcode)){
    stop("[predictECDNACells] Some matrix cell barcodes are missing from object@cellColData.", call. = FALSE)
  }

  if(!"content" %in% colnames(meta)){
    stop("[predictECDNACells] The column 'content' is missing from object@cellColData.", call. = FALSE)
  }

  ## Define putative training cells from content thresholds.
  negBarcodes <- meta %>% dplyr::filter(content <= extrema$left) %>% dplyr::select(barcode) %>% pull(barcode)
  posBarcodes <- meta %>% dplyr::filter(content > extrema$right) %>% dplyr::select(barcode) %>% pull(barcode)

  if(length(negBarcodes) == 0){
    stop("[predictECDNACells] No putative ecDNA-negative training cells were identified")
  }
  if(length(posBarcodes) == 0){
    stop("[predictECDNACells] No putative ecDNA-positive training cells were identified")
  }

  ## Build the supervised training dataset.
  X_train <- rbind(X_raw[negBarcodes, , drop = FALSE], X_raw[posBarcodes, , drop = FALSE])
  y_train <- c(rep(0L, length(negBarcodes)), rep(1L, length(posBarcodes)))

  ## Convert training data to Python-compatible arrays.
  X_train <- numpy$array(X_train, dtype = "float32")
  y_train <- numpy$array(y_train, dtype = "int64")

  ## Allow TabPFN to process large datasets on CPU when needed.
  Sys.setenv(TABPFN_ALLOW_CPU_LARGE_DATASET = "1")

  ## Initialize and train the TabPFN classifier.
  resource <- .detectResource(torch, cpuForce = cpuForce)

  msg("[predictECDNACells] Running TabPFN using %s ...", toupper(as.character(resource$device)))

  classifier <- tabpfn$TabPFNClassifier(
    device = as.character(resource$device),
    memory_saving_mode = TRUE,
    ignore_pretraining_limits = TRUE,
    random_state = as.integer(seed)
  )

  msg("[predictECDNACells] Training TabPFN with putative %d ecDNA-positive and %d ecDNA-negative cells ...", length(posBarcodes), length(negBarcodes))

  classifier$fit(X_train, y_train)

  ## Predict ecDNA-positivity probabilities for all cells.
  msg("[predictECDNACells] Performing probability inference for all %d cells ...", nrow(X_raw))

  proba_all <- .predict_proba_batched(
    classifier = classifier,
    X = X_raw,
    batchSize = batchSize,
    verbose = verbose
  )

  probability <- as.numeric(proba_all[, 2])
  names(probability) <- rownames(X_raw)

  msg("[predictECDNACells] Classifying cells as ecDNA-positive or ecDNA-negative using probability threshold: %.2f.", proThreshold)
  status <- ifelse(probability > proThreshold, "ecDNA+", "ecDNA-")
  names(status) <- rownames(X_raw)

  ## Store probabilities and classifications in cell metadata.
  object@cellColData$probability <- probability[rownames(object@cellColData)]
  object@cellColData$status <- status[rownames(object@cellColData)]

  ## Write cell classification outputs to disk.
  positive <- sort(probability[probability > proThreshold], decreasing = TRUE)
  negative <- sort(probability[probability <= proThreshold], decreasing = TRUE)

  posfile <- file.path(outputDir, sprintf("%s-%s_ecDNA_positive_cell.txt", sampleName, nth))
  negfile <- file.path(outputDir,sprintf("%s-%s_ecDNA_negative_cell.txt", sampleName, nth))

  utils::write.table(
    names(positive),
    file = posfile,
    row.names = FALSE, col.names = FALSE,
    sep = "\t", quote = FALSE
  )

  utils::write.table(
    names(negative),
    file = negfile,
    row.names = FALSE, col.names = FALSE,
    sep = "\t", quote = FALSE
  )

  cellColData_df <- object@cellColData %>% as.data.frame() %>% tibble::rownames_to_column(var = "barcode")
  cellColDatafile <- file.path(outputDir, sprintf("%s-%s_cellColData.txt", sampleName, nth))

  utils::write.table(
    cellColData_df,
    file = cellColDatafile,
    row.names = FALSE, col.names = TRUE,
    sep = "\t", quote = FALSE
  )

  msg("[predictECDNACells] Probability estimation completed. Cell-level ecDNA probabilities and status have been added to 'object@cellColData'.")

  return(object)
}
