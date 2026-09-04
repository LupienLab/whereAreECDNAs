"%ni%" <- Negate("%in%")
'%bcin%' <- function(x, table) S4Vectors::match(x, table, nomatch = 0) > 0
'%bcni%' <- function(x, table) !(S4Vectors::match(x, table, nomatch = 0) > 0)

## Apply a list of arguments through the package-safe lapply wrapper.
.batchlapply <- function(args = NULL, sequential = FALSE){

  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }
  args$subThreads <- 1
  args$threads <- 1

  args <- args[names(args) %ni% c("registryDir", "parallelParam", "subThreading")]
  outlist <- do.call(.safelapply, args)
  return(outlist)
}

.safelapply <- function(..., threads = 1, preschedule = FALSE){
  threads <- 1
  o <- lapply(...)
  return(o)
}

## Extract the file extension from a path-like string.
.fileExtension <- function (x = NULL){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

## Return the package-level verbose option, defaulting to TRUE when unset.
getWhereAreECDNAsVerbose <- function(){
  whereAreECDNAsVerbose <- options()[["whereAreECDNAs.verbose"]]
  if(!is.logical(whereAreECDNAsVerbose)){
    options(whereAreECDNAs.verbose = TRUE)
    return(TRUE)
  }
  whereAreECDNAsVerbose
}

## Suppress package startup messages, warnings, and regular messages.
.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

## Create a temporary file path inside a managed temporary directory.
.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){

  if(!dir.exists(tmpdir)){
    if(file.exists(tmpdir)){
      stop(paste0("Attempted to create temporary directory ", tmpdir," but a file already exists with this name. Please remove this file and try again!"))
    }
  }

  dir.create(tmpdir, showWarnings = FALSE)

  if(!dir.exists(tmpdir)){
    stop(paste0("Unable to create temporary directory ", tmpdir,". Check file permissions!"))
  }

  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }

  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))

}


## Detect density extrema for ecDNA content thresholding.
.findDensityExtrema <- function(object = NULL,
                                manualForce = NULL,
                                proportion = NULL,
                                verbose = TRUE){

  ## Define a formatted progress-message helper.
  msg <- function(...) if(isTRUE(verbose)) message(sprintf(...))

  if(missing(object) || is.null(object)){
    stop("'object' is required.", call. = FALSE)
  }

  if(!is.logical(manualForce) || length(manualForce) != 1){
    stop("'manualForce' must be a single logical value.", call. = FALSE)
  }

  if(!is.numeric(proportion) || length(proportion) != 2){
    stop("'proportion' must be a numeric vector of length 2.", call. = FALSE)
  }

  if(any(is.na(proportion)) || any(proportion <= 0) || any(proportion >= 1)){
    stop("'proportion' values must be between 0 and 1.", call. = FALSE)
  }

  if(sum(proportion) >= 1){
    stop("The sum of 'proportion' values must be less than 1.", call. = FALSE)
  }

  sampleName <- object@sampleMetadata$sample
  nth <- object@sampleMetadata$nth
  content <- as.numeric(object@cellColData$content)

  if(is.null(content) || all(is.na(content))){
    stop("No valid ecDNA content values were found. Please run quantifyECDNAContent() first.", call. = FALSE)
  }

  content <- content[is.finite(content)]

  if(length(content) < 10){
    stop("At least 10 cells are required to estimate the ecDNA content density.", call. = FALSE)
  }

  if(length(unique(content)) < 2){
    stop("At least two distinct ecDNA content values are required to estimate density-based thresholds.", call. = FALSE)
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

  ## Estimate a moderately smoothed density to stabilize summit detection.
  msg("[findDensityExtrema] Estimating ecDNA content density for %s-%s using %d cells.", sampleName, nth, length(content))

  densityBandwidthFactor <- 1.25
  densityPilot <- stats::density(content, n = 512)
  densityBandwidth <- densityPilot$bw * densityBandwidthFactor
  dens <- stats::density(content, bw = densityBandwidth, n = 512)
  dens_df <- data.frame(
    x = dens$x,
    y = dens$y
  )

  density_x <- dens_df$x
  density_y <- dens_df$y
  n_density <- length(density_y)
  maxValleyRatio <- 0.995
  minSummitSeparationSd <- 0.5
  minSummitHeight <- 0.03

  ## Identify local maxima and select the best separated pair of density modes.
  summit_idx <- which(
    density_y[2:(n_density - 1L)] >= density_y[1:(n_density - 2L)] &
      density_y[2:(n_density - 1L)] >= density_y[3:n_density] &
      (
        density_y[2:(n_density - 1L)] > density_y[1:(n_density - 2L)] |
          density_y[2:(n_density - 1L)] > density_y[3:n_density]
      )
  ) + 1L

  ## Remove negligible local maxima that arise from sparse distribution tails.
  summit_idx <- summit_idx[
    density_y[summit_idx] >= max(density_y) * minSummitHeight
  ]

  candidateModes <- NULL

  if(length(summit_idx) >= 2L && !isTRUE(manualForce)){

    summitPairs <- utils::combn(summit_idx, 2L)
    content_sd <- stats::sd(content)

    candidateModes <- lapply(seq_len(ncol(summitPairs)), function(i){

      left_idx <- summitPairs[1L, i]
      right_idx <- summitPairs[2L, i]
      between_idx <- seq.int(left_idx, right_idx)
      valley_idx <- between_idx[which.min(density_y[between_idx])]

      if(valley_idx %in% c(left_idx, right_idx)){
        return(NULL)
      }

      lower_summit_y <- min(density_y[left_idx], density_y[right_idx])
      valley_ratio <- density_y[valley_idx] / lower_summit_y
      separation_sd <- abs(density_x[right_idx] - density_x[left_idx]) / content_sd

      data.frame(
        left_idx = left_idx,
        right_idx = right_idx,
        valley_idx = valley_idx,
        left_x = density_x[left_idx],
        left_y = density_y[left_idx],
        right_x = density_x[right_idx],
        right_y = density_y[right_idx],
        valley_x = density_x[valley_idx],
        valley_y = density_y[valley_idx],
        valley_ratio = valley_ratio,
        separation_sd = separation_sd,
        score = lower_summit_y * (1 - valley_ratio) * log1p(separation_sd)
      )
    })

    candidateModes <- do.call(rbind, candidateModes[!vapply(candidateModes, is.null, logical(1))])

    if(!is.null(candidateModes) && nrow(candidateModes) > 0L){
      ## Retain mode pairs with a local minimum and adequate separation.
      candidateModes <- candidateModes[
        candidateModes$valley_ratio <= maxValleyRatio &
          candidateModes$separation_sd >= minSummitSeparationSd,
        ,
        drop = FALSE
      ]
    }
  }

  has_valley <- !is.null(candidateModes) && nrow(candidateModes) > 0L

  valley_x <- NA_real_
  valley_y <- NA_real_

  if(has_valley){
    selectedMode <- candidateModes[which.max(candidateModes$score), , drop = FALSE]
    valley_x <- selectedMode$valley_x
    valley_y <- selectedMode$valley_y
    msg("[findDensityExtrema] Selected a bimodal threshold model with valley at x = %.3f.", valley_x)
  }else{
    msg("[findDensityExtrema] No well-separated bimodal density structure was detected. Quantile-based thresholding will be used.")
  }

  ## Initialize threshold and plotting variables.
  left_cutoff <- NA_real_
  right_cutoff <- NA_real_
  method <- NA_character_

  left_summit <- NULL
  right_summit <- NULL
  left_coordinate <- NULL
  right_coordinate <- NULL

  if(has_valley && !isTRUE(manualForce)){

    left_summit <- data.frame(x = selectedMode$left_x, y = selectedMode$left_y)
    right_summit <- data.frame(x = selectedMode$right_x, y = selectedMode$right_y)
    left_cutoff <- left_summit$x
    right_cutoff <- right_summit$x
    method <- "bimodal"
    msg(
      "[findDensityExtrema] Selected density summits for putative training cells: negative <= %.3f; positive > %.3f.",
      left_cutoff,
      right_cutoff
    )
  }

  if(!has_valley || isTRUE(manualForce)){

    ## Use quantile-based thresholds when valley-based thresholds are unavailable.
    left_cutoff_x <- as.numeric(stats::quantile(content, probs = proportion[1], na.rm = TRUE))
    right_cutoff_x <- as.numeric(stats::quantile(content, probs = 1 - proportion[2], na.rm = TRUE))

    ## Interpolate density values at the selected cutoff coordinates.
    left_cutoff_y <- stats::approx(x = dens_df$x, y = dens_df$y, xout = left_cutoff_x, rule = 2)$y
    right_cutoff_y <- stats::approx(x = dens_df$x, y = dens_df$y, xout = right_cutoff_x, rule = 2)$y

    ## Store cutoff coordinates for plotting.
    left_coordinate <- data.frame(x = left_cutoff_x, y = left_cutoff_y)
    right_coordinate <- data.frame(x = right_cutoff_x, y = right_cutoff_y)

    left_cutoff <- left_cutoff_x
    right_cutoff <- right_cutoff_x
    method <- "quantile"

    msg("[findDensityExtrema] Quantile-based thresholds selected: left = %.3f, right = %.3f.", left_cutoff, right_cutoff)
  }

  ## Define density regions for visualization.
  left_region <- subset(dens_df, x <= left_cutoff)
  middle_region <- subset(dens_df, x > left_cutoff & x < right_cutoff)
  right_region <- subset(dens_df, x >= right_cutoff)

  ## Generate the density-threshold plot.
  p <- ggplot2::ggplot(dens_df, ggplot2::aes(x = x, y = y)) +

    ggplot2::geom_line(color = "#1f77b4", linewidth = 1) +
    ggplot2::labs(x = paste0(sampleName, "-", nth, " ecDNA content"), y = "Density") +

    ggplot2::geom_area(data = left_region, ggplot2::aes(x = x, y = y), fill = "#A6CEE3", alpha = 0.7) +
    ggplot2::geom_area(data = middle_region, ggplot2::aes(x = x, y = y), fill = "#D9D9D9", alpha = 0.7) +
    ggplot2::geom_area(data = right_region, ggplot2::aes(x = x, y = y), fill = "#CAB2D6", alpha = 0.7) +

    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  label_y <- max(dens_df$y) * 0.15

  if(method == "bimodal"){
    p <- p +

      ggplot2::annotate("segment",
                        x = left_summit$x, xend = left_summit$x,
                        y = 0, yend = left_summit$y,
                        color = "black", linetype = "dashed", linewidth = 0.25) +

      ggplot2::annotate("segment",
                        x = right_summit$x, xend = right_summit$x,
                        y = 0, yend = right_summit$y,
                        color = "black", linetype = "dashed", linewidth = 0.25) +

      ## Mark the valley-bottom point.
      ggplot2::geom_point(
        data = data.frame(x = valley_x, y = valley_y),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = "blue",
        size = 1.2) +

      ## Mark the density summit points.
      ggplot2::geom_point(
        data = data.frame(x = c(left_summit$x, right_summit$x),
                          y = c(left_summit$y, right_summit$y)),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = "blue",
        size = 1.2) +

      ## Label the valley coordinates.
      ggplot2::annotate(
        "text",
        x = valley_x,
        y = valley_y,
        label = paste0("(", round(valley_x, 3), ", ", round(valley_y, 3), ")"),
        size = 3,
        hjust = -0.15,
        vjust = 0.5
      ) +

      ## Label the left summit coordinates.
      ggplot2::annotate(
        "text",
        x = left_summit$x,
        y = left_summit$y,
        label = paste0("(", round(left_summit$x, 3), ", ", round(left_summit$y, 3), ")"),
        size = 3,
        hjust = -0.15,
        vjust = 0.5
      ) +

      ## Label the right summit coordinates.
      ggplot2::annotate(
        "text",
        x = right_summit$x,
        y = right_summit$y,
        label = paste0("(", round(right_summit$x, 3), ", ", round(right_summit$y, 3), ")"),
        size = 3,
        hjust = -0.15,
        vjust = 0.5
      ) +

      ggplot2::annotate("text",
                        x = left_summit$x,
                        y = label_y,
                        label = "ecDNA(-)",
                        size = 3,
                        hjust = 1.2) +

      ggplot2::annotate("text",
                        x = (left_summit$x + right_summit$x) / 2,
                        y = label_y,
                        label = "grey-zone",
                        size = 3,
                        hjust = 0.5) +

      ggplot2::annotate("text",
                        x = right_summit$x,
                        y = label_y,
                        label = "ecDNA(+)",
                        size = 3,
                        hjust = -0.2)
  }

  if(method == "quantile"){
    p <- p +

      ggplot2::annotate(
        "segment",
        x = left_coordinate$x, xend = left_coordinate$x,
        y = 0, yend = left_coordinate$y,
        color = "black", linetype = "dashed", linewidth = 0.25) +

      ggplot2::annotate(
        "segment",
        x = right_coordinate$x, xend = right_coordinate$x,
        y = 0, yend = right_coordinate$y,
        color = "black", linetype = "dashed", linewidth = 0.25) +

      ## Mark the threshold points.
      ggplot2::geom_point(
        data = data.frame(x = c(left_coordinate$x, right_coordinate$x),
                          y = c(left_coordinate$y, right_coordinate$y)),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = "blue",
        size = 1.2) +

      ## Label the left threshold coordinates.
      ggplot2::annotate(
        "text",
        x = left_coordinate$x,
        y = left_coordinate$y,
        label = paste0("(", round(left_coordinate$x, 3), ", ", round(left_coordinate$y, 3), ")"),
        size = 3,
        hjust = -0.15,
        vjust = 0.5
      ) +

      ## Label the right threshold coordinates.
      ggplot2::annotate(
        "text",
        x = right_coordinate$x,
        y = right_coordinate$y,
        label = paste0("(", round(right_coordinate$x, 3), ", ", round(right_coordinate$y, 3), ")"),
        size = 3,
        hjust = -0.15,
        vjust = 0.5
      ) +

      ggplot2::annotate("text",
                        x = left_coordinate$x,
                        y = label_y,
                        label = "ecDNA(-)",
                        size = 3,
                        hjust = 1.2) +

      ggplot2::annotate("text",
                        x = (left_coordinate$x + right_coordinate$x) / 2,
                        y = label_y,
                        label = "grey-zone",
                        size = 3,
                        hjust = -0.5) +

      ggplot2::annotate("text",
                        x = right_coordinate$x,
                        y = label_y,
                        label = "ecDNA(+)",
                        size = 3,
                        hjust = -0.2)
  }

  print(p)

  ggplot2::ggsave(paste0(outputDir, "/", sampleName, "-", nth, "_ecDNA_content_distribution.pdf"), plot = p, width = 10, height = 4)

  thresholdSummary <- data.frame(
    method = method,
    left = left_cutoff,
    right = right_cutoff,
    valley = valley_x,
    valley_density = valley_y
  )

  return(thresholdSummary)
}


## Infer a MACS2 effective genome-size shortcut from the object annotation.
.inferMacs2Genome <- function(object = NULL){

  genome <- tryCatch(
    object@genomeAnnotation$genome,
    error = function(e) NULL
  )

  if(is.null(genome) || length(genome) != 1L || is.na(genome) || !nzchar(genome)){
    warning("Could not infer the MACS2 genome size from object@genomeAnnotation; using 'hs'.")
    return("hs")
  }

  genomeLower <- tolower(as.character(genome))

  if(grepl("hg|hsapiens|homo|human", genomeLower)){
    return("hs")
  }

  if(grepl("mm|mmusculus|musculus|mouse", genomeLower)){
    return("mm")
  }

  warning("Could not map genome annotation '", genome, "' to a MACS2 shortcut; using 'hs'.")
  "hs"
}


## Detect compute resources and report available hardware.
.detectResource <- function(torch = NULL,
                            cpuForce = FALSE,
                            verbose = TRUE){

  ## Validate resource-detection parameters.
  if(!is.logical(cpuForce) || length(cpuForce) != 1){
    stop("'cpuForce' must be a single logical value.", call. = FALSE)
  }

  if(!is.logical(verbose) || length(verbose) != 1){
    stop("'verbose' must be a single logical value.", call. = FALSE)
  }

  osType <- unname(Sys.info()[["sysname"]])

  ## Limit selected thread variables on macOS to improve stability.
  if(osType == "Darwin"){
    threadVars <- c("OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS")

    changedVars <- character(0)

    for(v in threadVars){
      if(!nzchar(Sys.getenv(v))){
        do.call(Sys.setenv, stats::setNames(list("1"), v))
        changedVars <- c(changedVars, v)
      }
    }

    if(length(changedVars) > 0 && isTRUE(verbose)){
      message("[detectResource] macOS detected. The following thread variables were set to 1 for improved stability: ", paste(changedVars, collapse = ", "), ".")
    }
  }

  ## Detect CPU resources.
  nLogical <- NA_integer_
  nPhysical <- NA_integer_

  if(requireNamespace("parallel", quietly = TRUE)){
    nLogical <- parallel::detectCores(logical = TRUE)
    nPhysical <- parallel::detectCores(logical = FALSE)
  }

  ## Detect system RAM.
  ram <- .getSystemRam()

  ## Initialize device information.
  device <- "cpu"
  gpuInfo <- list()

  ## Detect GPU or accelerator backends unless CPU usage is forced.
  if(!isTRUE(cpuForce) && !is.null(torch)){

    ## Detect the Apple Silicon MPS backend.
    if(identical(osType, "Darwin")){
      mpsAvailable <- tryCatch(
        as.logical(torch$backends$mps$is_available()),
        error = function(e) FALSE
      )

      if(isTRUE(mpsAvailable)){
        device <- "mps"

        chipInfo <- tryCatch(
          system("sysctl -n machdep.cpu.brand_string", intern = TRUE),
          error = function(e) "Unknown"
        )

        gpuInfo <- list(
          type = "Apple Silicon",
          backend = "MPS",
          name = ifelse(length(chipInfo) > 0, chipInfo[1], "Unknown")
        )
      }
    }

    ## Detect the NVIDIA CUDA backend.
    if(identical(device, "cpu")){
      cudaAvailable <- tryCatch(
        as.logical(torch$cuda$is_available()),
        error = function(e) FALSE
      )

      if(isTRUE(cudaAvailable)){
        device <- "cuda"

        gpuCount <- tryCatch(
          as.integer(torch$cuda$device_count()),
          error = function(e) NA_integer_
        )

        gpuName <- tryCatch(
          torch$cuda$get_device_name(0L),
          error = function(e) "Unknown"
        )

        gpuMem <- tryCatch({
          props <- torch$cuda$get_device_properties(0L)
          round(as.numeric(props$total_memory) / (1024^3), 2)
        }, error = function(e) NA_real_)

        gpuInfo <- list(
          type = "NVIDIA",
          backend = "CUDA",
          name = gpuName,
          count = gpuCount,
          memory = gpuMem
        )
      }
    }
  }

  ## Print the resource-detection report.
  if(isTRUE(verbose)){
    cat("\n")
    cat("========================================\n")
    cat("System resource detection report\n")
    cat("========================================\n")
    cat(sprintf("Operating system: %s\n", osType))

    if(!is.na(nPhysical) && !is.na(nLogical)){
      cat(sprintf("CPU cores: %d physical / %d logical\n", nPhysical, nLogical))
    } else if(!is.na(nLogical)){
      cat(sprintf("CPU cores: %d logical\n", nLogical))
    } else {
      cat("CPU cores: detection failed\n")
    }

    if(!is.na(ram)){
      cat(sprintf("System RAM: %.2f GB\n", ram))
    } else {
      cat("System RAM: detection failed\n")
    }

    cat("----------------------------------------\n")

    if(isTRUE(cpuForce)){
      cat("Compute mode: CPU only, as specified by the user.\n")
    } else if(identical(device, "mps")){
      cat("Compute mode: GPU acceleration\n")
      cat("GPU backend: Apple Metal Performance Shaders (MPS)\n")
      cat(sprintf("GPU type: %s\n", gpuInfo$type))
      cat(sprintf("Chip model: %s\n", gpuInfo$name))
    } else if(identical(device, "cuda")){
      cat("Compute mode: GPU acceleration\n")
      cat("GPU backend: NVIDIA CUDA\n")
      cat(sprintf("GPU model: %s\n", gpuInfo$name))

      if(!is.na(gpuInfo$count)){
        cat(sprintf("GPU count: %d\n", gpuInfo$count))
      }

      if(!is.na(gpuInfo$memory)){
        cat(sprintf("GPU memory: %.2f GB\n", gpuInfo$memory))
      }
    } else {
      cat("Compute mode: CPU\n")
      cat("Note: No supported GPU backend was detected; processing may be slower.\n")
    }

    cat("========================================\n\n")
  }

  return(
    list(
      device = device,
      os = osType,
      cpuPhysical = nPhysical,
      cpuLogical = nLogical,
      ram = ram,
      gpuInfo = gpuInfo
    )
  )
}


## Get total system RAM.
.getSystemRam <- function(){

  osType <- Sys.info()["sysname"]

  ramBytes <- tryCatch({

    if(osType == "Windows"){

      ## Windows: query physical memory capacity using WMIC.
      ramStr <- system("wmic MemoryChip get Capacity", intern = TRUE)
      ramStr <- ramStr[nzchar(trimws(ramStr))]

      if(length(ramStr) > 1){
        ramValues <- suppressWarnings(as.numeric(ramStr[-1]))
        sum(ramValues, na.rm = TRUE)
      }else{
        NA_real_
      }

    } else if(osType == "Linux"){

      ## Linux: read total memory from /proc/meminfo.
      ramLine <- system("grep MemTotal /proc/meminfo", intern = TRUE)

      if(length(ramLine) > 0){
        ramKb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ramLine[1])))
        ramKb * 1024
      } else {
        NA_real_
      }

    } else if(osType == "Darwin"){

      ## macOS: query physical memory using sysctl.
      ramLine <- system("sysctl -n hw.memsize", intern = TRUE)

      if(length(ramLine) > 0){
        suppressWarnings(as.numeric(ramLine[1]))
      } else{
        NA_real_
      }

    } else{
      NA_real_
    }

  }, error = function(e) NA_real_)

  if(!is.na(ramBytes)){
    round(ramBytes / (1024^3), 2)
  } else {
    NA_real_
  }
}


## Run batched probability prediction for TabPFN.
.predict_proba_batched <- function(classifier, X, batchSize = 64, verbose = TRUE){

  X <- as.matrix(X)
  n <- nrow(X)
  if(is.null(n) || n == 0){
    stop("X must be a non-empty cell-by-peak matrix.")
  }

  idx_list <- split(seq_len(n), ceiling(seq_len(n) / batchSize))

  if(verbose){
    message(sprintf("[predictECDNACells] Running prediction in %d batches (batchSize = %d)",length(idx_list), batchSize))
  }

  proba_list <- vector("list", length(idx_list))

  for(i in seq_along(idx_list)){
    idx <- idx_list[[i]]

    if(verbose){
      message(sprintf("[predictECDNACells] Batch %d/%d: cells %d-%d", i, length(idx_list), min(idx), max(idx)))
    }

    X_batch <- X[idx, , drop = FALSE]

    ## Run prediction for the current batch.
    proba_batch <- classifier$predict_proba(X_batch)
    proba_batch <- as.matrix(proba_batch)

    rownames(proba_batch) <- rownames(X_batch)
    proba_list[[i]] <- proba_batch

    ## Release batch-specific objects before processing the next batch.
    rm(X_batch, proba_batch)
    gc()
  }

  proba_all <- do.call(rbind, proba_list)
  proba_all <- proba_all[rownames(X), , drop = FALSE]

  return(proba_all)
}


## Save an RDS file after verifying that the target path is writable.
.safeSaveRDS <- function(
    object = NULL,
    file = "",
    ascii = FALSE,
    version = NULL,
    compress = TRUE,
    refhook = NULL
){
  ## Attempt to save a small test object at the target path.
  testDF <- data.frame(a = 1, b = 2)
  canSave <- suppressWarnings(tryCatch({
    saveRDS(object = testDF, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
    TRUE
  }, error = function(x){
    FALSE
  }))
  if(!canSave){
    dirExists <- dir.exists(dirname(file))
    if(dirExists){
      stop("Cannot saveRDS. File Path : ", file)
    }else{
      stop("Cannot saveRDS because directory does not exist (",dirname(file),"). File Path : ", file)
    }
  }else{
    saveRDS(object = object, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
  }
}
