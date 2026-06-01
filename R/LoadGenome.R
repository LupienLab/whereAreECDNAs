#' Set or query the default genome build
#'
#' This helper sets a global option used by `whereAreECDNAs` to determine
#' the default genome build (e.g., `"hg19"`, `"hg38"`, `"mm9"`, `"mm10"`).
#'
#' If `genome` is `NULL`, the function leaves the current option unchanged and
#' reports the active genome, when one is set, along with the supported genomes.
#'
#' @param genome Character(1) specifying the genome build. Must be one of
#'   `"hg19"`, `"hg38"`, `"mm9"`, or `"mm10"` (case-insensitive).
#'   If `NULL`, no change is made and the current setting is reported.
#'
#' @return Invisibly returns the canonical genome string (e.g., `"hg38"`)
#'   if a new value is set; otherwise invisibly returns the currently set
#'   genome (or `NULL` if unset).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   ## Query the current genome and supported builds
#'   addwhereAreGenomes()
#'
#'   ## Set genome to hg38
#'   addwhereAreGenomes("hg38")
#' }
#'
addwhereAreGenomes <- function(genome = NULL){

  supportedGenomes <- c("hg19", "hg38", "mm9", "mm10")
  optName <- "whereAreECDNAs.genome"

  ## Report the current and supported genomes when no new genome is supplied.
  if(is.null(genome)){
    current <- getOption(optName, default = NULL)
    if(is.null(current)){
      message("[addwhereAreGenomes] No default genome is currently set.")
    }else{
      message("[addwhereAreGenomes] Current default genome: ", current)
    }
    message("[addwhereAreGenomes] Supported genomes: ", paste(supportedGenomes, collapse = ", "))
    return(invisible(current))
    }

  ## Validate the requested genome build.
  if(!is.character(genome) || length(genome) != 1L || is.na(genome) || !nzchar(genome)){
    stop("'genome' must be a non-empty Character(1).", call. = FALSE)
    }

  genomeLower <- tolower(genome)
  if(!genomeLower %in% supportedGenomes){
    stop("Genome '", genome, "' is not currently supported by whereAreECDNAs.\n",
         "Supported genomes: ", paste(supportedGenomes, collapse = ", "), call. = FALSE)
    }
  options(structure(list(genomeLower), .Names = optName))

  message("[addwhereAreGenomes] Setting default genome to: ", genome, ".")
  invisible(genome)
}


#' Retrieve genome annotation
#'
#' If an object is provided, returns its `@genomeAnnotation` slot. If `object`
#' is `NULL`, retrieves the default annotation configured by
#' [addwhereAreGenomes()]. An error is raised when no default genome is set.
#'
#' @param object A `whereAreECDNAs` object, or `NULL`.
#'
#' @return The genome annotation object.
#'
#' @importFrom methods is
#' @export
#'
#' @examples
#' \dontrun{
#'   ## From an object
#'   annotation <- setGenomeAnnotation(object)
#'
#'   ## From the globally configured default genome
#'   annotation <- setGenomeAnnotation()
#' }
#'
setGenomeAnnotation <- function(object = NULL){
  ## Retrieve the globally configured annotation when no object is supplied.
  if(is.null(object)){
    genomeAnnotation <- getwhereAreGenomes(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation)
    }
    stop(
      "setGenomeAnnotation: no object provided and no default genome set.\n",
      "Please either:\n",
      " * provide a 'whereAreECDNAs' object with a valid '@genomeAnnotation', or\n",
      " * call addwhereAreGenomes(\"hg38\") (or another supported genome) first.",
      call. = FALSE
    )
  }

  ## Validate the object before returning its genome annotation.
  if(!methods::is(object, "whereAreECDNAs")){
    stop("'object' must be a 'whereAreECDNAs' object or NULL.", call. = FALSE)
  }
  object@genomeAnnotation
}

#' Get the default genome or its annotation
#'
#' This helper reads the global option set by [addwhereAreGenomes()] and returns
#' either the genome label (e.g., `"Hg38"`) or the corresponding genome
#' annotation object stored in package data (e.g., `genomeAnnoHg38`).
#'
#' @param genomeAnnotation Logical(1). If `FALSE` (default), return the genome
#'   label as a string (e.g., `"Hg38"`). If `TRUE`, loads and returns the
#'   corresponding genome annotation object.
#'
#' @return
#'   - If `genomeAnnotation = FALSE`: a character scalar with the genome label
#'     (e.g. `"Hg38"`), or `NULL` if no default genome is set.
#'   - If `genomeAnnotation = TRUE`: the genome annotation object stored in
#'     `genomeAnno*`, or `NULL` if no default genome is set.
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#'   # Return the genome label
#'   getwhereAreGenomes()
#'
#'   # Return the annotation object
#'   anno <- getwhereAreGenomes(genomeAnnotation = TRUE)
#' }
#'
getwhereAreGenomes <- function(genomeAnnotation = FALSE){

  supportedGenomes <- c("hg19", "hg38", "mm9", "mm10")
  optName <- "whereAreECDNAs.genome"

  ## Validate the return-type flag.
  if(!is.logical(genomeAnnotation) || length(genomeAnnotation) != 1L || is.na(genomeAnnotation)){
    stop("'genomeAnnotation' must be TRUE/FALSE.", call. = FALSE)
  }

  ## Read the configured genome option.
  ag <- getOption(optName, default = NULL)
  if(is.null(ag)){
    ## No default genome has been configured.
    return(NULL)
  }

  if(!is.character(ag) || length(ag) != 1L || is.na(ag) || !nzchar(ag)){
    stop(
      "Option '", optName, "' is not a valid non-empty string (found: ",
      paste0(deparse(ag), collapse = " "), ").",
      call. = FALSE
    )
  }

  genomeLower <- tolower(ag)
  if(!genomeLower %in% supportedGenomes){
    stop(
      "Option '", opt_name, "' is set to '", ag, "', which is not supported.\n",
      "Supported genomes: ", paste(supportedGenomes, collapse = ", "),
      call. = FALSE
    )
  }

  ## Format the genome label, for example "Hg38".
  genomeLabel <- paste0(
    toupper(substr(genomeLower, 1L, 1L)),
    substr(genomeLower, 2L, nchar(genomeLower))
  )

  if(!isTRUE(genomeAnnotation)){
    ## Return only the genome label when annotation data are not requested.
    return(genomeLabel)
  }

  ## Load the annotation dataset corresponding to the configured genome.
  datasetName <- paste0("genomeAnno", genomeLabel)

  # Load the dataset into this function's environment when needed.
  if(!exists(datasetName, envir = environment(), inherits = FALSE)){
    utils::data(list = datasetName, envir = environment())
  }

  if(!exists(datasetName, envir = environment(), inherits = FALSE)){
    stop(
      "Could not find genome annotation dataset '", dataset_name, "'. ",
      "Ensure it is installed and exported as package data.",
      call. = FALSE
    )
  }

  get(datasetName, envir = environment(), inherits = FALSE)
}

