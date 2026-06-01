.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Expected a valid GRanges object.")
  }
}

.validGenomeAnnotation <- function(genomeAnnotation = NULL){

  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list or SimpleList containing blacklist GRanges, chromSizes GRanges, and a genome BSgenome package string, such as hg38 or BSgenome.Hsapiens.UCSC.hg38.")
    }
  }

  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){

    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    gA$genome <- .validGenomeString(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])

  }else{
    stop("genomeAnnotation must be a list or SimpleList containing blacklist GRanges, chromSizes GRanges, and a genome BSgenome package string, such as hg38 or BSgenome.Hsapiens.UCSC.hg38.")
  }
  return(gA)
}

.validGenomeString <- function(genome = NULL){

  if(inherits(genome, "BSgenome")){
    return(genome@pkgname)
  }

  if(is.character(genome) && length(genome) == 1L && !is.na(genome) && nzchar(genome)){
    return(genome)
  }

  stop("'genome' must be a non-empty character string or a BSgenome object.", call. = FALSE)
}

validBSgenome <- function(genome = NULL, masked = FALSE){

  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    if(!requireNamespace("BSgenome", quietly = TRUE)){
      stop("The BSgenome package is required to load genome objects.", call. = FALSE)
    }
    return(BSgenome::getBSgenome(genome, masked = masked))
  }else{
    stop("'genome' must be a BSgenome object or a character value accepted by BSgenome::getBSgenome().")
  }
}
