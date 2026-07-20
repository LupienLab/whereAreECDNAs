rm(list = ls())
set.seed(seed = 17)

library(dplyr)
library(S4Vectors)
library(reticulate)
library(GenomicRanges)
library(whereAreECDNAs)

py_require(c("tabpfn","torch","numpy"))
tabpfn <- reticulate::import("tabpfn")
np <- reticulate::import("numpy")
torch <- reticulate::import("torch")

sample = "id" # Sample identifier
nth <- 1 # ecDNA species

# Step 0: Load input files
addwhereAreGenomes("hg38")

# Load ecDNA amplicon intervals, for example id_ecDNA_nth_intervals.bed
amplicon <- getAmplicon(ampliconPath = "/path/to/amplicon")

# Load high-quality cell barcode file, for example id.filtercell.txt
barcode <- getBarcode(barcodePath = "/path/to/barcode")

# Load fragment file, for example id.fragments.tsv.gz
fragment <- getFragment(fragmentPath = "/path/to/fragment", nThreads = 10)

# Step 1: Initialize a whereAreECDNAs object.
object <- initwhereAreObject(sample = sample,
                             nth = nth,
                             fragment = fragment,
                             barcode = barcode,
                             amplicon = amplicon)

# Step 2: Generate pseudobulk fragments filtered by high-quality cell barcodes.
object <- object %>% generatePseudobulk(.)

# Step 3: Call peaks from the pseudobulk BED file.
object <- object %>% callPeak(., macs2Path = macs2Path, nThreads = 10)

# Step 4: Build a weighted cell-by-peak accessibility matrix.
object <- object %>% buildCellPeakMatrix(., nThreads = 10)

# Step 5: Quantify ecDNA content for individual cells.
object <- object %>% quantifyECDNAContent(., tiesMethod = 'min')

# Step 6: Predict ecDNA-positive cells.
object <- object %>% predictECDNACells(.,
                                       tabpfn = tabpfn,
                                       numpy = numpy,
                                       torch = torch)

## Remove large intermediate slots before serialization.
object@fragment <- GRanges()
object@narrowPeak <- GRanges()
object@cellPeakMatrix <- NULL

# Step 7: Save the final object without fragment, narrowPeak, and cellPeakMatrix data.
saveRDS(object, file = paste0(sample, ".outputs/", sample, "-", nth, ".object.rds"))

