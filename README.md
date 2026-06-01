This is a tremendous software tool, termed whereAreECDNAs, designed to identify ecDNA-positive cells from single-cell ATAC-seq data.

## How to install, two ways:

1. The easiest way to run is to download the singularity image from singularity hub:
```
<--link coming soon--> 
```
2. You could also install it is through the devtools package:
   
```
install.packages("devtools")
library(devtools)
install_github("LupienLab/whereAreECDNAs")
```

## How to run multiple ways:

1. With GPU based system using Singularity(recommended):
   
a) Command line
```
singularity exec --nv whereAreECDNAs.sif Rscript testDemo.R
```
b) Submit as a job on a slurm cluster:
```
module load singularity
sbatch -p gpu -J whereAreECDNAs_gpu --account myaccount_gpu --export=ALL --gres=gpu:3  -N 1 -c 1 --mem 20G -t 1-0 --wrap "singularity exec --nv whereAreECDNAs.sif Rscript testDemo.R”
```
2. With CPU only based system using Singularity:
   
a) Command line 
```
singularity exec whereAreECDNAs.sif Rscript testDemo.R
```
b) Submit as a job on a slurm cluster:
```
module load singularity
sbatch -p superhimem -J whereAreECDNAs --account myaccount --export=ALL -c 12 -N 1 --mem 90G -t 1-0 --wrap "singularity exec whereAreECDNAs.sif Rscript testDemo.R"
```
3. With GPU based system without Singularity 
```
Rscript testDemo.R
```
4. With CPU based system without Singularity
```
Rscript testDemo.R
```

## Steps required:
```
#1.Load genome
addwhereAreGenomes("hg38")

#2.Load ecDNA amplicon intervals 
amplicon <- getAmplicon(ampliconPath = "/path/to/amplicon")

#3.Load high-quality cell barcode file 
barcode <- getBarcode(barcodePath = "/path/to/barcode")

#4.Load fragments file
fragment <- getFragment(fragmentPath = "/path/to/fragment", nThreads = 10)

#5.Create a whereAreECDNAs object
object <- whereAreECDNAs(sample = sample,
                         nth = nth,
                         fragment = fragment,
                         cellBarcode = barcode,
                         amplicon = amplicon)

#6.Generate a pseudo-bulk BED file filtered by high-quality cell barcodes
object <- object %>% generatePseudoBulk(.)

#7.Call peaks on the pseudo-bulk BED file
object <- object %>% callPeak(., macs2Path = "/path/to/MACS2", nThreads = 10)

#8.Build a weighted cell-by-peak score matrix
object <- object %>% buildCellPeakMatrix(., nThreads = 10)

#9.Quantify the ecDNA content for individual cells
object <- object %>% quantifyECDNAContent(., tiesMethod = 'min')

#10.Predict ecDNA-positive cells
object <- object %>% predictECDNACells(., tabpfn = tabpfn, numpy = numpy, torch = torch)

#11.Save results
saveRDS(object, file = paste0(sample, ".outputs/", sample, "-", nth, ".object.rds"))
```

