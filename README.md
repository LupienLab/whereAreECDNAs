# whereAreECDNAs

`whereAreECDNAs` is an R package designed to identify ecDNA-positive cells from single-cell ATAC-seq data.

The package uses chromatin accessibility signals within predefined ecDNA amplicon intervals to estimate the likelihood that individual cells harbor ecDNA.

## Installation

The package can be installed directly from GitHub using `devtools`:

```r
install.packages("devtools")
library(devtools)

install_github("LupienLab/whereAreECDNAs")
```

Alternatively, using `remotes`:

```r
install.packages("remotes")

remotes::install_github("LupienLab/whereAreECDNAs")
```

After installation, load the package with:

```r
library(whereAreECDNAs)
```

## Usage

A complete example workflow is provided in the `testDemo.R` script included in the repository.

The example demonstrates how to:

1. Prepare the single-cell ATAC-seq input data.
2. Define ecDNA amplicon regions.
3. Calculate ecDNA-associated chromatin accessibility signals for individual cells.
4. Estimate the probability that each cell is ecDNA-positive.
5. Assign ecDNA-positive and ecDNA-negative cell classifications.

To run the demonstration, download or clone the repository and execute:

```r
source("testDemo.R")
```

Before running the script, update the input file paths and analysis parameters according to your dataset.

## Input Data

The analysis generally requires:

* A single-cell ATAC-seq fragment file.
* A list of valid cell barcodes.
* Genomic coordinates defining the ecDNA amplicon intervals.
* Appropriate genome annotation and chromosome information.

All genomic coordinates should use the same reference genome assembly, such as `hg38`.

## Output

`whereAreECDNAs` generates cell-level results describing:

* The estimated ecDNA content of individual cells.
* The probability that each cell is ecDNA-positive.
* A binary ecDNA-positive or -negative classification for each cell.

These outputs can be integrated with cell-type annotations, malignant cell classifications, transcriptional states, copy-number profiles, or other single-cell metadata for downstream analyses.

## Citation

A formal citation for `whereAreECDNAs` will be provided upon publication.

Until then, please cite the GitHub repository:

> Lupien Lab. `whereAreECDNAs`: Identification of ecDNA-positive cells from single-cell ATAC-seq data. GitHub repository: https://github.com/LupienLab/whereAreECDNAs

## Issues and Support

Questions, bug reports, and feature requests can be submitted through the GitHub Issues page:

https://github.com/LupienLab/whereAreECDNAs/issues

When reporting an issue, please include:

* The R version.
* The package version or GitHub commit.
* A minimal reproducible example.
* The complete error message.
* Relevant session information generated using:

```r
sessionInfo()
```

## License

Please refer to the `LICENSE` file in the repository for licensing information.
