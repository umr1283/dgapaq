
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DNA Genotyping Arrays Processing And Quality-Control

<!-- badges: start -->

[![Lifecycle:
questioning](https://img.shields.io/badge/lifecycle-questioning-orange.svg)]()
[![GitHub
tag](https://img.shields.io/github/tag/umr1283/dgapaq.svg?label=latest%20tag&include_prereleases)](https://github.com/umr1283/dgapaq)
[![R build
status](https://github.com/umr1283/dgapaq/workflows/R-CMD-check/badge.svg)](https://github.com/umr1283/dgapaq/actions)
<!-- badges: end -->

## Installation

``` r
# Install dgapaq from CRAN:
install.packages("dgapaq")

# Or the the development version from GitHub:
# install.packages("remotes")
remotes::install_github("umr1283/dgapaq")
```

## Overview

-   `qc_plink()` allows to compute quality-control of genotyping array
    (PLINK format) using a [rmarkdown
    template](inst/rmarkdown/templates/qc_plink/skeleton/skeleton.Rmd).
-   `qc_vcf()` allows to compute post-imputation quality-control report
    using a default [rmarkdown
    template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).
-   `convert_assembly()` allows to convert VCFs to target genome
    assembly using the software
    [CrossMap](https://crossmap.readthedocs.io/en/latest/).
-   `compress_coverage()` allows to compress coverage file (output of
    `samtools depth`) into contiguous segments based on position.
-   `create_genotype_matrix()` allows to create a genotype matrix based
    on VCFs.
-   `check_genotype()` allows to check missing data in genotype matrix
    against coverage information.
-   `tidy_vcf()` allows to correct missing genotypes in VCF using
    corrected genotype matrix.

------------------------------------------------------------------------

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/umr1283/dgapaq/issues).  
For questions and other discussion, please contact the package
maintainer.

## Code of Conduct

Please note that the dgapaq project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
