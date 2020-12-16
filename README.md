
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DNA Genotyping Arrays Processing And Quality-Control <!--<img src="man/figures/dgapaq.png" align="right" width="120" />-->

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![GitHub
tag](https://img.shields.io/github/tag/mcanouil/dgapaq.svg?label=latest%20tag&include_prereleases)](https://github.com/mcanouil/dgapaq)
[![R build
status](https://github.com/umr1283/dgapaq/workflows/R-CMD-check/badge.svg)](https://github.com/umr1283/dgapaq/actions)
<!-- badges: end -->

## Installation

``` r
# Install dgapaq from CRAN:
install.packages("dgapaq")

# Or the the development version from GitHub:
# install.packages("remotes")
remotes::install_github("mcanouil/dgapaq")
```

## Overview

  - `qc_plink()` allows to compute quality-control of genotyping array
    (PLINK format) using a [rmarkdown
    template](inst/rmarkdown/templates/qc_plink/skeleton/skeleton.Rmd).
  - `qc_vcf()` allows to compute post-imputation quality-control report
    using a default [rmarkdown
    template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).
  - `convert_assembly()` allows to convert VCFs to target genome
    assembly using the software
    [CrossMap](https://crossmap.readthedocs.io/en/latest/).

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/omicsr/dgapaq/issues).  
For questions and other discussion, please contact the package
maintainer.

-----

Please note that this project is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md).  
By participating in this project you agree to abide by its terms.
