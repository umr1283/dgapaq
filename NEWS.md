# dgapaq 0.4.1

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Reorder code in setup chunk.

# dgapaq 0.4.0

## Major improvements

* In [rmarkdown templates](inst/rmarkdown/templates),
  + Now uses `data.table`.
  + Now uses `gt`.
  + Now uses `ggplot2`.
  + Improved text, figures and tables.
  
* In `R/qc_vcf.R` and `R/qc_plink.R`, 
  + Core update based on new [rmarkdown templates](inst/rmarkdown/templates).

# dgapaq 0.3.0

## Minor improvements and fixes

* In [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd),
  + Improve performance and memory consumption.
* In `R/qc_vcf.R`, 
  + Update imports based on [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.2.0

## Minor improvements and fixes

* In `DESCRIPTION`, 
  + Remove dependendy to libterm-readkey-perl.
  + Update all packages listed in `Imports`.
  + Now imports `ggplot2 v3.3.0`.
  + Now imports `gt`.
* Add the no read key version of [HRC perl script](inst/perl/HRC-1000G-check-bim-NoReadKey.pl).
* In [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd),
  + Reorder YAML header in `qc_vcf()`.
  + Fix condition when be files are not available.
  + Fix merge of files after check against reference panel.
  + Fix combined vcf files not created.
  + Fix issue with reference allele.
  + Fix a typos in the names of a parameters.
  + Fix gender check text printed when check was disabled.
  + Fix condition when fasta file is not provided.
  + Ensure missing values are moved before computing mean.
  + Ensure scale for shape has always the right number of values.
  + Add missing new line in summary.
  + Add condition when no fasta file is provided.
* In `R/qc_vcf.R`, 
  + Add filenames and URLs to download needed files.
* In `R/qc_plink.R`, 
  + Fix a typos in the names of a parameters.
  + Now also work on non-binary PLINK files (*i.e.*, `.map` and `.ped`).

# dgapaq 0.1.0

## New features

* `qc_plink()` allows to compute quality-control of genotyping array (PLINK format) 
    using a [rmarkdown template](inst/rmarkdown/templates/qc_plink/skeleton/skeleton.Rmd).
* `qc_vcf()` allows to compute post-imputation quality-control report 
    using a default [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

