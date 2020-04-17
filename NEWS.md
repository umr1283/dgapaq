# dgapaq 0.1.8

## Minor improvements and fixes

* Fix issue with reference allele [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.1.7

## Minor improvements and fixes

* Add missing new line in summary [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.1.6

## Minor improvements and fixes

* Fix combined vcf files not created [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.1.5

## Minor improvements and fixes

* Fix merge of files after check against reference panel [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.1.4

## Minor improvements and fixes

* Fix condition when be files are not available [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).
* Add the no read key version of [HRC perl script](inst/perl/HRC-1000G-check-bim-NoReadKey.pl).

# dgapaq 0.1.3

## Minor improvements and fixes

* Reorder YAML header in `qc_vcf()` [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

# dgapaq 0.1.2

## Minor improvements and fixes

* Add `term::readkey` module in system requirement.
* `qc_plink()` now also work on non-binary PLINK files (*i.e.*, `.map` and `.ped`).
* Now imports `ggplot2 v3.3.0`.
* Update all packages listed in `Imports`.

# dgapaq 0.1.1

## Minor improvements and fixes

* Use `gt` for tables.

# dgapaq 0.1.0

## New features

* `qc_plink()` allows to compute quality-control of genotyping array (PLINK format) 
    using a [rmarkdown template](inst/rmarkdown/templates/qc_plink/skeleton/skeleton.Rmd).
* `qc_vcf()` allows to compute post-imputation quality-control report 
    using a default [rmarkdown template](inst/rmarkdown/templates/qc_vcf/skeleton/skeleton.Rmd).

