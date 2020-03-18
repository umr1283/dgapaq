# dgapaq (development version)

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

