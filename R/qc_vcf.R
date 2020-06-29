#' Compute post-imputation quality-control report using a default rmarkdown template.
#'
#' @param input_directory A `character`. The path to the VCF files.
#' @param output_directory A `character`. The path to the output directory.
#' @param cohort_name A `character`. The name of the studied cohort / population.
#' @param output_file A `character`. The name of the html file produced.
#' @param title A `character`. The report's title. Default is `"Post-Imputation Quality-Control"`.
#' @param author_name A `character`. The author's name to be printed in the report.
#'     Default is `Unknown`.
#' @param author_affiliation A `character`. The affiliation to be printed in the report.
#'     Default is `NULL`.
#' @param author_email A `character`. The email to be printed in the report.
#'     Default is `NULL`.
#' @param encoding A `character`. The encoding to be used for the html report.
#'     Default is `"UTF-8"`.
#' @param vcftools_path A `character`. Path to vcftools software binary.
#' @param ... Parameters to pass to `rmarkdown::render()`.
#'
#' @return NULL
#'
#' @importFrom bookdown html_document2
#'
#' @import data.table
#' @import ggplot2
#' @import gt
#' @import patchwork
#'
#' @importFrom knitr opts_chunk
#' @importFrom parallel mclapply detectCores
#' @importFrom scales comma percent viridis_pal comma_format percent_format
#' @importFrom sessioninfo session_info
#'
#' @export
qc_vcf <- function(
  input_directory = NULL,
  output_directory = NULL,
  cohort_name = "COHORT",
  output_file = paste0(cohort_name, "_imputation_QC.html"),
  title = "Post-Imputation Quality-Control",
  author_name = "Unknown",
  author_affiliation = NULL,
  author_email = NULL,
  encoding = "UTF-8",
  vcftools_path = "/usr/bin/vcftools",
  ...
) {
  message_prefix <- "[dgapaq] "

  message(message_prefix, "Quality-Control started ...")

  file.copy(
    from = system.file("rmarkdown", "templates", "qc_impute", "skeleton.Rmd", package = "umr1283"),
    to = file.path(tempdir(), "qc_impute.Rmd"),
    overwrite = TRUE
  )

  rmarkdown::render(
    input = file.path(tempdir(), "qc_vcf.Rmd"),
    output_file = output_file,
    output_dir = output_directory,
    encoding = encoding,
    params = list(
      input_directory = input_directory,
      output_directory = output_directory,
      title = title,
      author_name = author_name,
      author_affiliation = author_affiliation,
      author_email = author_email,
      vcftools_path = vcftools_path
    ),
    ...
  )

  message(message_prefix, "Quality-Control ended!")

  invisible()
}
