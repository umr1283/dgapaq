#' tidy_vcf
#'
#' Correct missing genotype in VCF using corrected genotype matrix.
#'
#' @param raw_vcf A `character`. The path to the raw merged VCF (format .vcf or .vcf.gz).
#'     Default is `NULL`.
#' @param clean_geno A `character`. The path to the corrected genotype matrix (output of `check_genotype()`)
#'     with variants in rows and samples in column. The first column should be the unique ID of variants,
#'     which is composed of variant's chromosome, position, reference allele and alternative allele,
#'     separated by "_". Default is `NULL`.
#' @param col_variant A `character`. The column name of variant ID in `clean_geno`. Default is `var_id`.
#' @param check_dim A `logical`. Does the number of samples and variants should be checked in data from
#'     `raw_vcf` and `clean_geno`. Default is `TRUE`.
#' @param output_directory A `character`. The path to the output directory. Default is `NULL`.
#' @param output_name A `character`. The output VCF names. Default is `tidy.vcf`.
#' @param bin_path A `list(character)`. A list giving the binary path of `bcftools` and `tabix`.
#'     Final VCF will be ordered and compressed in gz format if `bcftools` is available, and will be
#'     indexed if `tabix` is available. Default is `NULL` for all binary paths.
#' @param nb_cores A `numeric`. The number of CPUs to use. Default is `1`.
#'
#' @importFrom parallel mclapply
#' @importFrom data.table `:=`
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
tidy_vcf <- function(
  raw_vcf = NULL,
  clean_geno = NULL,
  col_variant = "var_id",
  check_dim = TRUE,
  output_directory = NULL,
  output_name = "tidy.vcf",
  bin_path = list(
    bcftools = NULL,
    tabix = NULL
  ),
  nb_cores = 1
) {

  if (!file.exists(raw_vcf)) stop("Provided 'raw_vcf' file not exist!")
  if (!file.exists(clean_geno)) stop("Provided 'clean_geno' file not exist!")

  `#CHROM` = "ALT" = "CHROM" = "POS" = "REF" = "var_id" <- NULL

  ## get meta-info lines
  system(paste(
    ifelse(
      test = grepl("\\.gz", basename(raw_vcf)),
      yes = paste("head -c 10000", raw_vcf, "| zcat 2>/dev/null | grep '##'"),
      no = paste("grep '##'", raw_vcf)
    ),
    ">", file.path(output_directory, output_name)
  ))

  raw_vcf <- data.table::fread(raw_vcf, showProgress = FALSE, nThread = nb_cores)
  samples_vcf <- colnames(raw_vcf)[-c(1:9)]
  raw_vcf[, var_id := paste(`#CHROM`, POS, REF, ALT, sep = "_")]
  clean_geno <- data.table::fread(clean_geno, showProgress = FALSE, nThread = nb_cores)
  samples_geno <- setdiff(colnames(clean_geno), col_variant)

  if (check_dim) {
    ## check samples and variants in both data table
    if (!identical(samples_vcf, samples_geno)) {
      message(
        "The number of samples is different between 'raw_vcf' and 'clean_geno',",
        "only samples in common will be checked."
      )
    }
    if (nrow(raw_vcf) != nrow(clean_geno)) {
      variants <- intersect(raw_vcf[["var_id"]], clean_geno[[col_variant]])
      if (length(variants) == 0)
        stop("No variant in common! Does the clean_geno$col_variant column in format of 'CHR_POS_REF_ALT'?")
      message(
        "The number of variants is different between 'raw_vcf' and 'clean_geno',",
        "only variants in common will be checked."
      )
      raw_vcf <- raw_vcf[var_id %in% variants]
      variants <- raw_vcf[["var_id"]]
      clean_geno <- clean_geno[
        get(col_variant) %in% variants
      ][,
        (col_variant) := factor(get(col_variant), levels = variants)
      ][
        order(get(col_variant))
      ][,
        (col_variant) := as.character(get(col_variant))
      ]
    }
  }

  ## modify missing geno
  tmp <- do.call("cbind", parallel::mclapply(
    X = samples_vcf,
    mc.cores = min(nb_cores, length(samples_vcf)),
    function(iid) {
      if (iid %in% colnames(clean_geno)) {
        real_0 <- which(clean_geno[[iid]] == 0)
        iid_col <- raw_vcf[[iid]]
        iid_col[real_0] <- gsub(pattern = "(.{1}/.{1})", replacement = "0/0", 0, x = iid_col[real_0])
        return(iid_col)
      } else {
        return(raw_vcf[[iid]])
      }
    }))
  colnames(tmp) <- samples_vcf

  data.table::fwrite(
    cbind(raw_vcf[, 1:9], tmp), file = file.path(output_directory, output_name),
    sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  if (!is.null(bin_path[["bcftools"]])) {
    system(paste(
      bin_path[["bcftools"]], "sort",
      "--output-type z",
      "--output-file", paste0(file.path(output_directory, output_name), ".gz"),
      file.path(output_directory, output_name)
    ))
    if (!is.null(bin_path[["tabix"]])) {
      system(paste(
        bin_path[["tabix"]], "-f -p vcf", paste0(file.path(output_directory, output_name), ".gz")
      ))
    }
    unlink(file.path(output_directory, output_name))
  }
}
