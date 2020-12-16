#' create_genotype_matrix
#'
#' Creates a raw genotype matrix with provided VCF files.
#'
#' @param sample_sheet A `data.frame`. A data frame with two mandatory columns `id` and `vcf_file`.
#' @param output_name A `character`. The prefix for output genotype data. Default is "raw_merged"
#' @param output_directory A `character`. The path to the output directory. Default is `NULL`.
#' @param chromosomes A vector of `character`. Chromosomes to keep
#'  (same nomenclature as chromosome name in VCF).
#' If argument not specified (as default `NULL`), this filter will be turned off.
#' @param remove_multiallelic A `logical`. Does the multi-allelic variants (i.e.: genotype composed
#' of different alternative alleles) should be removed. Default is `TRUE`.
#' @param split_multiallelic A `logical`. Does the multi-allelic variants should be splitted into
#' bi-allelic variants. Default is `FALSE`.
#' @param split_type A `character`. The type of variants should be splitted by `bcftools norm`. Default is `both`.
#' @param ref_fasta A `character`. A path to the reference fasta file, it should be provided
#' if multiallelic INDELs need left-alignment after splitting. Default is `NULL` which lead to no left-alignment.
#' @param bin_path A `list(character)`. A list giving the binary path of `vcftools`, `bcftools`, `tabix` and `bgzip`.
#' @param nb_cores An `integer`. The number of CPUs to use. Default is `1`.
#'
#' @return
#' A genotype matrix coded NA, 1, or 2 (additive model), with variants in rows and individuals in
#' columns. The first column `var_id` is composed of variant's chromosome, position, reference allele
#' and alternative allele, separated by "_".
#'
#' @importFrom parallel mclapply
#' @importFrom data.table `:=`
#' @importFrom data.table `.SD`
#' @importFrom data.table setDTthreads
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
create_genotype_matrix <- function(
  sample_sheet,
  output_name = "raw_merged",
  output_directory = NULL,
  chromosomes = NULL,
  remove_multiallelic = TRUE,
  split_multiallelic = FALSE,
  split_type = "both",
  ref_fasta = NULL,
  bin_path = list(
    vcftools = "/usr/bin/vcftools",
    bcftools = "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix"
  ),
  nb_cores = 1L
) {

  if (! all(c("id", "vcf_file") %in% colnames(sample_sheet))) {
    stop("Either 'id' or 'vcf_file' is missing in the sample_sheet.")
  }

  if (is.null(output_directory)) {
    stop("The output_directory must be specified.")
  }

  if (all(remove_multiallelic, split_multiallelic)) {
    stop("'remove_multiallelic' and 'split_multiallelic' can't be both TRUE!")
  }

  if (split_multiallelic & ! split_type %in% c("both", "snps", "indels", "any")) {
    stop("'split_type' should be either 'both', 'snps' or 'indels'.")
  }

  `#CHROM` = "ALT" = "CHROM" = "FORMAT" = "POS" = "REF" = "chr" = "idx" = "var_id" <- NULL

  output_tmp_dir <- file.path(output_directory, "tmp")
  dir.create(path = output_tmp_dir, showWarnings = FALSE, recursive = TRUE)
  check_exist <- sapply(sample_sheet$vcf_file, file.exists)
  if (!all(check_exist)) {
    message("Following sample's vcfs do not exist: ", paste(sample_sheet$id[!check_exist], collapse = ", "))
    sample_sheet <- sample_sheet[check_exist, ]
  }

  invisible(
    parallel::mclapply(X = sample_sheet$id, mc.cores = nb_cores, function(iid){
      ivcf <- sample_sheet$vcf_file[sample_sheet$id %in% iid]
      iout <- file.path(output_tmp_dir, paste0(iid, "_tmp.vcf.gz"))
      is_gzvcf <- grepl("vcf\\.gz$", ivcf)

      if (remove_multiallelic) {
        ## remove variants with double mutated alleles, i.e.: geno 1/2
        system(paste(
          bin_path[["vcftools"]],
          paste0("--", if (is_gzvcf) "gvcf" else "vcf"), ivcf,
          "--min-alleles 2 --max-alleles 2 --recode --stdout |",
          bin_path[["bcftools"]], "sort --output-type z --output-file", iout
        ))
      } else if (split_multiallelic) {
        ## split multiallelic variants, e.g.: ref:A, alt:T,G, geno:1/2 => A/T (0/1) and A/G (0/1)
        system(paste(
          bin_path[["bcftools"]],
          "norm",
          if (is.null(ref_fasta)) "" else paste("--fasta-ref", ref_fasta),
          "-m", paste0("-", split_type), ivcf, "|",
          bin_path[["bcftools"]], "sort --output-type z --output-file", iout
        ))
      } else {
        if (!is_gzvcf) {
          system(paste(
            bin_path[["bcftools"]], "sort --output-type z --output-file", iout, ivcf
          ))
        } else {
          system(paste("cp", ivcf, iout))
        }
      }

      ## indexing vcf
      system(paste(bin_path[["tabix"]], "-f -p vcf", iout))
    })
  )

  ## merge vcf
  vcfs_to_merge <- list.files(output_tmp_dir, pattern = "_tmp.vcf.gz$", full.names = TRUE)
  merged_put <- file.path(output_tmp_dir, paste0(output_name, ".vcf.gz"))
  if (length(vcfs_to_merge) > 1) {
    system(paste(
      bin_path[["bcftools"]],
      "merge -m none",
      paste0(vcfs_to_merge, collapse = " "),
      "-O z >", merged_put
    ))
  } else {
    file.copy(from = vcfs_to_merge, to = merged_put, overwrite = TRUE)
  }

  ## prepare raw matrix
  data.table::setDTthreads(threads = nb_cores)
  geno_mat <- data.table::fread(merged_put, header = TRUE)
  if (!is.null(chromosomes)) {
    geno_mat <- geno_mat[`#CHROM` %in% chromosomes]
  }
  geno_mat <- geno_mat[, var_id := paste(`#CHROM`, POS, REF, ALT, sep = "_")]
  geno_mat <- geno_mat[, .SD, .SDcols = !`#CHROM`:FORMAT]
  patern_geno_sep <- ifelse(grepl(".{1}/.{1}:", geno_mat[[1]][1]), "/", "|")
  geno_mat <- geno_mat[
    ,
    (sample_sheet$id) := lapply(
      .SD,
      FUN = gsub,
      pattern = paste0("(.{1}\\", patern_geno_sep, ".{1}):.*"),
      replacement = "\\1"),
    .SDcols = sample_sheet$id
  ]

  ## possible genotypes
  possible_genotypes_chr <- apply(
    X = expand.grid(c(".", "0", "1"), c(".", "0", "1")),
    MARGIN = 1,
    FUN = paste,
    collapse = patern_geno_sep
  )
  possible_genotypes_num <- rep(1, length(possible_genotypes_chr))
  names(possible_genotypes_num) <- possible_genotypes_chr
  possible_genotypes_num[grepl(paste0("0\\", patern_geno_sep, "0"), names(possible_genotypes_num))] <- 0
  possible_genotypes_num[grepl(paste0("1\\", patern_geno_sep, "1"), names(possible_genotypes_num))] <- 2
  possible_genotypes_num[grepl("\\.", names(possible_genotypes_num))] <- NA
  possible_genotypes_num <- c(possible_genotypes_num, "0" = 0, "1" = 1, "2" = 2)

  ## recode genotype
  geno_mat <- geno_mat[
    ,
    (sample_sheet$id) := lapply(.SD, FUN = function(x) possible_genotypes_num[x]),
    .SDcols = sample_sheet$id
  ][, .SD, .SDcols = c("var_id", sample_sheet$id)]

  data.table::fwrite(
    geno_mat,
    file = file.path(output_directory, paste0(output_name, ".tsv.gz")),
    col.names = TRUE
  )

  unlink(output_tmp_dir, recursive = TRUE)
  invisible(geno_mat)
}


#' check_genotype
#'
#' Checks the missing genotype against coverage information and impute missing genotype to 0 if
#' coverage is sufficient, otherwise missing genotype is recoded to -1.
#'
#'
#' @param sample_sheet A `data.frame`. A data frame with one mandatory column `id`.
#' The column `vcf_file` (path to vcf files) should be provided if you need to create the
#' raw genotype matrix. Please see details in `create_genotype_matrix()`.
#' The column `cov_file` (path to coverage files) should be provided if you need to compress the
#' coverage file (output of `samtools depth`) into contiguous segments. Please see details in `compress_coverage()`.
#' @param genotype_matrix A `character` or `NULL`. The path to the genotype matrix to be imported,
#' with samples in columns and variants in rows, the first columns should be `var_id` which is
#' unique for each variant. `NULL` leading the creation of raw genotype matrix with `create_genotype_matrix()`.
#' @param min_depth An `integer`. The minimum depth to filter when compressing coverage file. Default is `8`.
#' @param ... Optional arguments to `create_genotype_matrix()`
#' @param chromosomes A vector of `character`. Chromosomes to check
#'  (same nomenclature as chromosome name in VCF).
#' If argument not specified (default is `NULL`), this filter will be turned off.
#' @param output_directory A `character`. The path to the output directory. Default is `NULL`.
#' @param compress_cov_directory A `character`. The path the to directory of compressed coverage files. Default is `NULL`.
#' @param output_name A `character`. The prefix for output genotype data. Default is "clean_merged".
#' @param nb_cores An `integer`. The number of CPUs to use. Default is `1`.
#'
#' @return
#' A genotype matrix coded 0, 1, or 2 (additive model), with variants in rows and individuals in
#' columns. Missing genotype is coded to -1. The first column `var_id` is composed of variant's
#' chromosome, position, reference allele and alternative allele, separated by "_".
#'
#' @importFrom parallel mclapply
#' @importFrom data.table setDTthreads
#' @importFrom data.table `:=`
#' @importFrom data.table `.SD`
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table rbindlist
#' @importFrom data.table tstrsplit
#' @importFrom data.table inrange
#'
#' @export
check_genotype <- function(
  sample_sheet,
  genotype_matrix = NULL,
  min_depth = 8,
  ...,
  chromosomes = NULL,
  output_directory = NULL,
  compress_cov_directory = NULL,
  output_name = "clean_merged",
  nb_cores = 1L
) {

  if (! "id" %in% colnames(sample_sheet)) {
    stop("'id' is missing in the sample_sheet.")
  }

  if (is.null(output_directory)) {
    stop("The output_directory must be specified.")
  }

  if (is.null(compress_cov_directory)) {
    stop("The 'compress_cov_directory' should be specified.")
  }

  "CHROM" = "POS" = "chr" = "idx" = "var_id" <- NULL

  dir.create(path = output_directory, showWarnings = FALSE, recursive = TRUE)
  output_name <- file.path(output_directory, paste0(output_name, ".tsv.gz"))
  check_log <- file.path(output_directory, "check_geno.log")
  cat(paste("Start at", Sys.time()), file = check_log, sep = "\n")

  samples_ids <- sample_sheet[["id"]]
  data.table::setDTthreads(threads = nb_cores)

  analysis_time <- system.time({
    ## check .cov.compress file
    cat("Checking samples' compressed coverage files...", file = check_log, append = TRUE, sep = "\n")

    not_compressed <- !file.exists(file.path(compress_cov_directory, paste0(samples_ids, ".cov.compress")))
    if (any(not_compressed)) {
      cat(paste(
        "Start compressing the coverage file for sample:",
        paste0(sample_sheet$id[not_compressed], collapse = ", ")
      ), file = check_log, append = TRUE, sep = "\n")

      to_compress <- sample_sheet[not_compressed, ]
      compress_coverage(
        sample_sheet = to_compress,
        output_directory = compress_cov_directory,
        min_depth = min_depth,
        chromosomes = chromosomes,
        nb_cores = nb_cores
      )
    }

    ## raw genotype matrix
    cat("Getting samples' genotype matrix...", file = check_log, append = TRUE, sep = "\n")
    if (is.null(genotype_matrix)) {
      cat(paste0(
        'Info: raw genotype matrix will be created using files in sample_sheet$vcf_file'
      ), file = check_log, append = TRUE, sep = "\n")

      genotype_matrix <- create_genotype_matrix(
        sample_sheet = sample_sheet,
        output_name = "raw_merged",
        output_directory = output_directory,
        chromosomes = chromosomes,
        ...,
        nb_cores = nb_cores
      )
    } else {
      genotype_matrix <- data.table::fread(genotype_matrix, header = TRUE)
    }

    genotype_matrix <- genotype_matrix[, c("CHROM", "POS") := data.table::tstrsplit(var_id, split = "_")[1:2]]
    genotype_matrix <- genotype_matrix[, POS := as.numeric(POS)]

    if (!is.null(chromosomes)) {
      genotype_matrix <- genotype_matrix[CHROM %in% chromosomes]
    }

    samples_available <- samples_ids %in% colnames(genotype_matrix)
    if (!all(samples_available)) {
      cat(paste(
        "Following samples do not exist in genotype matrix and will not be check against the coverage info:",
        paste(samples_ids[!samples_available], collapse = ", ")
      ), file = check_log, append = TRUE, sep = "\n")
    }
    samples_ids <- samples_ids[samples_available]
    genotype_matrix <- genotype_matrix[, .SD, .SDcols = c("var_id", "CHROM", "POS", samples_ids)]

    ## correct geno NA to 0 if position exists in coverage file, else to -1
    cat("Rectifying genotype against coverage info for sample:", file = check_log, append = TRUE, sep = "\n")

    res_i <- lapply(samples_ids, function(iid) { # sample i
      cat(paste0(iid, "\t"), file = check_log, append = TRUE)

      cov_path <- file.path(compress_cov_directory, paste0(iid, ".cov.compress"))

      if (file.size(cov_path) == 0) {
        cat(paste("\nInfo:", cov_path, "is empty."), file = check_log, append = TRUE)
        return(invisible())
      }

      tmp_cov <- data.table::fread(cov_path, header = TRUE)

      res_j <- data.table::rbindlist(parallel::mclapply(
        X = unique(genotype_matrix[["CHROM"]]),
        mc.cores = min(length(unique(tmp_cov[["chr"]])), nb_cores),
        mc_tmp_cov = tmp_cov,
        mc_all_snp_geno = genotype_matrix,
        mc_iid = iid,
        FUN = function(jchr, mc_tmp_cov, mc_all_snp_geno, mc_iid) {
          if (! jchr %in% tmp_cov[["chr"]]) {
            return(NULL)
          }

          tmp_cov_chr <- mc_tmp_cov[chr %in% jchr]
          tmp_ind_chr <- mc_all_snp_geno[CHROM %in% jchr, .SD, .SDcols = c("POS", mc_iid)]

          real_0 <- tmp_ind_chr[, idx := seq(nrow(tmp_ind_chr))][is.na(get(mc_iid))][data.table::inrange(POS, tmp_cov_chr$start, tmp_cov_chr$end)][["idx"]]
          tmp_ind_chr[real_0, (mc_iid) := 0]
          tmp_ind_chr[is.na(get(mc_iid)), (mc_iid) := -1][, !c("idx", "POS")]
        }
      ))
    })
    final_matrix <- cbind(genotype_matrix[, "var_id"], do.call("cbind", res_i))
  })

  cat(
    paste("\nAnalysis completed in", round(analysis_time[["elapsed"]]/60, 2), "minutes"),
    file = check_log, append = TRUE, sep = "\n"
  )

  data.table::fwrite(
    final_matrix, file = output_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  invisible(final_matrix)
}
