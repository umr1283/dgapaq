#' Convert VCF to different assembly with CrossMap
#'
#' @param ref_fasta A `character`. A path to the reference fasta file.
#' @param chain_file A `character`. A path to the chain file.
#' @param input_directory  A `character`. The path to the VCFs, files should be prefixed by
#'     chromosome name like: 1, 2, ..., 23 (X), 24 (Y), 25 (M or MT).
#' @param output_directory A `character`. The path to the output directory.
#' @param bin_path A `list(character)`. A list giving the binary path of `CrossMap`, `bcftools`,
#'     `tabix` and `bgzip`
#' @param nb_cores An `integer`. The number of CPUs to use.
#'
#' @return NULL
#'
#' @export
convert_assembly <- function(
  ref_fasta,
  chain_file,
  input_directory ,
  output_directory,
  bin_path = list(
    CrossMap = "/usr/local/bin/CrossMap.py",
    bcftools = "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
  ),
  nb_cores = 1
) {

  if (all(utils::file_test("-f", input_directory))) {
    vcf_files <- input_directory
  } else if (utils::file_test("-d", input_directory)) {
    vcf_files <- list.files(input_directory , pattern = "(.vcf.gz|.vcf)$", full.names = TRUE)
  } else {
    stop('"input_directory" contains file or directory which does not exist, please check!')
  }

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  }

  ## convert to target assembly
  invisible(parallel::mclapply(vcf_files, mc.cores = min(length(vcf_files), nb_cores), function(file_i) {
    out_name <- file.path(output_directory, basename(file_i))
    system(paste(
      bin_path[["CrossMap"]],
      "vcf",
      chain_file,
      file_i,
      ref_fasta,
      out_name
    ))
    system(paste(
      bin_path[["bcftools"]],
      "sort",
      "-Oz -o", out_name,
      out_name
    ))
    system(paste(bin_path[["tabix"]], "-f -p vcf", out_name))
  }))

  ## resolve changed chromosomes
  chrs <- c(1:25, "X", "Y", "M", "MT")
  chrs <- c(chrs, paste0("chr", chrs))
  list_vcfs <- list.files(output_directory, pattern = "vcf.gz$", full.names = TRUE)
  invisible(lapply(list_vcfs, function(vcf_i) {
    chr_i_orig <- gsub("([0-9]+|X|Y|M|MT).*", "\\1", basename(vcf_i))
    chr_list <- system(paste(
      bin_path[["tabix"]],
      "--list-chroms",
      vcf_i
    ), intern = TRUE)
    chr_list <- intersect(chr_list, chrs)

    parallel::mclapply(chr_list, mc.cores = min(length(chr_list), nb_cores), function(chr_i) {
      out_name <- file.path(
        dirname(vcf_i),
        gsub(chr_i_orig, paste0(chr_i, "_from_", chr_i_orig), basename(vcf_i))
      )
      system(paste(
        bin_path[["bcftools"]],
        "view",
        "--regions", chr_i,
        vcf_i,
        "-Oz -o", out_name
      ))
      system(paste(bin_path[["tabix"]], "-f -p vcf", out_name))
    })
  }))

  list_vcfs_split <- list.files(output_directory, pattern = "from.*vcf.gz$", full.names = TRUE)
  invisible(parallel::mclapply(
    list_vcfs, mc.cores = min(length(list_vcfs), nb_cores), mc.preschedule = FALSE, function(vcf_i) {
      chr_i <- gsub("([0-9]+|X|Y|M|MT).*", "\\1", basename(vcf_i))
      vcf_to_bind <- list_vcfs_split[grep(paste0("^", chr_i, "_from"), basename(list_vcfs_split))]
      vcf_tmp <- file.path(dirname(vcf_i), paste0("tempo_", basename(vcf_i)))
      system(paste(
        bin_path[["bcftools"]],
        "concat -a", paste(vcf_to_bind, collapse = " "),
        "-Oz -o", vcf_tmp
      ))
      unlink(paste0(vcf_i, ".tbi"))
      system(paste(
        bin_path[["bcftools"]],
        "sort",
        "-Oz -o", vcf_i,
        vcf_tmp
      ))
      unlink(vcf_tmp)
      system(paste(bin_path[["tabix"]], "-p vcf", vcf_i))
    }
  ))

  unlink(list.files(output_directory, pattern = "_from_", full.names = TRUE))
}
