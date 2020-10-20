#' Convert VCF to different assembly with CrossMap
#'
#' @param input_directory  A `character`. The path to the VCFs files or directory. Default is `NULL`.
#' @param output_directory A `character`. The path to the output directory. Default is `NULL`.
#' @param ref_fasta A `character`. A path to the reference fasta file. Default is `NULL`.
#' @param chain_file A `character`. A path to the chain file. Default is `NULL`.
#' @param bin_path A `list(character)`. A list giving the binary path of
#'     `CrossMap`, `bcftools`, `tabix` and `bgzip`.
#' @param nb_cores An `integer`. The number of CPUs to use.
#'
#' @return NULL
#'
#' @export
convert_assembly <- function(
  input_directory = NULL,
  output_directory = NULL,
  ref_fasta = NULL,
  chain_file = NULL,
  bin_path = list(
    CrossMap = "/usr/local/bin/CrossMap.py",
    bcftools = "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
  ),
  nb_cores = 1
) {

  if (length(input_directory) == 1 & length(list.files(input_directory, pattern = "(.vcf.gz|.vcf)$")) <= 1) {
    stop('"input_directory" must be a directory of VCF files splitted by chromosome!')
  }

  if (
    length(input_directory) > 1 &
      !all(file.exists(input_directory)) &
      !all(grepl("(.vcf.gz|.vcf)$", input_directory))
  ) {
    stop('"input_directory" must be a vector of VCF files splitted by chromosome!')
  }

  if (length(input_directory) > 1) {
    vcf_files <- input_directory
  } else {
    vcf_files <- list.files(input_directory, pattern = "(.vcf.gz|.vcf)$", full.names = TRUE)
  }

  temp_directory <- file.path(tempdir(), "convert_assembly")
  dir.create(path = temp_directory, showWarnings = FALSE)

  ## convert to target assembly
  invisible(parallel::mclapply(
    X = vcf_files,
    mc.cores = min(length(vcf_files), nb_cores),
    FUN = function(file_i) {
      output_file <- file.path(temp_directory, basename(file_i))
      system(
        intern = TRUE, wait = TRUE,
        command = paste(bin_path[["CrossMap"]], "vcf", chain_file, file_i, ref_fasta, output_file)
      )
      system(
        intern = TRUE, wait = TRUE,
        command = paste(bin_path[["bcftools"]], "sort", "-Oz -o", output_file, output_file)
      )
      system(intern = TRUE, wait = TRUE, command = paste(bin_path[["tabix"]], "-f -p vcf", output_file))
    }
  ))

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
      unlink(c(paste0(vcf_i, ".tbi"), vcf_to_bind))
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

}
