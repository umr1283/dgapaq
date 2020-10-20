#' Convert VCF to different assembly with CrossMap
#'
#' @param input_directory  A `character`. The path to the VCFs files or directory. Default is `NULL`.
#' @param output_directory A `character`. The path to the output directory. Default is `NULL`.
#' @param ref_fasta A `character`. A path to the reference fasta file. Default is `NULL`.
#' @param chain_file A `character`. A path to the chain file. Default is `NULL`.
#' @param bin_path A `list(character)`. A list giving the binary path of `CrossMap`, `bcftools`, `tabix` and `bgzip`.
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
    crossmap = "/usr/local/bin/CrossMap.py",
    bcftools = "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
  ),
  nb_cores = 1
) {

  names(bin_path) <- tolower(names(bin_path))

  if (
    length(intersect(names(bin_path), c("crossmap", "bcftools", "tabix", "bgzip"))) != 4 |
      !all(sapply(bin_path, file.exists))
  ) {
    stop('"bin_path" must contains valid named path to "crossmap", "bcftools", "tabix" and "bgzip"!')
  }

  if (
    length(input_directory) == 1 &
      length(list.files(input_directory, pattern = "(.vcf.gz|.vcf)$")) <= 1
  ) {
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
        command = paste(bin_path[["crossmap"]], "vcf", chain_file, file_i, ref_fasta, output_file),
        wait = TRUE
      )
      system(
        command = paste(bin_path[["bcftools"]], "sort", "-Oz -o", output_file, output_file),
        wait = TRUE
      )
      system(
        command = paste(bin_path[["tabix"]], "-f -p vcf", output_file),
        wait = TRUE
      )
    }
  ))

  ## resolve changed chromosomes
  chrs <- paste0(rep(c("", "chr"), each = 29), c(1:25, "X", "Y", "M", "MT"))
  list_vcfs <- list.files(temp_directory, pattern = "vcf.gz$", full.names = TRUE)
  invisible(lapply(list_vcfs, function(vcf_i) {
    chr_i_orig <- gsub(paste0("^(", paste(chrs, collapse = "|"), ").*"), "\\1", basename(vcf_i))
    chr_list <- system(
      command = paste(bin_path[["tabix"]], "--list-chroms", vcf_i),
      intern = TRUE,
      wait = TRUE
    )
    chr_list <- intersect(chr_list, chrs)

    parallel::mclapply(
      X = chr_list,
      mc.cores = min(length(chr_list), nb_cores),
      mc.preschedule = FALSE,
      FUN = function(chr_i) {
        output_file <- file.path(
          temp_directory,
          gsub(chr_i_orig, paste0(chr_i, "_from_", chr_i_orig), basename(vcf_i))
        )
        system(
          command = paste(bin_path[["bcftools"]], "view", "--regions", chr_i, vcf_i, "-Oz -o", output_file),
          wait = TRUE
        )
        system(
          command = paste(bin_path[["tabix"]], "-f -p vcf", output_file),
          wait = TRUE
        )
      }
    )
  }))

  invisible(parallel::mclapply(
    X = list_vcfs,
    mc.cores = min(length(list_vcfs), nb_cores),
    mc.preschedule = FALSE,
    FUN = function(vcf_i) {
      chr_i <- gsub(paste0("^(", paste(chrs, collapse = "|"), ").*"), "\\1", basename(vcf_i))
      vcf_to_bind <- paste(
        list.files(temp_directory, pattern = paste0(chr_i, "_from_.*vcf.gz$"), full.names = TRUE),
        collapse = " "
      )

      vcf_tmp <- file.path(temp_directory, paste0("temp_", basename(vcf_i)))
      vcf_out <- file.path(output_directory, basename(vcf_i))

      system(
        command = paste(bin_path[["bcftools"]], "concat -a", vcf_to_bind, "-Oz -o", vcf_tmp),
        wait = TRUE
      )
      system(
        command = paste(bin_path[["bcftools"]], "sort", "-Oz -o", vcf_out, vcf_tmp),
        wait = TRUE
      )
      system(
        command = paste(bin_path[["tabix"]], "-f -p  vcf", vcf_out),
        wait = TRUE
      )
    }
  ))

  invisible(unlink(temp_directory, recursive = TRUE, force = TRUE))
}
