#' compress_coverage
#'
#' Compresses coverage file into contiguous segments based on position,
#' this function works whole genome or speficied chromosome.
#' Attention, currently, the only supported raw coverage formats are .cov, .cov.gz or .cov.bz2.
#' Attention, if the *.cov.compress exists already, the function will skip to next sample.
#'
#' @importFrom parallel mclapply
#' @importFrom data.table fread
#' @importFrom R.utils seqToIntervals
#' @param sample_sheet A `data.frame` or a `character`. a data frame or the full path to sample sheet (.xlsx or .txt),
#'  which contains at least 2 columns: sample 'id' and sequencing 'run' name
#' @param output_directory A `character`. The directory path where to save compressed coverage files *.cov.compress
#' @param chromosome A vector of `character`. Chromosome names to filter
#'  (same nomenclature as chromosome name in coverage file, ex: "chr1", "chrX"),
#' by default filter on all chromosomes present in coverage file.
#' @param min_depth A `numeric`. The minimum depth to filter on each position, default is 8
#' @param nb_cores An `integer`. The number of CPUs to use.
#' @param ... optional arguments to function 'create_sample_sheet',
#' by default using create_cov=TRUE, create_vcf=TRUE and type_vcf="snps"
#'
#' @return NULL
#' @example
#'
#' # compress_coverage(
#   # sample_sheet = "/disks/DATATMP/SB_lning/test_compressCov.txt",
#   sample_sheet = imported_sheet,
#   output_directory = "/disks/DATATMP/SB_lning/test_cov/",
#   create_vcf = TRUE,
#   create_cov = TRUE,
#   type_vcf = "indels",
#   nb_core = 2
# )

compress_coverage <- function(
  sample_sheet,
  output_directory,
  chromosome,
  min_depth = 8,
  nb_cores = 1
) {

  if (!requireNamespace(c("parallel", "data.table", "utils", "R.utils"))) {
    stop('Packages "parallel", "data.table", "utils" and "R.utils" must be installed')
  }

  if (missing(output_directory)) {
    stop("The output_directory must be specified (path to store the compressed cov files).")
  }

  if (!grepl("/$", output_directory)) {
    output_directory <- paste0(output_directory, "/")
  }

  dir.create(path = output_directory, showWarnings = FALSE, recursive = TRUE)
  sample_sheet <- create_sample_sheet(sample_sheet = sample_sheet, ...)

  # temp <- file.path(output_directory, "tmp")
  # dir.create(file.path(temp, "1kgenome_pruned"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
  #

  if (missing(chromosome)) {
    filter_cond <- paste0("{if ($var >= ", min_depth, "){print}}")
  } else {
    filter_cond <- paste0(
      "{if ((",
      paste0(paste0("$1 == ", shQuote(chromosome, type = "cmd")), collapse = " || "),
      ") && $var >= ", min_depth, "){print}}"
    )
  }

  invisible(
    parallel::mclapply(sample_sheet$id, mc.cores = nb_core, function(iid) {
      log_file <- paste0(output_directory, iid, ".log")
      filtered_file <- paste0(output_directory, iid, ".cov.filtered")
      filter_script <- paste0(output_directory, "filter-script_", iid)
      compressed_file <- paste0(output_directory, iid, ".cov.compress")
      cov_file <- sample_sheet$cov_file[which(sample_sheet$id %in% iid)]
      isgz_cov <- grepl("\\.cov\\.gz$", cov_file)
      isbz2_cov <- grepl("\\.cov\\.bz2$", cov_file)
      isjust_cov <- grepl("\\.cov$", cov_file)

      if(!isjust_cov & !isgz_cov & !isbz2_cov) {
        message(paste0("WARNING : this cover file ",iid," do not end neither by .cov nor .cov.gz or .cov.bz2.\n"))
        message(cov_file)
        stop(paste0("WARNING : this cover file ",iid," do not end neither by .cov nor .cov.gz or .cov.bz2.\n"))
      }

      cat(paste0("Analyzing sample ", iid, ":"), file = log_file, sep = "\n")
      cat(
        paste0("Filtering position where coverage is lower than ", min_depth, "X..."),
        file = log_file,
        append = TRUE,
        sep = "\n"
      )

      ## 1/ decompress if needed and build filtering cmd
      filtering_time <- system.time({
        cat(paste0(
          "#!/bin/bash\n",
          ifelse(
            isgz_cov | isbz2_cov,
            paste0(
              "n_col=\"$(head -c100 ", cov_file, " | zcat 2>/dev/null | head -n1 | wc -w)\"", "\n",
              "awk -v var=$n_col '", filter_cond, "' ",
              "<(", ifelse(test = isgz_cov, yes = "gzip", no = "bzip2"), " -dc ", cov_file, ") ",
              "> ", filtered_file
            ),
            paste0(
              "n_col=\"$(head -n1 ", cov_file, "| wc -w)\"", "\n",
              "awk -v var=$n_col '", filter_cond, "' ",
              cov_file, "> ", filtered_file
            )
          )
        ), file = filter_script)

        ## 2/ exe cmd
        system(
          paste0("bash ", filter_script, " && rm ", filter_script),
          intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = TRUE
        )
      })
      cat(
        paste("Filtering time:", round(filtering_time[["elapsed"]]/60, 2), "minutes"),
        file = log_file,
        append = TRUE,
        sep = "\n"
      )

      compressing_time <- system.time({
        cat(paste0("Reading ", iid, ".cov.filtered..."), file = log_file, append = TRUE, sep = "\n")

        if (file.size(filtered_file) == 0) {
          cat(paste0("No position with coverage greater than", min_depth, "X!"), file = log_file, append = TRUE, sep = "\n")
          cat(paste0("Analysis is skipped for ", iid), file = log_file, append = TRUE, sep = "\n")
          return(invisible())
        }

        cov_file_filtered <- data.table::fread(filtered_file, data.table = FALSE, header = FALSE)
        trash <- file.remove(filtered_file)

        cat(paste0("Distinguishing contiguous segments..."), file = log_file, append = TRUE, sep = "\n")
        res <- lapply(unique(cov_file_filtered[["V1"]]), function(jchr) {
          cat(paste0("in chromosome ", jchr), file = log_file, append = TRUE, sep = "\n")

          sub_file <- cov_file_filtered[cov_file_filtered[["V1"]] %in% jchr, ]
          tmp <- split(
            x = sub_file[, ncol(sub_file) - 1],
            f = cumsum(c(1, diff(sub_file[, ncol(sub_file) - 1]) != 1))
          )
          res_intervals <- do.call("rbind", lapply(tmp, R.utils::seqToIntervals))
          cbind(jchr, res_intervals)
        })

        res_ok <- do.call(rbind, res)
        colnames(res_ok) <- c("chr", "start", "end")

        ref_chr <- c(1:26, 23:26)
        names(ref_chr) <- c(1:26, "X", "Y", "M", "MT")
        res_ok[, "chr"] <- ref_chr[gsub("chr", "", res_ok[, "chr"], ignore.case = TRUE)]

        utils::write.table(res_ok, file = compressed_file, sep = "\t", row.names = FALSE, quote = FALSE)

        cat(
          paste0("Result \"", iid, ".cov.compress\" is stocked in: ", output_directory),
          file = log_file, append = TRUE, sep = "\n"
        )
      })
      cat(
        paste("Compressing time:", round(compressing_time[["elapsed"]]/60, 2), "minutes"),
        file = log_file,
        append = TRUE,
        sep = "\n"
      )
    })
  )
}
