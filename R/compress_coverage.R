#' compress_coverage
#'
#' Compresses coverage file (output of `samtools depth`) into contiguous segments based on position.
#' This function works with whole genome or any specified chromosomes.
#' Each sample is supposed to have one coverage file.
#' The coverage file should not have header and contain 3 mandatory columns for chromosome, position
#' and depth, with no header.
#' Please note this function is designed for Linux environment, `awk` and `gzip` (or `bzip2`) are
#' required to process data.
#' Attention, currently only supports coverage file in raw text or .gz/.bz2 compressed format.
#'
#' @param sample_sheet A `data.frame`. A data frame containing samples `id` and
#' the full path to coverage files `cov_file`.
#' @param output_directory A `character`. The path to the output directory.
#' @param chromosome A vector of `character`. Chromosomes to keep
#'  (same nomenclature as chromosome name in coverage file).
#' By default use all chromosomes present in the coverage file.
#' @param min_depth An `integer`. The minimum depth to keep. Default is `8`.
#' @param nb_cores An `integer`. The number of CPUs to use. Default is `1`.
#'
#' @return NULL
#'
#' @importFrom parallel mclapply
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
compress_coverage <- function(
  sample_sheet,
  output_directory,
  chromosome,
  min_depth = 8,
  nb_cores = 1
) {

  if (! all(c("id", "cov_file") %in% colnames(sample_sheet))) {
    stop("Either 'id' or 'cov_file' is missing in the sample_sheet.")
  }

  if (missing(output_directory)) {
    stop("The output_directory must be specified.")
  }

  if (! dir.exists(output_directory)) {
    dir.create(path = output_directory, showWarnings = FALSE, recursive = TRUE)
  }

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
    parallel::mclapply(sample_sheet$id, mc.cores = nb_cores, function(iid) {
      log_file <- file.path(output_directory, paste0(iid, ".log"))
      filtered_file <- file.path(output_directory, paste0(iid, ".cov.filtered"))
      filter_script <- file.path(output_directory, paste0("filter-script_", iid))
      compressed_file <- file.path(output_directory, paste0(iid, ".cov.compress"))
      cov_file <- sample_sheet$cov_file[which(sample_sheet$id %in% iid)]

      cat(paste0("Analyzing sample ", iid, ":"), file = log_file, sep = "\n")

      if(! file.exists(cov_file)) {
        cat(
          "The coverage file does not exist, sample skipped!\n",
          file = log_file,
          append = TRUE,
          sep = "\n"
        )
        return(invisible())
      }

      isgz_cov <- grepl("\\.gz$", cov_file)
      isbz2_cov <- grepl("\\.bz2$", cov_file)

      cat(
        paste0("Filtering position where coverage is lower than ", min_depth, "X..."),
        file = log_file,
        append = TRUE,
        sep = "\n"
      )

      filtering_time <- system.time({
        ## decompress if needed and build filtering command
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

        ## execute command
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
          cat(
            paste0("No position with coverage greater than ", min_depth, "X!"),
            file = log_file,
            append = TRUE,
            sep = "\n"
          )
          cat(paste0("Analysis is skipped for ", iid), file = log_file, append = TRUE, sep = "\n")
          return(invisible())
        }

        V1 <- NULL
        cov_file_filtered <- data.table::fread(filtered_file, header = FALSE, nThread = nb_cores)
        invisible(file.remove(filtered_file))

        cat(
          paste0("Distinguishing contiguous segments..."),
          file = log_file,
          append = TRUE,
          sep = "\n"
        )
        res <- do.call("rbind", lapply(unique(cov_file_filtered[["V1"]]), function(jchr) {
          cat(paste0("in chromosome ", jchr), file = log_file, append = TRUE, sep = "\n")

          sub_file <- cov_file_filtered[V1 %in% jchr, ]
          tmp <- split(
            x = sub_file[[ncol(sub_file) - 1]],
            f = cumsum(c(1, diff(sub_file[[ncol(sub_file) - 1]]) != 1))
          )

          cbind(
            jchr,
            sapply(tmp, `[[`, 1L),
            sapply(tmp, utils::tail, 1L)
          )
        }))
        colnames(res) <- c("chr", "start", "end")

        invisible(data.table::fwrite(
          res, file = compressed_file, sep = "\t", nThread = nb_cores,
          row.names = FALSE, col.names = TRUE, quote = FALSE
        ))

        cat(
          paste0("Result \"", iid, ".cov.compress\" is stored in: ", output_directory),
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
