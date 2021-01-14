#' compress_coverage
#'
#' Compress coverage file (output of `samtools depth`) into contiguous segments based on position.
#' This function works with whole genome or specified chromosomes.
#' The coverage file should contain columns for chromosome, position and depth, with no header.
#' Attention, currently only supports coverage file in raw text or .gz/.bz2 compressed format.
#'
#' @param sample_sheet A `data.frame`. The data frame should contain two columns: samples `id` and
#'     the path to coverage files `cov_file`.
#' @param output_directory A `character`. The path to the output directory.
#' @param chromosomes A `character` vector. The chromosomes to keep (same nomenclature as chromosome
#'     in coverage file).
#'     If argument not specified (default is `NULL`), this filter will be turned off.
#' @param min_depth A `numeric`. The minimum depth to keep. Default is `8`.
#' @param bin_path A `list(character)`. A list giving the binary path of `awk` and `bgzip` (or `bgzip2`)
#'     in case of compressed coverage file.
#' @param nb_cores A `numeric`. The number of CPUs to use. Default is `1`.
#'
#' @return NULL
#'
#' @importFrom utils tail
#' @importFrom parallel mclapply
#' @importFrom data.table fread
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table fwrite
#'
#' @export
compress_coverage <- function(
  sample_sheet,
  output_directory = NULL,
  chromosomes = NULL,
  min_depth = 8,
  bin_path = list(
    awk = "/usr/bin/awk",
    bgzip = "/usr/bin/bgzip",
    bgzip2 = NULL
  ),
  nb_cores = 1
) {

  if (! all(c("id", "cov_file") %in% colnames(sample_sheet))) {
    stop("Either 'id' or 'cov_file' is missing in the sample_sheet.")
  }

  if (is.null(output_directory)) {
    stop("The output_directory must be specified.")
  }
  dir.create(path = output_directory, showWarnings = FALSE, recursive = TRUE)

  if (is.null(chromosomes)) {
    filter_cond <- paste0("{if ($3 >= ", min_depth, "){print}}")
  } else {
    filter_cond <- paste0(
      "{if ((",
      paste0(paste0("$1 == ", shQuote(chromosomes, type = "cmd")), collapse = " || "),
      ") && $3 >= ", min_depth, "){print}}"
    )
  }

  invisible(parallel::mclapply(
    X = seq_along(sample_sheet[["id"]]),
    mc.cores = min(nb_cores, nrow(sample_sheet)),
    mc_sample_sheet = sample_sheet,
    mc_filter_cond = filter_cond,
    FUN = function(i, mc_sample_sheet, mc_filter_cond) {
      iid <- mc_sample_sheet[["id"]][i]
      log_file <- file.path(output_directory, paste0(iid, ".log"))
      filtered_file <- file.path(output_directory, paste0(iid, ".cov.filtered"))
      filter_script <- file.path(output_directory, paste0("filter-script_", iid))
      cov_file <- mc_sample_sheet[["cov_file"]][i]

      cat(paste0("Analyzing sample ", iid, ":"), file = log_file, sep = "\n")

      if(! file.exists(cov_file)) {
        cat(
          "The coverage file does not exist, sample skipped!\n",
          file = log_file, append = TRUE, sep = "\n"
        )
        return(invisible())
      }

      isgz_cov <- grepl("\\.gz$", cov_file)
      isbz2_cov <- grepl("\\.bz2$", cov_file)

      cat(
        paste0("Filtering position where coverage is lower than ", min_depth, "X..."),
        file = log_file, append = TRUE, sep = "\n"
      )

      filtering_time <- system.time({
        ## decompress if needed and build filtering command
        cat(paste0(
          "#!/bin/bash\n",
          ifelse(
            isgz_cov | isbz2_cov,
            paste0(
              bin_path[["awk"]], " '", mc_filter_cond, "' ",
              "<(", ifelse(test = isgz_cov, yes = bin_path[["bgzip"]], no = bin_path[["bzip2"]]), " -dc ", cov_file, ") ",
              "> ", filtered_file
            ),
            paste0(
              bin_path[["awk"]], " '", mc_filter_cond, "' ",
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
    }
  ))

  invisible(
    lapply(X = sample_sheet[["id"]], FUN = function(iid) {
      log_file <- file.path(output_directory, paste0(iid, ".log"))
      filtered_file <- file.path(output_directory, paste0(iid, ".cov.filtered"))
      compressed_file <- file.path(output_directory, paste0(iid, ".cov.compress"))

      compressing_time <- system.time({
        cat(paste0("Reading ", iid, ".cov.filtered..."), file = log_file, append = TRUE, sep = "\n")

        if (file.size(filtered_file) == 0) {
          cat(
            paste0("No position with coverage greater than ", min_depth),
            file = log_file,
            append = TRUE,
            sep = "\n"
          )
          cat(paste0("Analysis is skipped for ", iid), file = log_file, append = TRUE, sep = "\n")
          return(invisible())
        }

        V1 <- NULL
        cov_file_filtered <- data.table::fread(
          filtered_file, header = FALSE, nThread = nb_cores, showProgress = FALSE
        )
        invisible(file.remove(filtered_file))

        cat(
          paste0("Distinguishing contiguous segments..."),
          file = log_file,
          append = TRUE,
          sep = "\n"
        )
        res <- data.table::rbindlist(lapply(unique(cov_file_filtered[["V1"]]), function(jchr) {
          cat(paste0("in chromosome ", jchr), file = log_file, append = TRUE, sep = "\n")

          sub_file <- cov_file_filtered[V1 %in% jchr]
          tmp <- split(
            x = sub_file[[ncol(sub_file) - 1]],
            f = cumsum(c(1, diff(sub_file[[ncol(sub_file) - 1]]) != 1))
          )

          res_chr <- data.table::data.table(
            "chr" = jchr,
            "start" = sapply(tmp, `[[`, 1L),
            "end" = sapply(tmp, utils::tail, 1L)
          )
          return(res_chr)
        }))

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
