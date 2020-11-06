#' compress_coverage
#'
#' Compresses coverage file into contiguous segments based on position,
#' this function works with whole genome or any specified chromosomes.
#' Attention, currently, the only supported raw coverage formats are .cov, .cov.gz or .cov.bz2.
#'
#' @importFrom parallel mclapply
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @param sample_sheet A `data.frame`. A data frame containing sample `id` and the full path to .cov files `cov_file`.
#' @param output_directory A `character`. The path to the output directory.
#' @param chromosome A vector of `character`. Chromosome names to filter
#'  (same nomenclature as chromosome name in coverage file),
#' by default filter on all chromosomes present in coverage file.
#' @param min_depth A `numeric`. The minimum depth to filter. Default is `8`.
#' @param nb_cores An `integer`. The number of CPUs to use. Default is `1`.
#'
#' @return NULL

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
    stop("The output_directory must be specified (path to store the compressed cov files).")
  }

  if (! dir.exists(output_directory)) {
    dir.create(path = normalizePath(output_directory), showWarnings = FALSE, recursive = TRUE)
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
      isgz_cov <- grepl("\\.cov\\.gz$", cov_file)
      isbz2_cov <- grepl("\\.cov\\.bz2$", cov_file)
      isjust_cov <- grepl("\\.cov$", cov_file)

      cat(paste0("Analyzing sample ", iid, ":"), file = log_file, sep = "\n")

      if(! any(isjust_cov, isgz_cov, isbz2_cov)) {
        cat(
          "The coverage file is neither '.cov' nor '.cov.gz' nor '.cov.bz2', sample skipped!\n",
          file = log_file,
          append = TRUE,
          sep = "\n"
        )
        return(invisible())
      }

      cat(
        paste0("Filtering position where coverage is lower than ", min_depth, "X..."),
        file = log_file,
        append = TRUE,
        sep = "\n"
      )

      filtering_time <- system.time({
        ## 1/ decompress if needed and build filtering cmd
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
          cat(
            paste0("No position with coverage greater than ", min_depth, "X!"),
            file = log_file,
            append = TRUE,
            sep = "\n"
          )
          cat(paste0("Analysis is skipped for ", iid), file = log_file, append = TRUE, sep = "\n")
          return(invisible())
        }

        cov_file_filtered <- data.table::fread(filtered_file, header = FALSE)
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
          ref_chr <- c(1:26, 23:26)
          names(ref_chr) <- c(1:26, "X", "Y", "M", "MT")

          cbind(
            as.numeric(ref_chr[gsub("chr", "", jchr, ignore.case = TRUE)]),
            sapply(tmp, `[[`, 1),
            sapply(tmp, tail, 1)
          )
        }))
        colnames(res) <- c("chr", "start", "end")

        invisible(data.table::fwrite(
          res, file = compressed_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
        ))

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
