# Ad-hoc diagnostic: scan one or every picked_peaks.rds under a directory
# for completely empty spectra (zero data points -- not "few peaks", a scan
# with nothing in it at all). This is the actual trigger for a known,
# unfixed xcms bug in profMat()'s .insertColumn() helper (used by
# adjustRtime(ObiwarpParam())): a trailing empty spectrum can make it
# compute a column-insert position past the matrix's current width,
# crashing with "Error in x[, pos[i]:ncol(x)]: subscript out of bounds"
# (see github.com/sneumann/xcms/issues/366, still unfixed on master as of
# this writing). Aggregate stats like raw peak count or RT range -- what
# scripts/check_file_health.R checks -- won't catch this: a file can look
# completely normal overall and still have one bad scan buried in it.
#
# Usage:
#   Rscript scripts/check_empty_spectra.R <path>
#
# <path> is either:
#   - a single picked_peaks.rds file -- prints full per-spectrum detail, or
#   - a directory -- recursively finds every picked_peaks.rds under it
#     (works for both batch scope, output/<group>/<batch>/picked_peaks.rds,
#     and global scope, output/<group>/picked_peaks.rds) and prints one
#     summary line per batch/group plus a final list of which ones need
#     their cache invalidated and re-picked.
#
# Examples:
#   Rscript scripts/check_empty_spectra.R /data/output/RP_NEG/B18W40/picked_peaks.rds
#   Rscript scripts/check_empty_spectra.R /data/output/RP_NEG

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/check_empty_spectra.R <picked_peaks.rds path or a directory containing them>", call. = FALSE)
}
target <- args[1]

library(xcms)

#' Check one picked_peaks.rds for completely empty spectra.
#' @return data.frame, one row per empty spectrum found (0 rows if none):
#'   file, spectrum_ordinal, total_spectra_in_file, is_trailing, rt.
check_cache <- function(cache_path) {
  cache <- readRDS(cache_path)
  xdata <- cache$xdata
  filepaths <- cache$filepaths
  n_files <- length(filepaths)

  n_peaks_per_spectrum <- MSnbase::peaksCount(xdata)
  file_idx <- xcms::fromFile(xdata)
  rt <- xcms::rtime(xdata)

  empty_idx <- which(n_peaks_per_spectrum == 0)
  if (length(empty_idx) == 0) {
    return(data.frame(
      file = character(0), spectrum_ordinal = integer(0), total_spectra_in_file = integer(0),
      is_trailing = logical(0), rt = numeric(0), stringsAsFactors = FALSE
    ))
  }

  empty_files <- file_idx[empty_idx]
  total_spectra_per_file <- tabulate(file_idx, nbins = n_files)

  # Ordinal position of each empty spectrum WITHIN its own file's spectrum
  # sequence -- the bug's classic trigger is a TRAILING empty spectrum (at
  # or near the end of a file's own run), so this is the most useful single
  # number for eyeballing whether that's what's going on.
  ordinal_in_file <- vapply(seq_along(empty_idx), function(k) {
    f <- empty_files[k]
    match(empty_idx[k], which(file_idx == f))
  }, integer(1))

  data.frame(
    file = basename(filepaths[empty_files]),
    spectrum_ordinal = ordinal_in_file,
    total_spectra_in_file = total_spectra_per_file[empty_files],
    is_trailing = ordinal_in_file >= total_spectra_per_file[empty_files] - 2, # last 3 scans
    rt = rt[empty_idx],
    stringsAsFactors = FALSE
  )
}

if (dir.exists(target)) {
  cache_paths <- list.files(target, pattern = "^picked_peaks\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(cache_paths) == 0) {
    stop(sprintf("No picked_peaks.rds found anywhere under: %s", target), call. = FALSE)
  }
  message(sprintf("Found %d picked_peaks.rds file(s) under %s -- checking each...", length(cache_paths), target))

  flagged <- character(0)
  for (cache_path in cache_paths) {
    label <- basename(dirname(cache_path)) # batch scope: <batch>; global scope: <group>
    detail <- tryCatch(check_cache(cache_path), error = function(e) {
      message(sprintf("  %s: FAILED to check (%s)", label, conditionMessage(e)))
      NULL
    })
    if (is.null(detail)) next
    if (nrow(detail) == 0) {
      message(sprintf("  %s: clean (no empty spectra)", label))
    } else {
      n_files_affected <- length(unique(detail$file))
      n_trailing <- sum(detail$is_trailing)
      message(sprintf(
        "  %s: %d empty spectrum/spectra across %d file(s) (%d trailing) -- %s",
        label, nrow(detail), n_files_affected, n_trailing, paste(unique(detail$file), collapse = ", ")
      ))
      flagged <- c(flagged, cache_path)
    }
  }

  cat("\n=== Summary ===\n")
  if (length(flagged) == 0) {
    cat("No picked_peaks.rds needs invalidating -- no empty spectra found anywhere.\n")
  } else {
    cat(sprintf("%d of %d picked_peaks.rds file(s) need to be deleted and re-picked:\n", length(flagged), length(cache_paths)))
    cat(paste(" ", flagged, collapse = "\n"), "\n")
  }
} else if (file.exists(target)) {
  detail <- check_cache(target)
  if (nrow(detail) == 0) {
    cat("No completely empty spectra found -- this bug's specific trigger isn't present in this file set.\n")
  } else {
    cat(sprintf(
      "\n%d empty spectrum/spectra found, across %d file(s):\n\n",
      nrow(detail), length(unique(detail$file))
    ))
    print(detail[order(-detail$is_trailing, detail$file), ])

    cat(sprintf(
      "\n%d of these are trailing (within the last 3 scans of their file) -- the pattern known to trigger the crash.\n",
      sum(detail$is_trailing)
    ))

    cat("\n=== Per-file summary ===\n")
    per_file <- aggregate(spectrum_ordinal ~ file, data = detail, FUN = length)
    names(per_file)[2] <- "n_empty_spectra"
    print(per_file[order(-per_file$n_empty_spectra), ])
  }
} else {
  stop(sprintf("Path not found: %s", target), call. = FALSE)
}
