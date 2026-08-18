# Ad-hoc diagnostic: scan every file in an already-cached picked_peaks.rds
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
#   Rscript scripts/check_empty_spectra.R <path-to-picked_peaks.rds>
#
# Example:
#   Rscript scripts/check_empty_spectra.R /data/output/RP_NEG/B18W40/picked_peaks.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/check_empty_spectra.R <picked_peaks.rds path>", call. = FALSE)
}
cache_path <- args[1]

library(xcms)

cache <- readRDS(cache_path)
xdata <- cache$xdata
filepaths <- cache$filepaths
n_files <- length(filepaths)

message(sprintf("Checking %d file(s) for completely empty spectra...", n_files))

n_peaks_per_spectrum <- MSnbase::peaksCount(xdata)
file_idx <- xcms::fromFile(xdata)
rt <- xcms::rtime(xdata)

empty_idx <- which(n_peaks_per_spectrum == 0)

if (length(empty_idx) == 0) {
  cat("No completely empty spectra found -- this bug's specific trigger isn't present in this file set.\n")
  quit(status = 0)
}

empty_files <- file_idx[empty_idx]
total_spectra_per_file <- tabulate(file_idx, nbins = n_files)

# Ordinal position of each empty spectrum WITHIN its own file's spectrum
# sequence -- the bug's classic trigger is a TRAILING empty spectrum (at or
# near the end of a file's own run), so this is the most useful single
# number for eyeballing whether that's what's going on here.
ordinal_in_file <- vapply(seq_along(empty_idx), function(k) {
  f <- empty_files[k]
  file_spectrum_positions <- which(file_idx == f)
  match(empty_idx[k], file_spectrum_positions)
}, integer(1))

detail <- data.frame(
  file = basename(filepaths[empty_files]),
  spectrum_ordinal = ordinal_in_file,
  total_spectra_in_file = total_spectra_per_file[empty_files],
  is_trailing = ordinal_in_file >= total_spectra_per_file[empty_files] - 2, # last 3 scans
  rt = rt[empty_idx],
  stringsAsFactors = FALSE
)

cat(sprintf(
  "\n%d empty spectrum/spectra found, across %d of %d file(s):\n\n",
  length(empty_idx), length(unique(empty_files)), n_files
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
