# Ad-hoc diagnostic: inspect one file's raw peak count and retention-time
# range from an already-cached picked_peaks.rds, to check whether it's a
# degenerate/near-empty injection (few peaks, abnormal RT range) -- the
# kind of file that can crash adjustRtime()'s obiwarp step (profile-matrix
# column-range slicing assumes every file has a normal, overlapping range).
#
# Usage:
#   Rscript scripts/check_file_health.R <path-to-picked_peaks.rds> <filename-pattern>
#
# Example:
#   Rscript scripts/check_file_health.R \
#     /data/output/RP_NEG/B18W40/picked_peaks.rds \
#     LV7005308613_339

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/check_file_health.R <picked_peaks.rds path> <filename pattern>", call. = FALSE)
}
cache_path <- args[1]
pattern <- args[2]

library(xcms)

cache <- readRDS(cache_path)
xdata <- cache$xdata
filepaths <- cache$filepaths

n_files <- length(filepaths)
raw_peak_count <- tabulate(xcms::chromPeaks(xdata)[, "sample"], nbins = n_files)

rt <- xcms::rtime(xdata)
file_idx <- xcms::fromFile(xdata)
rt_min <- vapply(seq_len(n_files), function(i) {
  vals <- rt[file_idx == i]
  if (length(vals) == 0) NA_real_ else min(vals)
}, numeric(1))
rt_max <- vapply(seq_len(n_files), function(i) {
  vals <- rt[file_idx == i]
  if (length(vals) == 0) NA_real_ else max(vals)
}, numeric(1))
n_spectra <- tabulate(file_idx, nbins = n_files)

summary_df <- data.frame(
  file = basename(filepaths),
  raw_peak_count = raw_peak_count,
  n_spectra = n_spectra,
  rt_min = rt_min,
  rt_max = rt_max,
  rt_range = rt_max - rt_min,
  stringsAsFactors = FALSE
)

target_idx <- grep(pattern, summary_df$file, fixed = TRUE)
if (length(target_idx) == 0) {
  stop(sprintf("No file matching pattern \"%s\" found in %s", pattern, cache_path), call. = FALSE)
}

cat("=== Target file(s) ===\n")
print(summary_df[target_idx, ])

cat("\n=== Batch-wide summary (for comparison) ===\n")
print(summary(summary_df[, c("raw_peak_count", "n_spectra", "rt_min", "rt_max", "rt_range")]))

cat("\n=== Files with suspiciously low peak count or spectra count (< 10% of median) ===\n")
low_peaks <- summary_df$raw_peak_count < 0.1 * median(summary_df$raw_peak_count)
low_spectra <- summary_df$n_spectra < 0.1 * median(summary_df$n_spectra)
suspicious <- summary_df[low_peaks | low_spectra, ]
if (nrow(suspicious) == 0) {
  cat("(none)\n")
} else {
  print(suspicious)
}
