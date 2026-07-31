# IPO-based optimization of xcms centWave parameters, and peak picking using
# the optimized parameters, run independently per column x polarity group.

#' Pick a small representative subset of files to run IPO optimization on.
#'
#' IPO's design-of-experiments search evaluates many parameter combinations,
#' so it's run on a handful of files rather than the full batch. Pooled QC
#' injections (sQC/ltQC) are preferred when present since they're designed to
#' be representative of the whole sample set; otherwise fall back to a random
#' subset of regular samples.
#'
#' @param group_sheet Sample sheet rows for a single column x polarity group.
#' @param n Number of files to select.
select_ipo_subset <- function(group_sheet, n = 4) {
  qc_rows <- group_sheet[group_sheet$is_qc, , drop = FALSE]

  if (nrow(qc_rows) >= 2) {
    n <- min(n, nrow(qc_rows))
    subset_rows <- qc_rows[seq_len(n), ]
  } else {
    regular_rows <- group_sheet[!group_sheet$is_qc, , drop = FALSE]
    n <- min(n, nrow(regular_rows))
    subset_rows <- regular_rows[sample(nrow(regular_rows), n), ]
  }

  subset_rows$filepath
}

#' Convert an IPO optimization result into an xcms::CentWaveParam.
#'
#' IPO targets the legacy xcmsSet() API; this adapts its `best_settings`
#' parameters to the modern findChromPeaks()/CentWaveParam workflow.
ipo_result_to_centwave_param <- function(ipo_result) {
  p <- ipo_result$best_settings$parameters

  xcms::CentWaveParam(
    ppm = p$ppm,
    peakwidth = p$peakwidth,
    snthresh = p$snthresh,
    prefilter = c(p$prefilter, p$value_of_prefilter),
    mzCenterFun = p$mzCenterFun,
    integrate = p$integrate,
    mzdiff = p$mzdiff,
    fitgauss = p$fitgauss,
    noise = p$noise,
    verboseColumn = p$verbose.column
  )
}

#' Run IPO optimization for one column x polarity group, caching the result
#' to disk so re-running peak picking doesn't re-optimize from scratch.
#'
#' @param group_sheet Sample sheet rows for this group.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @return An xcms::CentWaveParam with the optimized settings.
run_ipo_optimization <- function(group_sheet, out_dir) {
  params_path <- file.path(out_dir, "ipo_params.rds")

  if (file.exists(params_path)) {
    message(sprintf("Using cached IPO params: %s", params_path))
    return(readRDS(params_path))
  }

  subset_files <- select_ipo_subset(group_sheet)
  message(sprintf(
    "Running IPO optimization on %d file(s):\n  %s",
    length(subset_files), paste(subset_files, collapse = "\n  ")
  ))

  starting_params <- IPO::getDefaultXcmsSetStartingParams("centWave")
  ipo_result <- IPO::optimizeXcmsSet(
    files = subset_files,
    params = starting_params,
    nSlaves = 1,
    subdir = NULL
  )

  centwave_param <- ipo_result_to_centwave_param(ipo_result)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(centwave_param, params_path)
  message(sprintf("Saved optimized params to: %s", params_path))

  centwave_param
}

#' Run peak picking for one column x polarity group using the given params.
#'
#' @param group_sheet Sample sheet rows for this group.
#' @param centwave_param xcms::CentWaveParam to use for findChromPeaks().
#' @param out_dir Group's output directory (e.g. output/RP_POS).
run_peak_picking_group <- function(group_sheet, centwave_param, out_dir) {
  peaks_dir <- file.path(out_dir, "peaks")
  dir.create(peaks_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("Reading %d file(s) for peak picking...", nrow(group_sheet)))
  raw_data <- MSnbase::readMSData(group_sheet$filepath, mode = "onDisk")

  message("Running findChromPeaks()...")
  xdata <- xcms::findChromPeaks(raw_data, param = centwave_param)

  saveRDS(xdata, file.path(peaks_dir, "xdata.rds"))

  peak_table <- as.data.frame(xcms::chromPeaks(xdata))
  write.csv(peak_table, file.path(peaks_dir, "peak_table.csv"), row.names = FALSE)

  message(sprintf(
    "Peak picking complete: %d peaks across %d files. Saved to: %s",
    nrow(peak_table), nrow(group_sheet), peaks_dir
  ))

  invisible(xdata)
}
