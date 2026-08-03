# Detecting likely-faulty QC injections (missed injections, empty vials,
# degraded samples) via a quick default-parameters peak-picking pass, so they
# never get picked as the "representative" QC for IPO2 optimization.

#' Compute the value below which a point counts as a low outlier: a simple
#' fraction of the median. Only penalizes values that are lower than
#' expected — a high value never gets flagged, only low TIC/low feature
#' count does, which is the actual failure signature of a missed injection
#' or empty vial. Deliberately simple and always non-negative/well within
#' the plottable range, unlike a MAD-based z-score (which can go negative
#' for skewed data and end up off-chart).
#'
#' @param x Numeric vector.
#' @param min_fraction_of_median Values below this fraction of the median
#'   are flagged (e.g. 0.5 flags anything less than half the median).
#' @return The numeric threshold below which a value is flagged.
low_outlier_threshold <- function(x, min_fraction_of_median = 0.5) {
  min_fraction_of_median * median(x)
}

#' Flag values that fall below their peers' threshold (see
#' `low_outlier_threshold()`).
#'
#' @param x Numeric vector (e.g. TIC or aligned feature count per file).
#' @param min_fraction_of_median Passed to `low_outlier_threshold()`.
flag_low_outliers <- function(x, min_fraction_of_median = 0.5) {
  x < low_outlier_threshold(x, min_fraction_of_median)
}

#' Check QC injection quality within a group via TIC and aligned feature
#' count, both computed with default (unoptimized) centWave parameters —
#' this only needs to separate "clearly fine" from "clearly broken", not
#' produce publication-quality peak picking.
#'
#' A file is flagged if either its total TIC or its aligned feature count
#' (how many of the group's correspondence-grouped features it actually has
#' a peak for) falls well below the rest of the group's — the signature of a
#' missed injection, empty vial, or otherwise degraded run.
#'
#' Includes RT alignment (`adjustRtime()`) before correspondence, since this
#' check can span many batches/months and genuine RT drift could otherwise
#' make a perfectly fine QC look like it's missing features purely from
#' misalignment. Correspondence groups by `sample_type` (sQC vs ltQC are
#' different matrices and shouldn't be lumped into one group). Deliberately
#' does *not* run `fillChromPeaks()` — gap-filling would paper over exactly
#' the weak/missing signal this check exists to catch.
#'
#' @param qc_sheet Sample sheet rows for all QC (sQC/ltQC) files in a group.
#' @param min_fraction_of_median Threshold passed to `flag_low_outliers()`,
#'   used for the feature-count check always, and for the TIC check only
#'   when no instrument-specific absolute threshold is available (see
#'   below).
#' @return A data.frame with one row per file: filepath, tic,
#'   aligned_feature_count, flagged, reason. The actual threshold used for
#'   each metric is attached as attributes `tic_threshold` and
#'   `feature_threshold` (so reporting can draw the line that was actually
#'   used, not recompute a possibly-different one).
check_qc_quality <- function(qc_sheet, min_fraction_of_median = 0.5) {
  message(sprintf("Reading %d QC file(s) for quality check...", nrow(qc_sheet)))
  raw_data <- read_raw_data(qc_sheet$filepath, get_spectrum_modes(qc_sheet))

  tic_by_file <- split(MSnbase::tic(raw_data), MSnbase::fromFile(raw_data))
  tic <- vapply(tic_by_file, sum, numeric(1))
  tic <- tic[order(as.integer(names(tic)))]

  message("Running default-parameters findChromPeaks() for quality check...")
  xdata <- xcms::findChromPeaks(raw_data, param = xcms::CentWaveParam(), BPPARAM = bp_workers())

  message("Running adjustRtime()...")
  xdata <- xcms::adjustRtime(xdata, param = xcms::ObiwarpParam(), BPPARAM = bp_workers())

  xdata <- xcms::groupChromPeaks(
    xdata,
    param = xcms::PeakDensityParam(sampleGroups = qc_sheet$sample_type), BPPARAM = bp_workers()
  )

  feature_vals <- xcms::featureValues(xdata, value = "into")
  aligned_feature_count <- colSums(!is.na(feature_vals))

  # TIC prefers an absolute, instrument-specific floor over a threshold
  # computed from this batch's own distribution — median/MAD both have a
  # 50% breakdown point, so if a majority of this batch's QCs are actually
  # bad, nothing computed from the batch itself can reliably tell good from
  # bad. An absolute floor doesn't have that problem. No equivalent absolute
  # reference exists for aligned feature count (expected count is
  # method/sample-dependent, not just instrument-dependent), so that one
  # stays relative.
  instrument <- group_instrument(qc_sheet)
  instrument_config <- get_instrument_params(instrument, qc_sheet$column[1], qc_sheet$polarity[1])

  if (!is.null(instrument_config) && !is.null(instrument_config$int_threshold)) {
    tic_threshold <- instrument_config$int_threshold
    message(sprintf("Using instrument-specific absolute TIC threshold: %s", tic_threshold))
  } else {
    tic_threshold <- low_outlier_threshold(tic, min_fraction_of_median)
  }
  feature_threshold <- low_outlier_threshold(aligned_feature_count, min_fraction_of_median)

  tic_flag <- tic < tic_threshold
  feature_flag <- aligned_feature_count < feature_threshold

  reason <- ifelse(
    tic_flag & feature_flag, "low TIC and low aligned feature count",
    ifelse(tic_flag, "low TIC", ifelse(feature_flag, "low aligned feature count", NA_character_))
  )

  result <- data.frame(
    filepath = qc_sheet$filepath,
    sample_label = qc_sheet$sample_label,
    tic = tic,
    aligned_feature_count = aligned_feature_count,
    flagged = tic_flag | feature_flag,
    reason = reason,
    stringsAsFactors = FALSE
  )
  attr(result, "tic_threshold") <- tic_threshold
  attr(result, "feature_threshold") <- feature_threshold
  result
}

#' Barplot of one QC quality metric, colored by pass/fail status, with the
#' flagging threshold drawn as a horizontal line and labeled with its actual
#' value. The y-axis is explicitly forced to include the threshold (not just
#' the bar values), so the line is always on-page and visible even if it
#' sits far below the data (e.g. a strict instrument-specific floor) rather
#' than silently falling outside the default plot range.
#'
#' @param values Numeric vector to plot (e.g. tic or aligned_feature_count).
#' @param flagged Logical vector, same length/order as `values`.
#' @param threshold Numeric threshold to draw as a horizontal line.
#' @param labels Character vector of bar labels (e.g. sample_label).
#' @param title Plot title.
#' @param ylab Y-axis label.
plot_qc_metric <- function(values, flagged, threshold, labels, title, ylab) {
  colors <- ifelse(flagged, "firebrick2", "steelblue3")
  y_range <- range(c(values, threshold, 0))
  bar_x <- barplot(
    values, names.arg = labels, col = colors, las = 2, cex.names = 0.7,
    main = title, ylab = ylab, border = NA, ylim = y_range
  )
  abline(h = threshold, col = "gray30", lty = 2, lwd = 2)
  text(
    x = par("usr")[1], y = threshold, labels = sprintf("threshold = %.3g", threshold),
    pos = 3, cex = 0.7, col = "gray30", xpd = TRUE
  )
  legend(
    "topright", legend = c("Passed", "Flagged", "Threshold"),
    fill = c("steelblue3", "firebrick2", NA), border = c(NA, NA, NA),
    lty = c(NA, NA, 2), col = c(NA, NA, "gray30"), bty = "n", cex = 0.8
  )
  invisible(bar_x)
}

#' Generate a PDF report of QC quality check results: one page per group,
#' each with a TIC barplot and an aligned-feature-count barplot (blue =
#' passed, red = flagged), the pass/fail threshold marked on both.
#'
#' @param results_by_group Named list of `check_qc_quality()` result
#'   data.frames, one per group (each carrying `tic_threshold`/
#'   `feature_threshold` attributes — see `check_qc_quality()`).
#' @param out_path Path to write the PDF to.
generate_qc_report <- function(results_by_group, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(out_path, width = 10, height = 6)
  on.exit(grDevices::dev.off())

  for (group_name in names(results_by_group)) {
    result <- results_by_group[[group_name]]

    plot_qc_metric(
      result$tic, result$flagged, attr(result, "tic_threshold"),
      result$sample_label, sprintf("%s — Total ion current (TIC)", group_name), "TIC"
    )
    plot_qc_metric(
      result$aligned_feature_count, result$flagged, attr(result, "feature_threshold"),
      result$sample_label, sprintf("%s — Aligned feature count", group_name), "Feature count"
    )
  }

  message(sprintf("Saved QC quality report to: %s", out_path))
}
