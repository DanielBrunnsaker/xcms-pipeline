# Detecting likely-faulty QC injections (missed injections, empty vials,
# degraded samples) via a quick default-parameters peak-picking pass, so they
# never get picked as the "representative" QC for IPO2 optimization.

#' Compute the value below which a point counts as a low outlier: a simple
#' fraction of the median.
#' - Only penalizes lower-than-expected values (a high value never gets
#'   flagged) — matches the actual failure signature of a missed injection.
#' - Deliberately simple, always non-negative/plottable, unlike a MAD-based
#'   z-score (which can go negative for skewed data and end up off-chart).
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
#' count (default/unoptimized centWave params — just needs to separate
#' "clearly fine" from "clearly broken", not produce publication-quality
#' peak picking). A file is flagged if either metric falls well below the
#' rest of the group's — the signature of a missed injection, empty vial,
#' or degraded run.
#'
#' - Runs `adjustRtime()` before correspondence: this check can span many
#'   batches/months, and genuine RT drift could make a fine QC look like
#'   it's missing features purely from misalignment.
#' - Correspondence groups by `sample_type` (sQC vs ltQC are different
#'   matrices, shouldn't be lumped together).
#' - Deliberately skips `fillChromPeaks()` — gap-filling would paper over
#'   the exact weak/missing signal this check exists to catch.
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

  # TIC prefers an absolute, instrument-specific floor over a
  # batch-distribution threshold: median/MAD both break down once a
  # majority of a batch's QCs are actually bad. No equivalent absolute
  # reference exists for feature count (method/sample-dependent, not just
  # instrument-dependent), so that one stays relative.
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
    batch = qc_sheet$batch,
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

#' Convert a base color + alpha into an explicit "rgba(r,g,b,a)" string —
#' safer than relying on 8-digit hex alpha support across plotly renderers.
to_rgba <- function(color, alpha) {
  rgb <- grDevices::col2rgb(color)
  sprintf("rgba(%d,%d,%d,%.2f)", rgb[1], rgb[2], rgb[3], alpha)
}

#' Interactive plotly bar chart of one QC quality metric.
#' - Colored by batch identity: one trace per batch (fixed hue order, a
#'   qualitative HCL palette), so the legend is clickable to toggle batches.
#' - Flagged files are rendered at reduced alpha instead of a separate
#'   color — passed = opaque, flagged = faded.
#' - Hover shows the exact sample, batch, value, and status.
#' - Threshold drawn as a labeled horizontal line.
#'
#' @param values Numeric vector to plot (e.g. tic or aligned_feature_count).
#' @param flagged Logical vector, same length/order as `values`.
#' @param batches Character vector of batch labels, same length/order.
#' @param threshold Numeric threshold to draw as a horizontal line.
#' @param labels Character vector of bar labels (e.g. sample_label).
#' @param title Plot title.
#' @param ylab Y-axis label.
#' @return A plotly htmlwidget.
plot_qc_metric <- function(values, flagged, batches, threshold, labels, title, ylab) {
  unique_batches <- sort(unique(batches))
  batch_colors <- stats::setNames(
    grDevices::hcl.colors(length(unique_batches), palette = "Dynamic"), unique_batches
  )

  df <- data.frame(
    label = labels, value = values, batch = batches, flagged = flagged,
    stringsAsFactors = FALSE
  )
  df$color <- mapply(function(b, f) to_rgba(batch_colors[[b]], if (f) 0.3 else 1), df$batch, df$flagged)
  df$status <- ifelse(df$flagged, "Flagged", "Passed")
  df$hover <- sprintf(
    "%s<br>Batch: %s<br>%s: %.4g<br>Status: %s", df$label, df$batch, ylab, df$value, df$status
  )

  p <- plotly::plot_ly()
  for (b in unique_batches) {
    batch_df <- df[df$batch == b, ]
    p <- plotly::add_trace(
      p, data = batch_df, x = ~label, y = ~value, type = "bar", name = b,
      marker = list(color = ~color, line = list(width = 0)),
      text = ~hover, hoverinfo = "text"
    )
  }

  plotly::layout(
    p,
    title = title,
    xaxis = list(title = "", tickangle = -45, categoryorder = "array", categoryarray = labels),
    yaxis = list(title = ylab),
    barmode = "overlay",
    shapes = list(list(
      type = "line", x0 = 0, x1 = 1, xref = "paper",
      y0 = threshold, y1 = threshold, line = list(dash = "dash", color = "gray")
    )),
    annotations = list(list(
      x = 1, y = threshold, xref = "paper", yref = "y",
      text = sprintf("threshold = %.3g", threshold),
      showarrow = FALSE, xanchor = "right", yanchor = "bottom",
      font = list(size = 10, color = "gray")
    ))
  )
}

#' Generate an interactive HTML report of QC quality check results: one
#' TIC chart and one aligned-feature-count chart per group, colored by
#' batch (faded = flagged), pass/fail threshold marked on both, hover for
#' exact values.
#'
#' @param results_by_group Named list of `check_qc_quality()` result
#'   data.frames, one per group (each carrying `tic_threshold`/
#'   `feature_threshold` attributes — see `check_qc_quality()`).
#' @param out_path Path to write the HTML report to.
generate_qc_report <- function(results_by_group, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  plots <- list()
  for (group_name in names(results_by_group)) {
    result <- results_by_group[[group_name]]

    plots[[length(plots) + 1]] <- plot_qc_metric(
      result$tic, result$flagged, result$batch, attr(result, "tic_threshold"),
      result$sample_label, sprintf("%s — Total ion current (TIC)", group_name), "TIC"
    )
    plots[[length(plots) + 1]] <- plot_qc_metric(
      result$aligned_feature_count, result$flagged, result$batch, attr(result, "feature_threshold"),
      result$sample_label, sprintf("%s — Aligned feature count", group_name), "Feature count"
    )
  }

  htmltools::save_html(htmltools::tagList(plots), out_path)
  message(sprintf("Saved QC quality report to: %s", out_path))
}
