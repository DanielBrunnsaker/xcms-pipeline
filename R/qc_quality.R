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

#' Check QC injection quality via TIC and aligned feature count (default
#' centWave params -- just needs to separate "clearly fine" from "clearly
#' broken"). Flags a file if either metric falls well below its peers -- the
#' signature of a missed injection, empty vial, or degraded run.
#'
#' - Runs `adjustRtime()` first: RT drift across batches/months could make a
#'   fine QC look like it's missing features from misalignment alone.
#' - Groups by `sample_type` (sQC/ltQC are different matrices).
#' - Skips `fillChromPeaks()` -- gap-filling would paper over the weak
#'   signal this check exists to catch.
#'
#' @param qc_sheet Sample sheet rows for all QC (sQC/ltQC) files in a group.
#' @param min_fraction_of_median Threshold for `flag_low_outliers()` --
#'   always used for feature count, and for TIC when no instrument-specific
#'   absolute threshold exists.
#' @param min_batch_qc Minimum QC files for a batch to get its own local
#'   threshold (below this, relies on the global check alone).
#' @return A data.frame with one row per file: filepath, sample_label, batch,
#'   injection_order, tic, aligned_feature_count, flagged_global,
#'   flagged_batch, flagged, reason. `tic_threshold`/`feature_threshold`
#'   attributes carry the actual thresholds used (so reporting doesn't need
#'   to recompute them).
check_qc_quality <- function(qc_sheet, min_fraction_of_median = 0.5, min_batch_qc = 2) {
  message(sprintf("Reading %d QC file(s) for quality check...", nrow(qc_sheet)))
  raw_data <- read_raw_data(qc_sheet$filepath, get_spectrum_modes(qc_sheet))

  tic_by_file <- split(MSnbase::tic(raw_data), MSnbase::fromFile(raw_data))
  tic <- vapply(tic_by_file, sum, numeric(1))
  tic <- tic[order(as.integer(names(tic)))]

  message("Running default-parameters findChromPeaks() for quality check...")
  xdata <- with_bp_workers(xcms::findChromPeaks, raw_data, param = xcms::CentWaveParam())

  # A file with zero detected peaks has no computable m/z profile range,
  # which crashes adjustRtime()'s obiwarp alignment outright ("'from' must
  # be a finite number") rather than just aligning poorly -- and zero peaks
  # is already the strongest possible "this injection failed" signal, no
  # threshold comparison needed to know that. Pull these out before
  # alignment/correspondence so one dead file doesn't take down the whole
  # group's QC check, and auto-flag them directly.
  n_files <- nrow(qc_sheet)
  peak_counts <- tabulate(xcms::chromPeaks(xdata)[, "sample"], nbins = n_files)
  no_peaks <- peak_counts == 0
  good_idx <- which(!no_peaks)

  if (any(no_peaks)) {
    message(sprintf(
      "%d file(s) had zero chromatographic peaks under default params -- excluding from alignment, auto-flagging as likely failed injections:\n  %s",
      sum(no_peaks), paste(qc_sheet$sample_label[no_peaks], collapse = "\n  ")
    ))
  }

  aligned_feature_count <- rep(NA_real_, n_files)
  if (length(good_idx) >= 2) {
    xdata_good <- if (any(no_peaks)) xcms::filterFile(xdata, file = good_idx) else xdata
    qc_sheet_good <- qc_sheet[good_idx, , drop = FALSE]

    message("Running adjustRtime()...")
    xdata_good <- with_bp_workers(xcms::adjustRtime, xdata_good, param = xcms::ObiwarpParam())

    xdata_good <- with_bp_workers(
      xcms::groupChromPeaks, xdata_good, param = xcms::PeakDensityParam(sampleGroups = qc_sheet_good$sample_type)
    )

    feature_vals <- xcms::featureValues(xdata_good, value = "into")
    aligned_feature_count[good_idx] <- colSums(!is.na(feature_vals))
  } else {
    message("Fewer than 2 files with any peaks -- skipping alignment, feature-count check unavailable.")
  }

  # TIC prefers an absolute, instrument-specific floor over a
  # batch-distribution threshold: median/MAD both break down once a
  # majority of a batch's QCs are actually bad. No equivalent absolute
  # reference exists for feature count (method/sample-dependent, not just
  # instrument-dependent), so that one stays relative.
  instrument <- group_instrument(qc_sheet)
  instrument_config <- get_instrument_params(instrument, qc_sheet$column[1], qc_sheet$polarity[1])
  has_absolute_tic_floor <- !is.null(instrument_config) && !is.null(instrument_config$int_threshold)

  # Thresholds are computed only from files that actually got peak-picked
  # -- the auto-flagged zero-peak files would otherwise drag down the
  # median/floor for everyone else being compared against it.
  if (has_absolute_tic_floor) {
    tic_threshold <- instrument_config$int_threshold
    message(sprintf("Using instrument-specific absolute TIC threshold: %s", tic_threshold))
  } else {
    tic_threshold <- low_outlier_threshold(tic[good_idx], min_fraction_of_median)
  }
  feature_threshold <- if (length(good_idx) >= 2) {
    low_outlier_threshold(aligned_feature_count[good_idx], min_fraction_of_median)
  } else {
    NA_real_
  }

  tic_flag_global <- tic < tic_threshold
  feature_flag_global <- aligned_feature_count < feature_threshold
  tic_flag_global[no_peaks] <- TRUE
  feature_flag_global[no_peaks] <- TRUE

  # Per-batch check: a batch could drift as a whole relative to the rest of
  # the group, in which case a file can look "normal" globally while still
  # being an outlier within its own batch (or vice versa) — flag on either.
  # Skipped for batches with too few QC files for a reliable local
  # threshold; the absolute TIC floor (when available) applies per-batch
  # too, since it's the same physical floor regardless of scope. Restricted
  # to peak-picked files for the same reason thresholds above are.
  batches <- qc_sheet$batch
  tic_flag_batch <- rep(FALSE, length(tic))
  feature_flag_batch <- rep(FALSE, length(aligned_feature_count))

  for (b in unique(batches)) {
    idx <- which(batches == b & !no_peaks)
    if (length(idx) < min_batch_qc) next

    batch_tic_threshold <- if (has_absolute_tic_floor) {
      instrument_config$int_threshold
    } else {
      low_outlier_threshold(tic[idx], min_fraction_of_median)
    }
    batch_feature_threshold <- low_outlier_threshold(aligned_feature_count[idx], min_fraction_of_median)

    tic_flag_batch[idx] <- tic[idx] < batch_tic_threshold
    feature_flag_batch[idx] <- aligned_feature_count[idx] < batch_feature_threshold
  }

  tic_flag <- tic_flag_global | tic_flag_batch
  feature_flag <- feature_flag_global | feature_flag_batch

  scope_label <- function(flag_global, flag_batch) {
    ifelse(
      flag_global & flag_batch, "global+batch",
      ifelse(flag_global, "global", ifelse(flag_batch, "batch", NA_character_))
    )
  }
  tic_reason <- ifelse(tic_flag, sprintf("low TIC (%s)", scope_label(tic_flag_global, tic_flag_batch)), NA)
  feature_reason <- ifelse(
    feature_flag,
    sprintf("low aligned feature count (%s)", scope_label(feature_flag_global, feature_flag_batch)),
    NA
  )
  reason <- mapply(function(a, b) {
    parts <- c(a, b)[!is.na(c(a, b))]
    if (length(parts) == 0) NA_character_ else paste(parts, collapse = " and ")
  }, tic_reason, feature_reason)
  reason[no_peaks] <- "no chromatographic peaks detected (likely failed injection or empty vial)"

  result <- data.frame(
    filepath = qc_sheet$filepath,
    sample_label = qc_sheet$sample_label,
    batch = qc_sheet$batch,
    injection_order = qc_sheet$injection_order,
    tic = tic,
    aligned_feature_count = aligned_feature_count,
    flagged_global = tic_flag_global | feature_flag_global,
    flagged_batch = tic_flag_batch | feature_flag_batch,
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
#' - X-axis is `injection_order` (see `scan_mzml_files()`) -- a global,
#'   chronological rank across the whole dataset, preferring true mzML
#'   acquisition timestamps over the filename-encoded number when every
#'   file's timestamp is readable. Unlike batch+sample_label, this is
#'   guaranteed unique on its own: a batch commonly splits into 2+ runlist
#'   halves, so the same QC name (e.g. "ltQC01") can legitimately repeat
#'   within one batch without a plate token to disambiguate it, which
#'   otherwise both misrepresents run order and can crash plot_ly() outright
#'   (duplicate factor levels).
#' - Colored by batch (fixed HCL hue order). One real trace with precomputed
#'   per-bar colors -- one trace per batch is a known plotly footgun that
#'   mis-aligns/stacks bars. Invisible per-batch "legend only" traces supply
#'   the clickable legend instead.
#' - Flagged = faded (reduced alpha), passed = opaque.
#' - Hover shows sample/batch/injection order/value/status/reason.
#' - Threshold drawn as a labeled line (global only; per-batch scope shown
#'   in hover `reason` instead).
#'
#' @param values Numeric vector to plot (e.g. tic or aligned_feature_count).
#' @param flagged Logical vector, same length/order as `values`.
#' @param batches Character vector of batch labels, same length/order.
#' @param reasons Character vector (or NA), same length/order — why a file
#'   was flagged, if it was.
#' @param threshold Numeric threshold to draw as a horizontal line.
#' @param labels Character vector of bar labels (e.g. sample_label).
#' @param injection_order Numeric vector, same length/order as `values` —
#'   the x-axis position (see above).
#' @param title Plot title.
#' @param ylab Y-axis label.
#' @return A plotly htmlwidget.
plot_qc_metric <- function(values, flagged, batches, reasons, threshold, labels, injection_order, title, ylab) {
  unique_batches <- sort(unique(batches))
  batch_colors <- stats::setNames(
    grDevices::hcl.colors(length(unique_batches), palette = "Dynamic"), unique_batches
  )

  bar_colors <- mapply(
    function(b, f) to_rgba(batch_colors[[b]], if (f) 0.3 else 1), batches, flagged
  )
  status <- ifelse(flagged, sprintf("Flagged (%s)", reasons), "Passed")
  hover <- sprintf(
    "%s<br>Batch: %s<br>Injection order: %d<br>%s: %.4g<br>%s",
    labels, batches, injection_order, ylab, values, status
  )

  # plot_ly()'s bar trace silently drops any point whose y is NA ("Ignoring
  # N observations") rather than drawing a zero-height bar -- e.g. the
  # aligned-feature-count chart for a file check_qc_quality() excluded from
  # correspondence entirely (zero peaks). Hover text already renders NA
  # values as literal "NA" (see sprintf("%.4g", NA) above), so only the
  # plotted height needs substituting -- the underlying data stays NA.
  plot_values <- ifelse(is.na(values), 0, values)

  p <- plotly::plot_ly(
    x = injection_order, y = plot_values, type = "bar",
    marker = list(color = bar_colors, line = list(width = 0)),
    text = hover, hoverinfo = "text", showlegend = FALSE
  )

  # Invisible dummy traces just to populate a clickable batch-color legend.
  for (b in unique_batches) {
    p <- plotly::add_trace(
      p, x = injection_order[1], y = 0, type = "bar", name = b,
      marker = list(color = to_rgba(batch_colors[[b]], 1)),
      visible = "legendonly", hoverinfo = "skip", showlegend = TRUE
    )
  }

  plotly::layout(
    p,
    title = title,
    xaxis = list(title = "Injection order", type = "linear"),
    yaxis = list(title = ylab),
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
#' exact values and flag reason.
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
      result$tic, result$flagged, result$batch, result$reason, attr(result, "tic_threshold"),
      result$sample_label, result$injection_order,
      sprintf("%s — Total ion current (TIC)", group_name), "TIC"
    )
    plots[[length(plots) + 1]] <- plot_qc_metric(
      result$aligned_feature_count, result$flagged, result$batch, result$reason,
      attr(result, "feature_threshold"),
      result$sample_label, result$injection_order,
      sprintf("%s — Aligned feature count", group_name), "Feature count"
    )
  }

  htmltools::save_html(htmltools::tagList(plots), out_path)
  message(sprintf("Saved QC quality report to: %s", out_path))
}
