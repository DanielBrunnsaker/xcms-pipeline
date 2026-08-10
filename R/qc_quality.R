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

#' Compute a low-outlier threshold separately for each sample_type present,
#' from only the "usable" subset of rows (e.g. files that actually got
#' peak-picked) -- sQC and ltQC are different matrices that can have
#' systematically different baselines for either metric (e.g. ltQC often
#' more concentrated than sQC, which would otherwise inflate a raw-peak-
#' count/TIC-style threshold), so pooling them into one threshold would
#' misjudge whichever type sits off that shared median. See
#' `check_qc_quality()`'s docstring for the full reasoning.
#'
#' @param x Numeric vector (e.g. raw_peak_count or aligned_feature_count).
#' @param sample_types Character vector, same length as `x`.
#' @param usable_idx Integer indices into `x`/`sample_types` to compute
#'   from (e.g. files that got peak-picked).
#' @param min_fraction_of_median Passed to `low_outlier_threshold()`.
#' @return Named numeric vector, one entry per unique value in
#'   `sample_types`, NA for any type with zero usable rows.
compute_type_thresholds <- function(x, sample_types, usable_idx, min_fraction_of_median = 0.5) {
  types_present <- unique(sample_types)
  thresholds <- stats::setNames(rep(NA_real_, length(types_present)), types_present)

  for (type in types_present) {
    type_idx <- usable_idx[sample_types[usable_idx] == type]
    if (length(type_idx) > 0) {
      thresholds[[type]] <- low_outlier_threshold(x[type_idx], min_fraction_of_median)
    }
  }

  thresholds
}

#' Check QC injection quality via raw and aligned peak counts (default
#' centWave params -- just needs to separate "clearly fine" from "clearly
#' broken"). Flags a file if either metric falls well below its peers -- the
#' signature of a missed injection, empty vial, or degraded run.
#'
#' - Raw peak count (per file, straight from `findChromPeaks()`, before any
#'   alignment) is the primary signal. TIC (a straight sum of intensity)
#'   isn't used at all -- it scales ~linearly with concentration, so two
#'   perfectly healthy QC types (or even two healthy injections of the same
#'   type) can differ a lot in TIC alone (e.g. ltQC often more concentrated
#'   than sQC). Peak count is far less sensitive to that: a real compound
#'   above the detection floor gets counted whether the sample is somewhat
#'   more or less concentrated, and only drops once concentration falls far
#'   enough to push peaks below that floor -- a better match for "did this
#'   injection clearly fail" than "is the signal magnitude a bit different."
#' - Aligned feature count (per file, after `adjustRtime()`/
#'   `groupChromPeaks()`) is a secondary, much more lenient check
#'   (`aligned_min_fraction_of_median`, far lower than
#'   `min_fraction_of_median`): a QC whose peaks don't correspond well with
#'   its peers is a poor anchor for downstream normalization/batch
#'   correction even when the raw injection itself was fine, so it's still
#'   worth flagging -- just not conflated with "this injection failed."
#' - Groups by `sample_type` for correspondence (`groupChromPeaks()`) AND
#'   for both metrics' outlier thresholds: sQC files are judged against
#'   other sQC, ltQC against other ltQC, never pooled together.
#' - Skips `fillChromPeaks()` -- gap-filling would paper over the weak
#'   signal this check exists to catch.
#'
#' @param qc_sheet Sample sheet rows for all QC (sQC/ltQC) files in a group.
#' @param min_fraction_of_median Threshold for the primary raw-peak-count
#'   check (e.g. 0.5 flags anything under half its type's median).
#' @param aligned_min_fraction_of_median Threshold for the secondary,
#'   much-more-lenient aligned-feature-count check.
#' @param min_batch_qc Minimum QC files (of the same sample_type) for a
#'   batch to get its own local threshold (below this, relies on the
#'   global check alone).
#' @return A data.frame with one row per file: filepath, sample_label,
#'   batch, sample_type, injection_order, raw_peak_count,
#'   aligned_feature_count, flagged_global, flagged_batch, flagged, reason.
#'   `raw_peak_count_threshold`/`feature_threshold` attributes carry the
#'   actual thresholds used, as named vectors keyed by sample_type (so
#'   reporting doesn't need to recompute them).
check_qc_quality <- function(qc_sheet, min_fraction_of_median = 0.5,
                              aligned_min_fraction_of_median = 0.3, min_batch_qc = 2) {
  message(sprintf("Reading %d QC file(s) for quality check...", nrow(qc_sheet)))
  raw_data <- read_raw_data(qc_sheet$filepath, get_spectrum_modes(qc_sheet))

  # Looked up once, up front -- used for the centWave params below.
  instrument <- group_instrument(qc_sheet)
  instrument_config <- get_instrument_params(instrument, qc_sheet$column[1], qc_sheet$polarity[1])
  has_instrument_config <- !is.null(instrument_config)

  # xcms::CentWaveParam()'s bare defaults (ppm=25 among others) are generic
  # library defaults, not tuned to any real instrument -- reuse the same
  # instrument-characterized starting_values used to seed the IPO2 search
  # instead, when available.
  centwave_param <- if (has_instrument_config) {
    message("Using instrument-specific centWave params for quality check.")
    build_centwave_param(instrument_config$starting_values)
  } else {
    xcms::CentWaveParam()
  }

  message("Running findChromPeaks() for quality check...")
  xdata <- with_bp_workers(xcms::findChromPeaks, raw_data, param = centwave_param)

  n_files <- nrow(qc_sheet)
  raw_peak_count <- tabulate(xcms::chromPeaks(xdata)[, "sample"], nbins = n_files)
  no_peaks <- raw_peak_count == 0
  good_idx <- which(!no_peaks)

  if (any(no_peaks)) {
    message(sprintf(
      "%d file(s) had zero chromatographic peaks under default params -- excluding from alignment, auto-flagging as likely failed injections:\n  %s",
      sum(no_peaks), paste(qc_sheet$sample_label[no_peaks], collapse = "\n  ")
    ))
  }

  # A file with zero detected peaks has no computable m/z profile range,
  # which crashes adjustRtime()'s obiwarp alignment outright ("'from' must
  # be a finite number") rather than just aligning poorly. Excluded before
  # alignment/correspondence so one dead file doesn't take down the whole
  # group's QC check -- it's already caught by the raw-peak-count check
  # above regardless.
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
    message("Fewer than 2 files with any peaks -- skipping alignment, aligned-feature-count check unavailable.")
  }

  # Thresholds computed separately per sample_type -- sQC and ltQC can have
  # systematically different baselines (different matrices), so pooling
  # them would judge one type against a threshold that's really typical of
  # the other. Restricted to files that actually got peak-picked -- the
  # auto-flagged zero-peak files would otherwise drag down the median for
  # everyone else being compared against it.
  sample_types <- qc_sheet$sample_type
  types_present <- unique(sample_types)

  raw_peak_count_threshold <- compute_type_thresholds(raw_peak_count, sample_types, good_idx, min_fraction_of_median)
  feature_threshold <- compute_type_thresholds(
    aligned_feature_count, sample_types, good_idx, aligned_min_fraction_of_median
  )

  raw_flag_global <- raw_peak_count < raw_peak_count_threshold[sample_types]
  feature_flag_global <- aligned_feature_count < feature_threshold[sample_types]
  # NA only when a type had zero peak-picked files -- every row of that
  # type is, by construction, a no_peaks row, so the next two lines already
  # force it TRUE regardless; this is just a safety net against NA leaking
  # into `flagged` below.
  raw_flag_global[is.na(raw_flag_global)] <- FALSE
  feature_flag_global[is.na(feature_flag_global)] <- FALSE
  raw_flag_global[no_peaks] <- TRUE
  feature_flag_global[no_peaks] <- TRUE

  # Per-batch check: a batch could drift as a whole relative to the rest of
  # the group, in which case a file can look "normal" globally while still
  # being an outlier within its own batch (or vice versa) — flag on either.
  # Local threshold computed within the same (batch, sample_type) pair, for
  # the same reason the global one is split by type above. Skipped for
  # batch/type pairs with too few QC files for a reliable local threshold.
  batches <- qc_sheet$batch
  raw_flag_batch <- rep(FALSE, length(raw_peak_count))
  feature_flag_batch <- rep(FALSE, length(aligned_feature_count))

  for (b in unique(batches)) {
    for (type in types_present) {
      idx <- which(batches == b & sample_types == type & !no_peaks)
      if (length(idx) < min_batch_qc) next

      batch_raw_threshold <- low_outlier_threshold(raw_peak_count[idx], min_fraction_of_median)
      batch_feature_threshold <- low_outlier_threshold(aligned_feature_count[idx], aligned_min_fraction_of_median)

      raw_flag_batch[idx] <- raw_peak_count[idx] < batch_raw_threshold
      feature_flag_batch[idx] <- aligned_feature_count[idx] < batch_feature_threshold
    }
  }

  raw_flag <- raw_flag_global | raw_flag_batch
  feature_flag <- feature_flag_global | feature_flag_batch

  scope_label <- function(flag_global, flag_batch) {
    ifelse(
      flag_global & flag_batch, "global+batch",
      ifelse(flag_global, "global", ifelse(flag_batch, "batch", NA_character_))
    )
  }
  raw_reason <- ifelse(raw_flag, sprintf("low raw peak count (%s)", scope_label(raw_flag_global, raw_flag_batch)), NA)
  feature_reason <- ifelse(
    feature_flag,
    sprintf("low aligned feature count (%s)", scope_label(feature_flag_global, feature_flag_batch)),
    NA
  )
  reason <- mapply(function(a, b) {
    parts <- c(a, b)[!is.na(c(a, b))]
    if (length(parts) == 0) NA_character_ else paste(parts, collapse = " and ")
  }, raw_reason, feature_reason)
  reason[no_peaks] <- "no chromatographic peaks detected (likely failed injection or empty vial)"

  result <- data.frame(
    filepath = qc_sheet$filepath,
    sample_label = qc_sheet$sample_label,
    batch = qc_sheet$batch,
    sample_type = qc_sheet$sample_type,
    injection_order = qc_sheet$injection_order,
    raw_peak_count = raw_peak_count,
    aligned_feature_count = aligned_feature_count,
    flagged_global = raw_flag_global | feature_flag_global,
    flagged_batch = raw_flag_batch | feature_flag_batch,
    flagged = raw_flag | feature_flag,
    reason = reason,
    stringsAsFactors = FALSE
  )
  attr(result, "raw_peak_count_threshold") <- raw_peak_count_threshold
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
#' - Hover shows sample/batch/type/injection order/value/status/reason.
#' - Threshold drawn as one labeled line per sample_type (global scope
#'   only; per-batch scope shown in hover `reason` instead) -- thresholds
#'   are computed separately per type (see `check_qc_quality()`), so a
#'   single shared line would misrepresent whichever type it wasn't
#'   computed from.
#'
#' @param values Numeric vector to plot (e.g. raw_peak_count or
#'   aligned_feature_count).
#' @param flagged Logical vector, same length/order as `values`.
#' @param batches Character vector of batch labels, same length/order.
#' @param reasons Character vector (or NA), same length/order — why a file
#'   was flagged, if it was.
#' @param sample_types Character vector (e.g. "sQC"/"ltQC"), same
#'   length/order as `values`.
#' @param thresholds Named numeric vector of thresholds to draw as
#'   horizontal lines, one per `sample_type` value present (see
#'   `check_qc_quality()`'s `raw_peak_count_threshold`/`feature_threshold`
#'   attributes).
#' @param labels Character vector of bar labels (e.g. sample_label).
#' @param injection_order Numeric vector, same length/order as `values` —
#'   the x-axis position (see above).
#' @param title Plot title.
#' @param ylab Y-axis label.
#' @return A plotly htmlwidget.
plot_qc_metric <- function(values, flagged, batches, reasons, sample_types, thresholds,
                            labels, injection_order, title, ylab) {
  unique_batches <- sort(unique(batches))
  batch_colors <- stats::setNames(
    grDevices::hcl.colors(length(unique_batches), palette = "Dynamic"), unique_batches
  )

  bar_colors <- mapply(
    function(b, f) to_rgba(batch_colors[[b]], if (f) 0.3 else 1), batches, flagged
  )
  status <- ifelse(flagged, sprintf("Flagged (%s)", reasons), "Passed")
  hover <- sprintf(
    "%s<br>Batch: %s<br>Type: %s<br>Injection order: %d<br>%s: %.4g<br>%s",
    labels, batches, sample_types, injection_order, ylab, values, status
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

  # A type with zero peak-picked files has no threshold to draw (NA) --
  # skip it rather than plot a line at NA.
  threshold_types <- names(thresholds)[!is.na(thresholds)]

  plotly::layout(
    p,
    title = title,
    xaxis = list(title = "Injection order", type = "linear"),
    yaxis = list(title = ylab),
    shapes = lapply(threshold_types, function(type) {
      list(
        type = "line", x0 = 0, x1 = 1, xref = "paper",
        y0 = thresholds[[type]], y1 = thresholds[[type]], line = list(dash = "dash", color = "gray")
      )
    }),
    annotations = lapply(threshold_types, function(type) {
      list(
        x = 1, y = thresholds[[type]], xref = "paper", yref = "y",
        text = sprintf("%s threshold = %.3g", type, thresholds[[type]]),
        showarrow = FALSE, xanchor = "right", yanchor = "bottom",
        font = list(size = 10, color = "gray")
      )
    })
  )
}

#' Generate an interactive HTML report of QC quality check results: one
#' raw-peak-count chart and one aligned-feature-count chart per group,
#' colored by batch (faded = flagged), pass/fail threshold(s) marked on
#' both, hover for exact values and flag reason.
#'
#' @param results_by_group Named list of `check_qc_quality()` result
#'   data.frames, one per group (each carrying `raw_peak_count_threshold`/
#'   `feature_threshold` attributes — see `check_qc_quality()`).
#' @param out_path Path to write the HTML report to.
generate_qc_report <- function(results_by_group, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  plots <- list()
  for (group_name in names(results_by_group)) {
    result <- results_by_group[[group_name]]

    plots[[length(plots) + 1]] <- plot_qc_metric(
      result$raw_peak_count, result$flagged, result$batch, result$reason, result$sample_type,
      attr(result, "raw_peak_count_threshold"),
      result$sample_label, result$injection_order,
      sprintf("%s — Raw peak count", group_name), "Raw peak count"
    )
    plots[[length(plots) + 1]] <- plot_qc_metric(
      result$aligned_feature_count, result$flagged, result$batch, result$reason, result$sample_type,
      attr(result, "feature_threshold"),
      result$sample_label, result$injection_order,
      sprintf("%s — Aligned feature count", group_name), "Feature count"
    )
  }

  htmltools::save_html(htmltools::tagList(plots), out_path)
  message(sprintf("Saved QC quality report to: %s", out_path))
}
