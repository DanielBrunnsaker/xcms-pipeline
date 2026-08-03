# IPO2-based optimization of xcms centWave parameters, and peak picking using
# the optimized parameters, run independently per column x polarity group.
#
# Uses the IPO2 package (github.com/wmoldham/IPO2) rather than the original
# Bioconductor IPO package — IPO2 is a from-scratch reimplementation targeting
# the modern findChromPeaks()/CentWaveParam workflow directly, whereas IPO
# routes through xcms's old, unmaintained xcmsSet() API, which broke against
# current xcms/BiocParallel (bpstopOnError could not be found).
# Install with: renv::install("wmoldham/IPO2")
#
# IPO2's own internal parameter-scoring loop calls BiocParallel::bplapply()
# without specifying a backend, so it uses whatever the default registered
# one is — typically MulticoreParam (fork-based) on macOS/Linux, which hits
# a known BiocParallel bug ("wrong args for environment subassignment",
# https://github.com/Bioconductor/BiocParallel/issues/206) that fails every
# worker in a batch identically. Registering SerialParam as the default here
# avoids that broken path; it doesn't affect our own findChromPeaks() calls
# elsewhere, since those already explicitly pass their own SerialParam.
BiocParallel::register(BiocParallel::SerialParam())

#' Pick a small representative subset of files to run parameter optimization
#' on.
#'
#' The design-of-experiments search evaluates many parameter combinations,
#' so it's run on a handful of files rather than the full batch. Preference
#' order:
#'   1. sQC — pooled from the study's own samples, so it's the same matrix
#'      injected repeatedly, giving the optimizer a stable, reproducible,
#'      well-populated signal to tune against.
#'   2. ltQC — a long-term reference matrix from outside the study; still
#'      QC-like reproducibility, but not representative of this study's
#'      chemistry, so only used when sQC isn't available in enough numbers.
#'   3. Regular study samples — heterogeneous by design, which works against
#'      the optimizer's reproducibility-based scoring, so this is a last
#'      resort.
#'
#' Picks are spread evenly across `injection_order` (time) rather than just
#' taking the first `n` rows — a group can span many batches (e.g. in
#' "global" IPO scope), and the sheet is sorted by batch, so naively taking
#' the first `n` QC rows would only ever sample the earliest batch instead
#' of being representative of the whole group.
#'
#' QC rows flagged as faulty by `scripts/check_qc_quality.R` (a `qc_flagged`
#' column with value `TRUE`) are excluded from consideration entirely — a
#' missed injection or empty vial should never end up as the "representative"
#' QC. If that check was never run, the sheet won't have the column, and
#' every QC row is treated as usable.
#'
#' @param group_sheet Sample sheet rows for a single column x polarity group.
#' @param n Number of files to select.
#' @param min_qc Minimum number of QC rows required before they're considered
#'   "enough" to use (below this, fall through to the next tier).
select_ipo_subset <- function(group_sheet, n = 4, min_qc = 2) {
  pick_spread <- function(rows, n) {
    rows <- rows[order(rows$injection_order), ]
    idx <- unique(round(seq(1, nrow(rows), length.out = n)))
    rows$filepath[idx]
  }

  if (!"qc_flagged" %in% names(group_sheet)) {
    group_sheet$qc_flagged <- FALSE
  }
  group_sheet$qc_flagged[is.na(group_sheet$qc_flagged)] <- FALSE
  group_sheet <- group_sheet[!group_sheet$qc_flagged, , drop = FALSE]

  sqc_rows <- group_sheet[group_sheet$sample_type == "sQC", , drop = FALSE]
  if (nrow(sqc_rows) >= min_qc) {
    return(pick_spread(sqc_rows, min(n, nrow(sqc_rows))))
  }

  ltqc_rows <- group_sheet[group_sheet$sample_type == "ltQC", , drop = FALSE]
  if (nrow(ltqc_rows) >= min_qc) {
    return(pick_spread(ltqc_rows, min(n, nrow(ltqc_rows))))
  }

  regular_rows <- group_sheet[group_sheet$sample_type == "Sample", , drop = FALSE]
  n <- min(n, nrow(regular_rows))
  regular_rows$filepath[sample(nrow(regular_rows), n)]
}

#' Run IPO2 optimization for one column x polarity group, caching the result
#' to disk so re-running peak picking doesn't re-optimize from scratch.
#'
#' Uses instrument/method-specific starting values and search-space bounds
#' from `R/instrument_params.R` when the sheet's `instrument` column (and a
#' matching config entry) is available, falling back to IPO2's own defaults
#' otherwise.
#'
#' @param group_sheet Sample sheet rows for this group.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @return An xcms::CentWaveParam with the optimized settings.
run_ipo_optimization <- function(group_sheet, out_dir) {
  params_path <- file.path(out_dir, "ipo_params.rds")

  if (file.exists(params_path)) {
    message(sprintf("Using cached IPO2 params: %s", params_path))
    return(readRDS(params_path))
  }

  subset_files <- select_ipo_subset(group_sheet)
  message(sprintf(
    "Running IPO2 optimization on %d file(s):\n  %s",
    length(subset_files), paste(subset_files, collapse = "\n  ")
  ))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  instrument <- group_instrument(group_sheet)
  instrument_config <- get_instrument_params(
    instrument, group_sheet$column[1], group_sheet$polarity[1]
  )

  subset_modes <- get_spectrum_modes(group_sheet)[match(subset_files, group_sheet$filepath)]
  raw_data <- read_raw_data(subset_files, subset_modes)
  ipo_result <- if (!is.null(instrument_config)) {
    message(sprintf(
      "Using instrument-specific search space for %s / %s / %s",
      instrument, group_sheet$column[1], group_sheet$polarity[1]
    ))
    IPO2::optimize_centwave(
      raw_data = raw_data,
      parameter_list = build_ipo2_parameter_list(instrument_config),
      out_dir = out_dir
    )
  } else {
    message("No instrument-specific config found; using IPO2 defaults.")
    IPO2::optimize_centwave(raw_data = raw_data, out_dir = out_dir)
  }
  centwave_param <- ipo_result$best_cwp

  saveRDS(centwave_param, params_path)

  # Keep the full search history too, not just the winning params — useful
  # for reviewing how the optimizer got there. Saved as both .rds (full
  # fidelity, including the CentWaveParam object tried at each step) and
  # .csv (the `cwp` column dropped, since a CentWaveParam object isn't
  # CSV-representable, for quick eyeballing without loading R).
  history <- ipo_result$history
  saveRDS(history, file.path(out_dir, "ipo_history.rds"))
  history_csv <- as.data.frame(history)
  history_csv$cwp <- NULL
  write.csv(history_csv, file.path(out_dir, "ipo_history.csv"), row.names = FALSE)

  message(sprintf("Saved optimized params to: %s", params_path))

  centwave_param
}

#' Run findChromPeaks() for a set of files with given params. A building
#' block reused both for a whole group (global IPO scope) and for a single
#' batch (batch IPO scope, before the per-batch results get combined).
#'
#' @param filepaths Character vector of mzML file paths.
#' @param centwave_param xcms::CentWaveParam to use for findChromPeaks().
#' @param spectrum_modes Character vector ("profile"/"centroid"/NA), one per
#'   file in `filepaths`, used to centroid profile-mode files before peak
#'   picking. Defaults to NA (no centroiding) for every file if omitted.
#' @return An XCMSnExp with peaks picked (no alignment/correspondence yet).
pick_peaks <- function(filepaths, centwave_param, spectrum_modes = NULL) {
  if (is.null(spectrum_modes)) {
    spectrum_modes <- rep(NA_character_, length(filepaths))
  }

  message(sprintf("Reading %d file(s) for peak picking...", length(filepaths)))
  raw_data <- read_raw_data(filepaths, spectrum_modes)

  message("Running findChromPeaks()...")
  xcms::findChromPeaks(
    raw_data, param = centwave_param,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
  )
}

#' Align retention times, group peaks into cross-sample features, and fill
#' gaps — turning individually-picked peaks into one aligned feature table.
#' Follows the same steps as github.com/MetaboComp/xcms_pipeline
#' (pure_xcms_pipeline.R): adjustRtime(ObiwarpParam) ->
#' groupChromPeaks(PeakDensityParam) -> fillChromPeaks(ChromPeakAreaParam).
#' Uses xcms's built-in defaults for all three params for now (untuned).
#'
#' `sampleGroups` mirrors the reference pipeline's actual usage: despite the
#' name, their "sample_group" is QC/sample TYPE classification (sQC/ltQC/
#' Blank/Sample), not biological condition — so this uses our own
#' already-populated `sample_type` column for the same purpose, rather than
#' the sheet's manually-curated (and possibly unfilled) `sample_group`.
#'
#' @param xdata XCMSnExp with peaks already picked (via `pick_peaks()`, or
#'   several combined with `xcms::c()`).
#' @param sample_types Character vector of `sample_type` values, one per
#'   file, in the same order as files in `xdata`.
#' @return The aligned, corresponded, gap-filled XCMSnExp.
align_and_correspond <- function(xdata, sample_types) {
  message("Running adjustRtime()...")
  xdata <- xcms::adjustRtime(xdata, param = xcms::ObiwarpParam())

  message("Running groupChromPeaks()...")
  xdata <- xcms::groupChromPeaks(
    xdata, param = xcms::PeakDensityParam(sampleGroups = sample_types)
  )

  message("Running fillChromPeaks()...")
  xcms::fillChromPeaks(xdata, param = xcms::ChromPeakAreaParam())
}

#' Combine an aligned XCMSnExp's feature definitions and per-sample values
#' into one flat table, with real sample names instead of raw file paths.
#'
#' @param xdata Aligned XCMSnExp (after `align_and_correspond()`).
#' @param sample_names Character vector of sample names, one per file, in
#'   the same order as files in `xdata`.
#' @return A data.frame: one row per aligned feature, one column per sample.
build_feature_table <- function(xdata, sample_names) {
  feature_defs <- as.data.frame(xcms::featureDefinitions(xdata))
  feature_defs$peakidx <- NULL

  feature_values <- as.data.frame(xcms::featureValues(xdata, value = "into"))
  colnames(feature_values) <- sample_names

  cbind(feature_defs, feature_values)
}

#' Save all peak-picking + alignment outputs for a group: the full aligned
#' XCMSnExp object, the flat per-peak table (including gap-filled entries),
#' and the aligned feature table.
#'
#' @param xdata Aligned XCMSnExp (after `align_and_correspond()`).
#' @param feature_table Data.frame from `build_feature_table()`.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
save_peak_picking_outputs <- function(xdata, feature_table, out_dir) {
  peaks_dir <- file.path(out_dir, "peaks")
  dir.create(peaks_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(xdata, file.path(peaks_dir, "xdata.rds"))

  peak_table <- as.data.frame(xcms::chromPeaks(xdata))
  write.csv(peak_table, file.path(peaks_dir, "peak_table.csv"), row.names = FALSE)

  write.csv(feature_table, file.path(out_dir, "feature_table.csv"), row.names = FALSE)

  message(sprintf(
    "Saved %d peaks and %d aligned features to: %s",
    nrow(peak_table), nrow(feature_table), out_dir
  ))
}
