# IPO2-based optimization of xcms centWave parameters, and peak picking using
# the optimized parameters, run independently per column x polarity group.
#
# Uses IPO2 (github.com/wmoldham/IPO2), not the original Bioconductor IPO:
# - IPO2 targets the modern findChromPeaks()/CentWaveParam workflow directly.
# - IPO routes through xcms's old xcmsSet() API, which broke against current
#   xcms/BiocParallel (bpstopOnError could not be found).
# Install with: renv::install("wmoldham/IPO2")
#
# IPO2's own scoring loop calls BiocParallel::bplapply() without a backend,
# so it uses whatever's registered as default:
# - MulticoreParam (fork-based) hit a known BiocParallel bug here ("wrong
#   args for environment subassignment", github.com/Bioconductor/
#   BiocParallel/issues/206), failing every worker in a batch identically.
# - fork() doesn't exist on Windows anyway, so MulticoreParam would degrade
#   to single-core there regardless.
# - SnowParam (separate processes over local sockets) was tried next, but
#   hangs indefinitely for IPO2's loop specifically under Docker/WSL2 on
#   Windows (workers spawn fine, CPU stays busy, but nothing ever reports
#   back — a socket-communication hang, not an error, so nothing here can
#   catch or recover from it). Back to SerialParam: proven stable, just
#   slower.
# Our own findChromPeaks()/adjustRtime()/groupChromPeaks()/fillChromPeaks()
# calls elsewhere are unaffected — they ran fine under SnowParam already,
# and explicitly pass their own BPPARAM (bp_workers() in R/parallel.R),
# overriding this default regardless of what's registered here.
BiocParallel::register(BiocParallel::SerialParam(progressbar = TRUE))

#' Pick a small representative subset of files to run parameter optimization
#' on.
#'
#' Runs on a handful of files, not the full batch (the DoE search evaluates
#' many parameter combinations). Preference order:
#' 1. sQC — pooled from the study's own samples, same matrix injected
#'    repeatedly: a stable, reproducible signal to tune against.
#' 2. ltQC — external long-term reference matrix; still QC-like
#'    reproducibility, but not this study's chemistry. Used only if sQC is
#'    insufficient.
#' 3. Regular study samples — heterogeneous by design, works against the
#'    optimizer's reproducibility-based scoring. Last resort.
#'
#' Batch handling:
#' - More batches than `n` (e.g. "global" scope with many batches): the
#'   batches themselves are subsampled evenly across the timespan, one
#'   file from each selected batch.
#' - Fewer (or equal) batches than `n`: `n` is distributed as evenly as
#'   possible across all batches (at least 1 each, more if `n` allows),
#'   each batch's own picks spread across its own injection timing — so a
#'   low batch count never silently hands the optimizer fewer files than
#'   requested (a handful of batches with only 1 file apiece gives IPO2
#'   very little to work with for its reproducibility-based scoring).
#' - Single batch (e.g. "batch" scope): picks spread evenly across
#'   `injection_order` within that batch, as before.
#'
#' QC rows flagged by `scripts/check_qc_quality.R` (`qc_flagged == TRUE`)
#' are excluded — a missed injection/empty vial should never be the
#' "representative" QC. If that check was never run, every QC row counts.
#'
#' @param group_sheet Sample sheet rows for a single column x polarity group.
#' @param n Number of files to select (or, when multiple batches are
#'   present, the max number of batches to draw one file from each).
#' @param min_qc Minimum number of QC rows required before they're considered
#'   "enough" to use (below this, fall through to the next tier).
select_ipo_subset <- function(group_sheet, n = 4, min_qc = 2) {
  pick_spread <- function(rows, k) {
    k <- min(k, nrow(rows))
    idx <- unique(round(seq(1, nrow(rows), length.out = k)))
    rows$filepath[idx]
  }

  pick_representative <- function(rows, n) {
    rows <- rows[order(rows$injection_order), ]
    batches <- unique(rows$batch)

    if (length(batches) <= 1) {
      return(pick_spread(rows, n))
    }

    if (length(batches) > n) {
      # More batches than requested files: subsample the batches
      # themselves, evenly across the timespan, one file from each.
      batch_idx <- unique(round(seq(1, length(batches), length.out = n)))
      batches <- batches[batch_idx]
      return(vapply(batches, function(b) {
        batch_rows <- rows[rows$batch == b, ]
        mid <- round((nrow(batch_rows) + 1) / 2)
        batch_rows$filepath[mid]
      }, character(1)))
    }

    # Fewer (or equal) batches than requested files: distribute n across
    # all batches as evenly as possible, each batch's share spread across
    # its own timeline.
    per_batch_n <- rep(n %/% length(batches), length(batches))
    remainder <- n %% length(batches)
    if (remainder > 0) {
      per_batch_n[seq_len(remainder)] <- per_batch_n[seq_len(remainder)] + 1
    }
    names(per_batch_n) <- batches

    unlist(lapply(batches, function(b) {
      pick_spread(rows[rows$batch == b, ], per_batch_n[[b]])
    }), use.names = FALSE)
  }

  if (!"qc_flagged" %in% names(group_sheet)) {
    group_sheet$qc_flagged <- FALSE
  }
  group_sheet$qc_flagged[is.na(group_sheet$qc_flagged)] <- FALSE
  group_sheet <- group_sheet[!group_sheet$qc_flagged, , drop = FALSE]

  sqc_rows <- group_sheet[group_sheet$sample_type == "sQC", , drop = FALSE]
  if (nrow(sqc_rows) >= min_qc) {
    return(pick_representative(sqc_rows, n))
  }

  ltqc_rows <- group_sheet[group_sheet$sample_type == "ltQC", , drop = FALSE]
  if (nrow(ltqc_rows) >= min_qc) {
    return(pick_representative(ltqc_rows, n))
  }

  regular_rows <- group_sheet[group_sheet$sample_type == "Sample", , drop = FALSE]
  n_regular <- min(n, nrow(regular_rows))
  subset_files <- regular_rows$filepath[sample(nrow(regular_rows), n_regular)]

  if (length(subset_files) == 0) {
    batches <- unique(group_sheet$batch)
    batch_desc <- if (length(batches) == 1) sprintf("batch: %s, ", batches) else ""
    stop(sprintf(
      paste(
        "No usable files for IPO optimization (%scolumn: %s, polarity: %s)",
        "— every QC is flagged/absent and there are no regular samples either.",
        "Review qc_flagged in the sample sheet, or exclude this batch."
      ),
      batch_desc, group_sheet$column[1], group_sheet$polarity[1]
    ), call. = FALSE)
  }

  subset_files
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
#' @param ipo_subset_size Passed to `select_ipo_subset()` as `n` (files, or
#'   batches to draw one file from each). Larger values improve batch
#'   coverage in "global" scope but multiply IPO2's per-trial cost (see the
#'   top of this file for IPO2's parallel-backend caveats).
#' @return An xcms::CentWaveParam with the optimized settings.
run_ipo_optimization <- function(group_sheet, out_dir, ipo_subset_size = 4) {
  params_path <- file.path(out_dir, "ipo_params.rds")

  if (file.exists(params_path)) {
    message(sprintf("Using cached IPO2 params: %s", params_path))
    return(readRDS(params_path))
  }

  subset_files <- select_ipo_subset(group_sheet, n = ipo_subset_size)
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
  with_bp_workers(xcms::findChromPeaks, raw_data, param = centwave_param)
}

#' Turns individually-picked peaks into one aligned feature table: aligns
#' retention times, groups peaks into cross-sample features, fills gaps.
#' Same steps as github.com/MetaboComp/xcms_pipeline (pure_xcms_pipeline.R):
#' adjustRtime(ObiwarpParam) -> groupChromPeaks(PeakDensityParam) ->
#' fillChromPeaks(ChromPeakAreaParam). Uses xcms's built-in defaults for all
#' three (untuned, for now).
#'
#' `sampleGroups` uses our `sample_type` column (sQC/ltQC/Blank/Sample) —
#' mirrors the reference pipeline's actual usage: despite the name, their
#' "sample_group" is QC/sample TYPE classification, not biological
#' condition, and is a better fit than the sheet's manually-curated
#' (possibly unfilled) `sample_group`.
#'
#' @param xdata XCMSnExp with peaks already picked (via `pick_peaks()`, or
#'   several combined with `xcms::c()`).
#' @param sample_types Character vector of `sample_type` values, one per
#'   file, in the same order as files in `xdata`.
#' @return The aligned, corresponded, gap-filled XCMSnExp.
align_and_correspond <- function(xdata, sample_types) {
  message("Running adjustRtime()...")
  xdata <- with_bp_workers(xcms::adjustRtime, xdata, param = xcms::ObiwarpParam())

  message("Running groupChromPeaks()...")
  xdata <- with_bp_workers(
    xcms::groupChromPeaks, xdata, param = xcms::PeakDensityParam(sampleGroups = sample_types)
  )

  message("Running fillChromPeaks()...")
  with_bp_workers(xcms::fillChromPeaks, xdata, param = xcms::ChromPeakAreaParam())
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
