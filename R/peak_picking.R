# IPO2-optimized centWave peak picking, per column x polarity group.
#
# The IPO2/Bioconductor session setup (library(xcms), the optimPP()
# monkeypatch) lives at the BOTTOM of this file, not the top -- it has no
# effect on any function defined above it, and keeping it out of the way
# means this file's pure logic (select_ipo_subset() etc.) can be sourced
# and unit-tested without Bioconductor installed. See tests/testthat/.

#' Whether at least one QC tier (sQC or ltQC, non-flagged) has enough files
#' to be usable on its own for IPO optimization -- same tiering/exclusion
#' rules as `select_ipo_subset()`'s first two tiers, factored out so
#' `run_ipo_optimization()` can decide whether to run IPO2 at all *before*
#' calling `select_ipo_subset()`, rather than let it silently fall through
#' to regular study samples as the optimization target.
#'
#' @param group_sheet Sample sheet rows (same shape `select_ipo_subset()`
#'   takes).
#' @param min_qc Minimum non-flagged files for a QC tier to count as usable.
qc_tier_available <- function(group_sheet, min_qc = 2) {
  if (!"qc_flagged" %in% names(group_sheet)) {
    group_sheet$qc_flagged <- FALSE
  }
  group_sheet$qc_flagged[is.na(group_sheet$qc_flagged)] <- FALSE
  group_sheet <- group_sheet[!group_sheet$qc_flagged, , drop = FALSE]

  sum(group_sheet$sample_type == "sQC") >= min_qc || sum(group_sheet$sample_type == "ltQC") >= min_qc
}

#' Pick a small representative subset of files for parameter optimization.
#'
#' Preference order: sQC (stable, repeated matrix) > ltQC (repeated but not
#' this study's chemistry) > regular samples (heterogeneous, last resort).
#' Flagged QC rows (`qc_flagged == TRUE`) are excluded.
#'
#' Batch handling:
#' - More batches than `n`: subsample batches evenly across time, one file
#'   each.
#' - Fewer/equal batches than `n`: spread `n` evenly across batches (each
#'   batch's share spread across its own timeline), so a low batch count
#'   never hands the optimizer fewer files than requested.
#' - Single batch: spread evenly across `injection_order`.
#'
#' Keeps the selection to a single spectrum mode (centroid preferred --
#' optimXCMS() has no centroiding step, so a mixed-mode subset would need
#' all-or-nothing centroiding that distorts the already-centroid files) --
#' but only when one mode alone can fill the request; mixed as a last
#' resort rather than under-filling.
#'
#' @param group_sheet Sample sheet rows for a single column x polarity group.
#' @param n Number of files to select (or, when multiple batches are
#'   present, the max number of batches to draw one file from each).
#' @param min_qc Minimum QC rows for a tier to count as "enough" (below
#'   this, fall through to the next tier).
select_ipo_subset <- function(group_sheet, n = 4, min_qc = 2) {
  # Prefer whichever spectrum mode has enough candidates to fill the
  # request on its own; mix only if neither alone does. See docstring above.
  prefer_uniform_mode <- function(rows, n) {
    is_profile <- get_spectrum_modes(rows) %in% "profile"
    if (sum(!is_profile) >= n) {
      rows[!is_profile, , drop = FALSE]
    } else if (sum(is_profile) >= n) {
      rows[is_profile, , drop = FALSE]
    } else {
      rows
    }
  }

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

  # pick_representative()/pick_spread() cap silently at however many rows a
  # tier actually has -- clearing min_qc only means the tier is usable at
  # all, not that it has n files. Warn here so a thin tier (e.g. 3 files
  # when 4 were requested) is visible rather than quietly optimizing
  # against fewer files than the caller asked for.
  return_subset <- function(rows, tier_label) {
    picked <- pick_representative(prefer_uniform_mode(rows, n), n)
    if (length(picked) < n) {
      message(sprintf(
        "Only %d %s file(s) available (of %d requested) for column: %s, polarity: %s%s -- optimizing against fewer files than ipo_subset_size.",
        length(picked), tier_label, n, group_sheet$column[1], group_sheet$polarity[1],
        if (length(unique(group_sheet$batch)) == 1) sprintf(", batch: %s", group_sheet$batch[1]) else ""
      ))
    }
    picked
  }

  sqc_rows <- group_sheet[group_sheet$sample_type == "sQC", , drop = FALSE]
  if (nrow(sqc_rows) >= min_qc) {
    return(return_subset(sqc_rows, "sQC"))
  }

  ltqc_rows <- group_sheet[group_sheet$sample_type == "ltQC", , drop = FALSE]
  if (nrow(ltqc_rows) >= min_qc) {
    return(return_subset(ltqc_rows, "ltQC"))
  }

  regular_rows <- group_sheet[group_sheet$sample_type == "Sample", , drop = FALSE]
  if (nrow(regular_rows) == 0) {
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

  return_subset(regular_rows, "regular sample")
}

#' Run IPO2 optimization for one column x polarity group, caching the result
#' to disk so re-running peak picking doesn't re-optimize from scratch.
#'
#' Uses instrument/method-specific starting values and search-space bounds
#' from `R/instrument_params.R` when the sheet's `instrument` column (and a
#' matching config entry) is available, falling back to
#' `default_ipo2_search_space()` otherwise.
#'
#' If neither QC tier (sQC or ltQC alone) has enough non-flagged files (see
#' `qc_tier_available()`), skips the IPO2 search entirely and returns that
#' starting-point CentWaveParam as-is, un-optimized -- rather than let
#' `select_ipo_subset()` fall through to regular study samples as the
#' optimization target.
#'
#' Sets the `xcms_pipeline.centroid_ipo_input` option before calling
#' `optimXCMS()` if any subset file is profile-mode -- read by
#' `patched_optimPP()` (bottom of this file) to centroid in-place, since
#' `optimXCMS()` otherwise reads its files with no centroiding step at all.
#'
#' @param group_sheet Sample sheet rows for this group.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @param ipo_subset_size Passed to `select_ipo_subset()` as `n` (files, or
#'   batches to draw one file from each). Larger values improve batch
#'   coverage in "global" scope but multiply each evaluation's cost.
#' @param fresh If TRUE, ignore (and delete) any cached final result or
#'   mid-search checkpoint from a prior run and re-optimize from scratch,
#'   rather than reusing/resuming from them.
#' @param min_qc Passed to `qc_tier_available()`/`select_ipo_subset()` --
#'   minimum non-flagged files for a QC tier to count as usable.
#' @return An xcms::CentWaveParam with the optimized settings.
run_ipo_optimization <- function(group_sheet, out_dir, ipo_subset_size = 4, fresh = FALSE, min_qc = 2) {
  params_path <- file.path(out_dir, "ipo_params.rds")
  checkpoint_path <- file.path(out_dir, "ipo_checkpoint.rds")

  if (fresh) {
    if (file.exists(params_path) || file.exists(checkpoint_path)) {
      message("fresh = TRUE: ignoring any cached IPO2 result or checkpoint, re-optimizing from scratch.")
    }
    unlink(params_path)
    unlink(checkpoint_path)
  }

  if (file.exists(params_path)) {
    message(sprintf("Using cached IPO2 params: %s", params_path))
    return(readRDS(params_path))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  instrument <- group_instrument(group_sheet)
  instrument_config <- get_instrument_params(
    instrument, group_sheet$column[1], group_sheet$polarity[1]
  )

  search_space <- if (!is.null(instrument_config)) {
    message(sprintf(
      "Using instrument-specific search space for %s / %s / %s",
      instrument, group_sheet$column[1], group_sheet$polarity[1]
    ))
    build_ipo2_search_space(instrument_config)
  } else {
    message("No instrument-specific config found; using generic default search space.")
    default_ipo2_search_space()
  }

  # No usable QC (sQC or ltQC alone) means select_ipo_subset() would fall
  # all the way through to regular study samples as IPO2's optimization
  # target -- heterogeneous biological material, not what optimXCMS()'s
  # search is meant to tune against. Rather than run IPO2 on that, skip the
  # search entirely and use the (already instrument-characterized, when
  # available) starting point as-is.
  if (!qc_tier_available(group_sheet, min_qc = min_qc)) {
    batches <- unique(group_sheet$batch)
    batch_desc <- if (length(batches) == 1) sprintf("batch: %s, ", batches) else ""
    message(sprintf(
      "No usable QC (sQC or ltQC alone, non-flagged) for IPO optimization (%scolumn: %s, polarity: %s) -- skipping the IPO2 search (won't run it against regular study samples) and using the %s centWave starting point as-is.",
      batch_desc, group_sheet$column[1], group_sheet$polarity[1],
      if (!is.null(instrument_config)) "instrument-specific" else "generic default"
    ))
    centwave_param <- search_space$cwParam
    saveRDS(centwave_param, params_path)
    message(sprintf("Saved (un-optimized) params to: %s", params_path))
    return(centwave_param)
  }

  subset_files <- select_ipo_subset(group_sheet, n = ipo_subset_size, min_qc = min_qc)
  message(sprintf(
    "Running IPO2 optimization on %d file(s):\n  %s",
    length(subset_files), paste(subset_files, collapse = "\n  ")
  ))

  subset_modes <- get_spectrum_modes(group_sheet)[match(subset_files, group_sheet$filepath)]
  n_profile <- sum(subset_modes == "profile", na.rm = TRUE)
  needs_centroiding <- n_profile > 0
  if (needs_centroiding) {
    if (n_profile < length(subset_modes)) {
      # select_ipo_subset() only mixes modes as a last resort (neither mode
      # alone had enough candidates) -- centroiding is applied to the whole
      # subset regardless, which would distort the already-centroid files.
      warning(
        "IPO subset mixes profile- and centroid-mode files (neither mode ",
        "alone had enough candidates); centroiding the whole subset would ",
        "distort the already-centroid ones.",
        call. = FALSE
      )
    } else {
      message("IPO subset is profile-mode; centroiding it for the search (see patched_optimPP() at the bottom of this file).")
    }
  }
  old_opts <- options(
    xcms_pipeline.centroid_ipo_input = needs_centroiding,
    xcms_pipeline.ipo_checkpoint_path = checkpoint_path
  )
  on.exit(options(old_opts), add = TRUE)

  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    message(sprintf(
      "Resuming from a checkpoint left by an interrupted prior run (score: %.3f) -- using it as the starting point instead of the default. Delete %s to discard it and start fresh.",
      checkpoint$score, checkpoint_path
    ))
    search_space$cwParam <- checkpoint$cwParam
  }

  optimum <- IPO2::optimXCMS(
    files = subset_files,
    cwParam = search_space$cwParam,
    optimVars = search_space$optimVars,
    upper = search_space$upper,
    lower = search_space$lower,
    algorithm = "NLOPT_LN_NELDERMEAD",
    verbose = TRUE
  )

  centwave_param <- IPO2:::getCentWaveParams(
    par = optimum$solution,
    optimVars = search_space$optimVars,
    customParam = IPO2:::cw2customParam(search_space$cwParam)
  )

  saveRDS(centwave_param, params_path)

  # optimXCMS()/nloptr doesn't expose a per-iteration trace the way the
  # earlier (wmoldham/IPO2) search did -- this is the optimizer's final
  # result summary, not a full search history, but kept under the same
  # filenames for continuity.
  result_summary <- stats::setNames(as.list(optimum$solution), search_space$optimVars)
  result_summary$objective <- -optimum$objective # optimPP() minimizes -score
  result_summary$iterations <- optimum$iterations
  result_summary$status <- optimum$status
  result_summary$message <- optimum$message

  saveRDS(result_summary, file.path(out_dir, "ipo_history.rds"))
  write.csv(as.data.frame(result_summary), file.path(out_dir, "ipo_history.csv"), row.names = FALSE)

  # Search finished successfully and its real result is saved above -- the
  # mid-search checkpoint has served its purpose, don't leave it around to
  # be mistaken for a still-in-progress run.
  unlink(checkpoint_path)

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
#' retention times, groups into cross-sample features, fills gaps. Same
#' steps as the reference pipeline (adjustRtime -> groupChromPeaks ->
#' fillChromPeaks); ChromPeakAreaParam stays at xcms defaults (untuned).
#'
#' - `sampleGroups` uses `sample_type` (sQC/ltQC/Blank/Sample), not the
#'   sheet's `sample_group` -- matches the reference's actual usage and
#'   doesn't depend on that column being filled in.
#' - `minFraction = 0.2` always (xcms default 0.5; matches the reference's
#'   own final `minFrac`), regardless of `retgroup_params`: pooling all
#'   real samples into one "Sample" group means a 50% bar would silently
#'   drop condition-specific features, and filtering is a deliberate later
#'   stage, not implicit here. `run_retgroup_optimization()`'s fixed
#'   `minfrac` (0.8) is a separate, internal-only scoring value -- not used
#'   for this.
#'
#' @param xdata XCMSnExp with peaks already picked (via `pick_peaks()`, or
#'   several combined with `c()`).
#' @param sample_types Character vector of `sample_type` values, one per
#'   file, in the same order as files in `xdata`.
#' @param retgroup_params Optional result of `run_retgroup_optimization()`
#'   (obiwarp params + density bandwidth/bin-size) -- falls back to xcms's
#'   plain defaults for both if omitted.
#' @return The aligned, corresponded, gap-filled XCMSnExp.
align_and_correspond <- function(xdata, sample_types, retgroup_params = NULL) {
  obiwarp_param <- if (!is.null(retgroup_params)) retgroup_params$obiwarp else xcms::ObiwarpParam()

  message("Running adjustRtime()...")
  xdata <- with_bp_workers(xcms::adjustRtime, xdata, param = obiwarp_param)

  density_param <- if (!is.null(retgroup_params)) {
    xcms::PeakDensityParam(
      sampleGroups = sample_types,
      bw = retgroup_params$density_bounds$bw,
      binSize = retgroup_params$density_bounds$binSize,
      minFraction = 0.2
    )
  } else {
    xcms::PeakDensityParam(sampleGroups = sample_types, minFraction = 0.2)
  }

  message("Running groupChromPeaks()...")
  xdata <- with_bp_workers(xcms::groupChromPeaks, xdata, param = density_param)

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

  # xcms names every feature (FT0001, FT0002, ...) as row names rather than
  # a regular column -- save_peak_picking_outputs() writes with
  # row.names = FALSE, so without pulling this out explicitly there'd be no
  # stable ID to trace a row back to peak_table.csv/xdata.rds by.
  feature_defs <- cbind(feature = rownames(feature_defs), feature_defs)

  feature_values <- as.data.frame(xcms::featureValues(xdata, value = "into"))
  colnames(feature_values) <- sample_names

  cbind(feature_defs, feature_values)
}

#' Save all peak-picking + alignment outputs for a group: the full aligned
#' XCMSnExp object, the flat per-peak table (including gap-filled entries),
#' and the aligned feature table.
#'
#' `peak_table.csv` is enriched beyond xcms::chromPeaks()'s raw columns:
#' - `sample_name`: chromPeaks()'s numeric `sample` index, translated.
#' - `is_filled`/`ms_level`: from the separate, parallel-indexed
#'   chromPeakData() -- `is_filled` matters for external match confidence
#'   (a gap-filled value is weaker evidence than a real detection).
#' - `feature`: aligned feature ID (FT0001, ...) this peak belongs to, if
#'   any -- lets you drill from a feature down to its raw per-sample peaks
#'   without reloading xdata.rds.
#'
#' @param xdata Aligned XCMSnExp (after `align_and_correspond()`).
#' @param feature_table Data.frame from `build_feature_table()`.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @param sample_names Character vector of sample names, one per file, in
#'   the same order as files in `xdata` -- same vector passed to
#'   `build_feature_table()`.
save_peak_picking_outputs <- function(xdata, feature_table, out_dir, sample_names) {
  peaks_dir <- file.path(out_dir, "peaks")
  dir.create(peaks_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(xdata, file.path(peaks_dir, "xdata.rds"))

  peak_table <- as.data.frame(xcms::chromPeaks(xdata))
  peak_table$sample_name <- sample_names[peak_table$sample]

  peak_data <- as.data.frame(xcms::chromPeakData(xdata))
  peak_table$is_filled <- peak_data$is_filled
  peak_table$ms_level <- peak_data$ms_level

  feature_defs_raw <- xcms::featureDefinitions(xdata)
  peak_table$feature <- NA_character_
  for (i in seq_len(nrow(feature_defs_raw))) {
    peak_table$feature[feature_defs_raw$peakidx[[i]]] <- rownames(feature_defs_raw)[i]
  }

  write.csv(peak_table, file.path(peaks_dir, "peak_table.csv"), row.names = FALSE)

  write.csv(feature_table, file.path(out_dir, "feature_table.csv"), row.names = FALSE)

  message(sprintf(
    "Saved %d peaks and %d aligned features to: %s",
    nrow(peak_table), nrow(feature_table), out_dir
  ))
}

# --- IPO2/Bioconductor session setup (deliberately last -- see file header) ---
#
# Uses gitlab.com/CarlBrunius/IPO2 (what the reference pipeline itself
# depends on), NOT github.com/wmoldham/IPO2 (same name, unrelated package).
# - wmoldham/IPO2: central-composite design, ~25 evaluations/iteration.
# - CarlBrunius/IPO2: nloptr Nelder-Mead, one evaluation at a time -- fewer
#   evaluations, reuses IPO's original isotope-pair scoring (0 for a
#   degenerate peak table rather than erroring).
# Install: renv::install("gitlab::CarlBrunius/IPO2") + nloptr (DESCRIPTION
# doesn't declare it despite needing it).
#
# No cross-evaluation parallelism needed here (evaluations are sequential).
# The only parallelism is within one evaluation, across the subset's files,
# inside optimPP()'s findChromPeaks() call (no BPPARAM -> session default) --
# same shape as with_bp_workers() elsewhere, so it rides the SnowParam
# default registered in R/parallel.R rather than needing its own.

# IPO2's internals call CentWaveParam()/findChromPeaks()/readMSData()
# unqualified, expecting xcms/MSnbase attached to the search path. True for
# library(IPO2); NOT true for IPO2::optimXCMS() via `::` (loads the
# namespace only). MSnbase is attached in R/spectrum_mode.R for the same
# reason; xcms needs it here too.
library(xcms)

# optimPP() reads its files itself with no centroiding step -- fine for
# centroid data, wrong otherwise. Tried writing a centroided temp mzML and
# pointing optimXCMS() at that instead, but the written files came back
# with zero MS1 spectra on reread (a real MSnbase write/read bug, not
# chased further without a live environment to debug it in). Patched here
# instead: reuse the same pickPeaks()-on-OnDiskMSnExp lazy-centroiding
# read_raw_data() already does elsewhere, applied conditionally via a
# session option run_ipo_optimization() sets before each optimXCMS() call.
patched_optimPP <- function(x0, files, optimVars, customParam, verbose = TRUE) {
  if (verbose) cat("\nParameter values:", x0, "\n")

  cwParam <- getCentWaveParams(par = x0, optimVars = optimVars, customParam = customParam)

  raw_data <- readMSData(files = files, mode = "onDisk")
  if (isTRUE(getOption("xcms_pipeline.centroid_ipo_input", FALSE))) {
    raw_data <- pickPeaks(raw_data)
  }
  xchr <- findChromPeaks(raw_data, param = cwParam)

  xset <- as(xchr, "xcmsSet")
  score <- suppressMessages(IPO2:::calcPPS(xset, "IPO")[5])

  if (verbose) cat("\nIPO score:", score, "\n")

  # Checkpoint best-so-far to disk so an interrupted run resumes from here
  # instead of the generic default (see run_ipo_optimization()). Not a true
  # resume -- nloptr doesn't expose the simplex state -- but avoids
  # restarting the whole search from scratch.
  checkpoint_path <- getOption("xcms_pipeline.ipo_checkpoint_path")
  if (!is.null(checkpoint_path) && is.finite(score)) {
    best_so_far <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path)$score else -Inf
    if (score > best_so_far) {
      saveRDS(list(cwParam = cwParam, score = score), checkpoint_path)
    }
  }

  -score
}
environment(patched_optimPP) <- asNamespace("IPO2")
assignInNamespace("optimPP", patched_optimPP, ns = "IPO2")
