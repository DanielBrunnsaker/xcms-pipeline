# IPO2-based optimization of xcms centWave parameters, and peak picking using
# the optimized parameters, run independently per column x polarity group.
#
# Uses IPO2 from gitlab.com/CarlBrunius/IPO2 — the package the reference
# pipeline (github.com/MetaboComp/xcms_pipeline) itself depends on for its
# optimXCMS() calls — not github.com/wmoldham/IPO2 (a same-named but
# unrelated reimplementation we used earlier). The two differ fundamentally:
# - wmoldham/IPO2 repeats a ~25-point central-composite design (response
#   surface methodology) per iteration, for up to 50 iterations.
# - CarlBrunius/IPO2 wraps nloptr::nloptr() (Nelder-Mead simplex by
#   default), evaluating one candidate parameter set at a time based on the
#   previous result — far fewer evaluations to converge, and it reuses
#   IPO's original (pre-XCMS3) isotope-pair scoring, which degrades to a
#   score of 0 for a degenerate/empty peak table rather than erroring.
# Install with: renv::install("gitlab::CarlBrunius/IPO2") (also needs
# nloptr, which this package's DESCRIPTION doesn't declare despite needing
# it internally).
#
# Because evaluations are one-at-a-time (not a parallel batch), there's no
# cross-evaluation parallelism to register here. The only parallelism
# available is *within* one evaluation, across the subset's files, inside
# its optimPP()'s plain findChromPeaks() call (no BPPARAM passed, so it
# uses whatever's registered as session default) — the same shape as our
# own with_bp_workers() calls elsewhere, which already run fine under
# SnowParam. So register that here too, instead of forcing SerialParam.
BiocParallel::register(bp_workers())

# IPO2's internals (getCentWaveParams(), optimPP()) call CentWaveParam(),
# findChromPeaks(), and readMSData() unqualified, expecting xcms/MSnbase to
# be attached to the search path -- true for a normal library(IPO2) session
# (Depends: xcms gets attached automatically), but NOT true here: we only
# ever call IPO2::optimXCMS() via `::`, which loads xcms as a namespace
# only, without attaching it (Depends only auto-attaches on a real
# library()/require() of the depending package, not on `::` access).
# MSnbase already gets attached this way via R/spectrum_mode.R for the same
# reason; xcms needs the same treatment here.
library(xcms)

# IPO2's optimPP() reads its subset files itself (readMSData() +
# findChromPeaks()) with no centroiding step in between -- fine for
# centroid-mode data, but the search runs against uncentroided data
# otherwise. A temp-mzML round trip was tried first (centroid via
# read_raw_data(), write out, point optimXCMS() at the temp files) but the
# written files came back with zero MS1 spectra on reread -- a real bug
# somewhere in MSnbase's write/read cycle for this object shape, not
# something to chase blind without a live environment. Patched here
# instead: reuses the exact pickPeaks()-on-OnDiskMSnExp lazy-centroiding
# this pipeline already relies on successfully elsewhere (read_raw_data()),
# applied conditionally via a session option run_ipo_optimization() sets
# before each optimXCMS() call -- no disk round trip, no new mechanism.
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

  # Checkpoint the best result seen so far to disk, so an interrupted run
  # (killed container, crash, power loss) never loses all progress -- a
  # fresh run picks this up as its starting point instead of the generic
  # default (see run_ipo_optimization()). Not a true resume (nloptr doesn't
  # expose the simplex's internal state to save/restore), but it never
  # starts the ~100-evaluation search over from scratch either.
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
#' Within whichever tier gets used, the selection is kept to a single
#' spectrum mode (centroid preferred, since IPO2::optimXCMS() reads files
#' itself with no centroiding step -- run_ipo_optimization() centroids the
#' whole subset only if it's uniformly profile-mode) -- but only when one
#' mode alone has enough candidates to fill the request; mixed modes are
#' used as a last resort rather than handing the optimizer fewer files than
#' asked for.
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
  # run_ipo_optimization() centroids the whole subset (all-or-nothing) if
  # any file in it is profile-mode -- correct only if the subset is
  # actually uniform, since applying pickPeaks() to an already-centroided
  # file can distort it. So keep whichever tier gets used to a single
  # spectrum mode: prefer whichever mode (centroid or profile) has enough
  # candidates to fill the request on its own; only actually mix modes if
  # neither alone does.
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

  sqc_rows <- group_sheet[group_sheet$sample_type == "sQC", , drop = FALSE]
  if (nrow(sqc_rows) >= min_qc) {
    return(pick_representative(prefer_uniform_mode(sqc_rows, n), n))
  }

  ltqc_rows <- group_sheet[group_sheet$sample_type == "ltQC", , drop = FALSE]
  if (nrow(ltqc_rows) >= min_qc) {
    return(pick_representative(prefer_uniform_mode(ltqc_rows, n), n))
  }

  regular_rows <- group_sheet[group_sheet$sample_type == "Sample", , drop = FALSE]
  regular_rows <- prefer_uniform_mode(regular_rows, n)
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
#' matching config entry) is available, falling back to
#' `default_ipo2_search_space()` otherwise.
#'
#' Sets the `xcms_pipeline.centroid_ipo_input` option before calling
#' `optimXCMS()` if any subset file is profile-mode -- read by
#' `patched_optimPP()` (top of this file) to centroid in-place, since
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
#' @return An xcms::CentWaveParam with the optimized settings.
run_ipo_optimization <- function(group_sheet, out_dir, ipo_subset_size = 4, fresh = FALSE) {
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

  subset_files <- select_ipo_subset(group_sheet, n = ipo_subset_size)
  message(sprintf(
    "Running IPO2 optimization on %d file(s):\n  %s",
    length(subset_files), paste(subset_files, collapse = "\n  ")
  ))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
      message("IPO subset is profile-mode; centroiding it for the search (see patched_optimPP() at the top of this file).")
    }
  }
  old_opts <- options(
    xcms_pipeline.centroid_ipo_input = needs_centroiding,
    xcms_pipeline.ipo_checkpoint_path = checkpoint_path
  )
  on.exit(options(old_opts), add = TRUE)

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
#' retention times, groups peaks into cross-sample features, fills gaps.
#' Same steps as github.com/MetaboComp/xcms_pipeline (pure_xcms_pipeline.R):
#' adjustRtime(ObiwarpParam) -> groupChromPeaks(PeakDensityParam) ->
#' fillChromPeaks(ChromPeakAreaParam). ChromPeakAreaParam uses xcms's
#' built-in defaults (untuned, for now).
#'
#' `sampleGroups` uses our `sample_type` column (sQC/ltQC/Blank/Sample) —
#' mirrors the reference pipeline's actual usage: despite the name, their
#' "sample_group" is QC/sample TYPE classification, not biological
#' condition, and is a better fit than the sheet's manually-curated
#' (possibly unfilled) `sample_group`.
#'
#' `minFraction = 0.2` deliberately (xcms's default is 0.5; matches the
#' reference pipeline's own final `minFrac`), always, regardless of
#' `retgroup_params`: with real study samples pooled into one big "Sample"
#' group regardless of actual biological condition, a 50% presence
#' requirement silently drops any feature specific to a subset of
#' conditions. This pipeline stops at the aligned feature table on purpose
#' (blank filtering, frequency filtering, etc. are a later stage), so the
#' bar here is kept low rather than losing signal now that can't be
#' recovered later. `run_retgroup_optimization()`'s own fixed `minfrac`
#' (0.8, also matching the reference) is a separate, internal-only value
#' used just for scoring during that search -- deliberately NOT used here,
#' see that file's header comment.
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
