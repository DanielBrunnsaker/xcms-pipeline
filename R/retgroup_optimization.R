# Legacy Bioconductor `IPO` package (not CarlBrunius/IPO2, which has no
# alignment/correspondence capability) -- optimizes retention-time
# alignment + correspondence, matching Part 3/4 of the reference pipeline.
# Separately-named package, installs alongside IPO2 with no conflict:
# BiocManager::install("IPO").
#
# Risk profile:
# - optimizeRetGroup() parallelizes via parallel::makeCluster(), not
#   BiocParallel -- shouldn't hit the bpstopOnError bug that broke IPO's
#   centWave optimizer for us earlier (that was BiocParallel-specific).
# - Still depends on legacy xcms::group()/retcor() (pre-XCMS3 API) --
#   confirmed present/not deprecated in current xcms, but otherwise
#   untested in this environment. First attempt.
#
# `minfrac` fixed at 0.8 (matches reference's `minfrac_IPO`) -- a strict
# bar used only internally so the search has a meaningful reproducibility
# signal to score against, NOT the value used for the final table (that's
# align_and_correspond()'s independent `minFraction = 0.2`). Only the
# obiwarp + density bandwidth/bin-size settings from this search get
# adopted.
#
# Runs on QC files only (sQC, falling back to ltQC if too few -- never
# pooled), not the full dataset -- matches the reference (Part 3 filters to
# just the combined sQC files from Part 1). Not just a perf shortcut: the
# scoring ("maximize features showing exactly one peak per injection of a
# pooled sample") only makes sense for QC, which really is the same sample
# repeated -- real study samples SHOULD vary, so scoring them against
# "should look identical" is meaningless.
#
# sQC/ltQC are never pooled here, unlike check_qc_quality()'s outlier
# flagging (which pools deliberately for a bigger population, while still
# grouping separately). They're chemically different matrices -- each
# internally consistent but not with each other, so pooling would look like
# noise to the algorithm. Pick one exclusively.

#' Report, per group, how many distinct batches each QC type (sQC/ltQC)
#' actually has non-flagged coverage in. Purely diagnostic -- doesn't change
#' which type `run_retgroup_optimization()` ends up using (that's still
#' decided by pooled count across the whole group, not batch coverage, see
#' this file's header comment for why pooling the two types is avoided).
#' Meant to be called once per group at the start of peak picking, so a
#' coverage gap is visible up front rather than discovered only after
#' interpreting a group's alignment results.
#'
#' @param group_sheet Sample sheet rows for one column x polarity group,
#'   every batch (this checks `qc_flagged` itself -- don't pre-filter).
report_qc_batch_coverage <- function(group_sheet) {
  qc_flagged <- group_sheet$qc_flagged
  if (is.null(qc_flagged)) {
    qc_flagged <- rep(FALSE, nrow(group_sheet))
  }
  qc_flagged[is.na(qc_flagged)] <- FALSE

  all_batches <- sort(unique(group_sheet$batch))

  for (qc_type in c("sQC", "ltQC")) {
    covered_batches <- sort(unique(group_sheet$batch[group_sheet$sample_type == qc_type & !qc_flagged]))
    missing_batches <- setdiff(all_batches, covered_batches)

    message(sprintf(
      "%s batch coverage: %d/%d batches have at least one non-flagged file.",
      qc_type, length(covered_batches), length(all_batches)
    ))
    if (length(missing_batches) > 0) {
      # message(), not warning() -- this loop is one iteration of a much
      # bigger for loop in scripts/run_peak_picking.R that processes every
      # group; a real warning() gets buffered by R and only dumped after
      # the WHOLE loop finishes (every group, all the expensive centWave/
      # alignment work included), not right here where it's actually
      # useful. message() prints immediately instead.
      message(sprintf(
        "WARNING: %s has no non-flagged files at all for %d batch(es): %s -- if %s ends up used for retention-time/correspondence optimization, these batches won't be represented in that search (still aligned using the group's shared params, just not part of what tuned them).",
        qc_type, length(missing_batches), paste(missing_batches, collapse = ", "), qc_type
      ))
    }
  }
}

#' Select which QC files feed the retention-time/correspondence
#' optimization search below: either the automatic sQC-first/ltQC-fallback
#' pick, or a caller-forced type -- sQC and ltQC are still never pooled
#' together either way (see this file's header comment for why).
#'
#' @param sample_types Character vector of `sample_type` values.
#' @param qc_flagged Logical vector, same length as `sample_types` --
#'   excluded from selection. NULL or NA per-element means "not flagged".
#' @param qc_type Force `"sQC"` or `"ltQC"` (case-insensitive), bypassing
#'   the automatic fallback below. NULL (default) keeps the automatic
#'   behavior. Forcing doesn't waive the `min_qc` floor -- still errors if
#'   the forced type doesn't have enough non-flagged files.
#' @param min_qc Minimum non-flagged files required for a type to count as
#'   usable.
#' @return list(idx = integer indices into `sample_types`, type = "sQC" or
#'   "ltQC" -- whichever was actually selected, forced = TRUE if `qc_type`
#'   drove the choice rather than the automatic fallback).
select_retgroup_qc_idx <- function(sample_types, qc_flagged = NULL, qc_type = NULL, min_qc = 2) {
  if (is.null(qc_flagged)) {
    qc_flagged <- rep(FALSE, length(sample_types))
  }
  qc_flagged[is.na(qc_flagged)] <- FALSE

  if (!is.null(qc_type)) {
    canonical <- c(sqc = "sQC", ltqc = "ltQC")
    key <- tolower(qc_type)
    if (!key %in% names(canonical)) {
      stop(sprintf('qc_type must be "sQC" or "ltQC" (got "%s").', qc_type), call. = FALSE)
    }
    forced_type <- canonical[[key]]
    idx <- which(sample_types == forced_type & !qc_flagged)
    if (length(idx) < min_qc) {
      stop(sprintf(
        "qc_type forced to %s, but only %d non-flagged %s file(s) available (need >= %d) for retention-time/correspondence optimization.",
        forced_type, length(idx), forced_type, min_qc
      ), call. = FALSE)
    }
    return(list(idx = idx, type = forced_type, forced = TRUE))
  }

  idx <- which(sample_types == "sQC" & !qc_flagged)
  if (length(idx) >= min_qc) {
    return(list(idx = idx, type = "sQC", forced = FALSE))
  }

  ltqc_idx <- which(sample_types == "ltQC" & !qc_flagged)
  if (length(ltqc_idx) >= min_qc) {
    return(list(idx = ltqc_idx, type = "ltQC", forced = FALSE))
  }

  stop(
    "Not enough non-flagged QC files (sQC or ltQC alone) to run retention-time/",
    "correspondence optimization for this group. Review qc_flagged in the sample sheet.",
    call. = FALSE
  )
}

#' Run IPO's legacy retention-time-alignment/correspondence parameter
#' search for one column x polarity group, caching the result to disk.
#'
#' @param xdata XCMSnExp with peaks already picked (before alignment) --
#'   the full, combined object for the group (same input
#'   align_and_correspond() itself expects). Filtered down to QC files
#'   only internally -- see this file's header comment for why.
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @param sample_types Character vector of `sample_type` values, one per
#'   file, in the same order as files in `xdata`.
#' @param qc_flagged Logical vector, same length/order as `sample_types` --
#'   `check_qc_quality()`'s `qc_flagged` column, if available. Excluded from
#'   the sQC/ltQC selection below (same reasoning as
#'   `select_ipo_subset()`'s: `optimizeRetGroup()`'s scoring assumes "this
#'   is the same sample repeated, should look identical", which a
#'   known-bad QC violates). NULL (default) or NA per-element means "not
#'   flagged" -- matches the sheet's own opt-out convention
#'   (`get_included()`), so this stays a no-op for a hand-curated sheet
#'   that never went through check_qc_quality().
#' @param qc_type Force `"sQC"` or `"ltQC"` (case-insensitive) for this
#'   search, bypassing the automatic sQC-first/ltQC-fallback selection --
#'   see `select_retgroup_qc_idx()`. NULL (default) keeps the automatic
#'   behavior.
#' @param n_slaves Number of parallel workers for optimizeRetGroup()'s own
#'   parallel::makeCluster() (a separate mechanism from BiocParallel/
#'   bp_workers() used elsewhere in this pipeline).
#' @param min_qc Minimum QC files (of whichever single type ends up used)
#'   required for the search to run at all -- also the bar for sQC to count
#'   as "enough" on its own before falling back to ltQC instead. Same
#'   meaning/convention as select_ipo_subset()'s.
#' @param fresh If TRUE, ignore (and delete) any cached result from a prior
#'   run and re-optimize from scratch -- same flag/meaning as
#'   run_ipo_optimization()'s.
#' @return A list(obiwarp = xcms::ObiwarpParam, density_bounds = list(bw=,
#'   binSize=)) -- see align_and_correspond() for how these get applied.
run_retgroup_optimization <- function(xdata, out_dir, sample_types, qc_flagged = NULL, qc_type = NULL,
                                       n_slaves = default_worker_count(),
                                       min_qc = 2, fresh = FALSE) {
  params_path <- file.path(out_dir, "retgroup_params.rds")

  if (fresh) {
    if (file.exists(params_path)) {
      message("fresh = TRUE: ignoring cached retention-time/correspondence result, re-optimizing from scratch.")
    }
    unlink(params_path)
  }

  if (file.exists(params_path)) {
    message(sprintf("Using cached retention-time/correspondence params: %s", params_path))
    return(readRDS(params_path))
  }

  if (!requireNamespace("IPO", quietly = TRUE)) {
    stop(
      "IPO package not installed (needed for run_retgroup_optimization()) -- ",
      "see docker/install_packages.R.",
      call. = FALSE
    )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  qc_selection <- select_retgroup_qc_idx(sample_types, qc_flagged, qc_type, min_qc)
  qc_idx <- qc_selection$idx
  if (qc_selection$forced) {
    message(sprintf("qc_type forced to %s for retention-time/correspondence optimization.", qc_selection$type))
  } else if (qc_selection$type == "ltQC") {
    message(sprintf(
      "Fewer than %d non-flagged sQC files; using ltQC instead (not pooled with sQC) for retention-time/correspondence optimization.",
      min_qc
    ))
  }

  message(sprintf(
    "Restricting retention-time/correspondence search to %d QC file(s) (out of %d total)...",
    length(qc_idx), length(sample_types)
  ))
  qc_xdata <- xcms::filterFile(xdata, file = qc_idx)
  qc_sample_types <- sample_types[qc_idx]

  message("Converting to legacy xcmsSet for IPO::optimizeRetGroup()...")
  xset <- as(qc_xdata, "xcmsSet")
  xcms::sampclass(xset) <- qc_sample_types

  retgroup_starting_params <- IPO::getDefaultRetGroupStartingParams()
  retgroup_starting_params$minfrac <- 0.8

  message("Running IPO::optimizeRetGroup() (retention time + correspondence search)...")
  result <- IPO::optimizeRetGroup(
    xset = xset,
    params = retgroup_starting_params,
    nSlaves = n_slaves,
    subdir = NULL,
    plot = FALSE
  )

  best <- result$best_settings

  retgroup_params <- list(
    obiwarp = xcms::ObiwarpParam(
      binSize = best$profStep,
      response = best$response,
      gapInit = best$gapInit,
      gapExtend = best$gapExtend
    ),
    density_bounds = list(bw = best$bw, binSize = best$mzwid)
  )

  saveRDS(retgroup_params, params_path)
  saveRDS(result, file.path(out_dir, "retgroup_history.rds"))

  message(sprintf("Saved optimized retention-time/correspondence params to: %s", params_path))

  retgroup_params
}
