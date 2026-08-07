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
run_retgroup_optimization <- function(xdata, out_dir, sample_types, n_slaves = default_worker_count(),
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

  qc_idx <- which(sample_types == "sQC")
  if (length(qc_idx) < min_qc) {
    ltqc_idx <- which(sample_types == "ltQC")
    if (length(ltqc_idx) >= min_qc) {
      message(sprintf(
        "Fewer than %d sQC files; using ltQC instead (not pooled with sQC) for retention-time/correspondence optimization.",
        min_qc
      ))
      qc_idx <- ltqc_idx
    } else {
      stop(
        "Not enough QC files (sQC or ltQC alone) to run retention-time/",
        "correspondence optimization for this group.",
        call. = FALSE
      )
    }
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
