# IPO (the legacy Bioconductor package -- NOT CarlBrunius/IPO2, which only
# optimizes centWave peak-picking and has no alignment/correspondence
# capability at all) -based optimization of retention-time alignment and
# correspondence parameters, matching Part 3/4 of the reference pipeline
# (github.com/MetaboComp/xcms_pipeline, pure_xcms_pipeline.R).
#
# IPO and IPO2 are separately-named packages (no install conflict). Install
# IPO with: BiocManager::install("IPO").
#
# Risk profile, checked before attempting this:
# - IPO::optimizeRetGroup() parallelizes via plain parallel::makeCluster(),
#   not BiocParallel -- so it should NOT hit the "bpstopOnError could not
#   be found" bug that broke IPO's own xcmsSet-based centWave optimizer for
#   us earlier in this project (that bug was specific to that other
#   function's BiocParallel usage).
# - It does still depend on legacy xcms::group()/xcms::retcor() (pre-XCMS3
#   API, operating on the xcmsSet class) -- confirmed these still exist and
#   aren't deprecated in current xcms source, but this whole integration is
#   otherwise untested in this environment. Treat as a first attempt.
#
# minfrac is fixed at 0.8 below (matching the reference's own
# minfrac_IPO), same idea as the reference: a strict, fixed bar used only
# internally so optimizeRetGroup()'s own search has a meaningful
# reproducibility signal to score alignment/grouping settings against --
# not the value actually used for the final feature table. That's set
# separately in align_and_correspond() (minFraction = 0.2, matching the
# reference's own final `minFrac`), independent of whatever this search
# finds. Only the alignment (obiwarp) and density bandwidth/bin-size
# settings from this search are actually adopted.

#' Run IPO's legacy retention-time-alignment/correspondence parameter
#' search for one column x polarity group, caching the result to disk.
#'
#' @param xdata XCMSnExp with peaks already picked (before alignment) --
#'   the full, combined object for the group (same input
#'   align_and_correspond() itself expects).
#' @param out_dir Group's output directory (e.g. output/RP_POS).
#' @param sample_types Character vector of `sample_type` values, one per
#'   file, in the same order as files in `xdata` -- used as the xcmsSet's
#'   sample class for optimizeRetGroup()'s own scoring.
#' @param n_slaves Number of parallel workers for optimizeRetGroup()'s own
#'   parallel::makeCluster() (a separate mechanism from BiocParallel/
#'   bp_workers() used elsewhere in this pipeline).
#' @return A list(obiwarp = xcms::ObiwarpParam, density_bounds = list(bw=,
#'   binSize=)) -- see align_and_correspond() for how these get applied.
run_retgroup_optimization <- function(xdata, out_dir, sample_types, n_slaves = default_worker_count()) {
  params_path <- file.path(out_dir, "retgroup_params.rds")

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

  message("Converting to legacy xcmsSet for IPO::optimizeRetGroup()...")
  xset <- as(xdata, "xcmsSet")
  xcms::sampclass(xset) <- sample_types

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
