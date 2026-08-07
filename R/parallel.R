# Shared parallel-backend configuration for BiocParallel-based steps.
#
# Uses SnowParam (separate processes over local sockets), not MulticoreParam
# (fork-based): MulticoreParam previously hit a BiocParallel bug ("wrong
# args for environment subassignment"), and fork() doesn't exist on
# Windows anyway -- SnowParam is what actually works cross-platform.
#
# IPO2::optimXCMS()'s internal findChromPeaks() doesn't take a BPPARAM
# argument, so it's controlled by the session default backend instead (see
# R/peak_picking.R, which registers bp_workers() for that reason).

#' Number of parallel workers to use for peak-picking/alignment steps (and
#' IPO2::optimXCMS()'s internal findChromPeaks() call, via the session
#' default registered in R/peak_picking.R).
#'
#' Reads the `XCMS_PIPELINE_CORES` env var if set (e.g.
#' `docker run -e XCMS_PIPELINE_CORES=8 ...`); otherwise defaults to all
#' detected cores minus `headroom`, leaving some free for the OS/other
#' processes rather than claiming every core.
#'
#' @param headroom Cores to leave free when no override is set.
default_worker_count <- function(headroom = 2) {
  override <- Sys.getenv("XCMS_PIPELINE_CORES", unset = NA)
  if (!is.na(override)) {
    n <- suppressWarnings(as.integer(override))
    if (!is.na(n) && n >= 1) {
      return(n)
    }
    warning(sprintf(
      "XCMS_PIPELINE_CORES=%s is not a valid positive integer; ignoring.", override
    ), call. = FALSE)
  }

  max(1, parallel::detectCores() - headroom)
}

#' Build the BiocParallel backend for peak-picking/alignment steps —
#' separate from IPO2's own loop, which is controlled by the session-wide
#' default registered in R/peak_picking.R.
bp_workers <- function() {
  BiocParallel::SnowParam(workers = default_worker_count(), progressbar = TRUE)
}

#' Call an xcms function with `BPPARAM = bp_workers()`, falling back to no
#' BPPARAM if that method doesn't accept one -- some xcms methods in this
#' Bioconductor devel snapshot don't take BPPARAM despite the docs (e.g.
#' adjustRtime's ObiwarpParam method). Adapts at call time rather than
#' hardcoding which ones do.
#'
#' @param fn The xcms function to call (e.g. xcms::adjustRtime).
#' @param ... Arguments to pass (not including BPPARAM).
with_bp_workers <- function(fn, ...) {
  tryCatch(
    fn(..., BPPARAM = bp_workers()),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("unused argument", msg, fixed = TRUE) && grepl("BPPARAM", msg, fixed = TRUE)) {
        warning(
          "This xcms call doesn't accept BPPARAM in this environment; running without it (single-threaded).",
          call. = FALSE
        )
        fn(...)
      } else {
        stop(e)
      }
    }
  )
}
