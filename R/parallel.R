# Shared parallel-backend configuration for BiocParallel-based steps.
#
# Uses SnowParam (separate processes over local sockets), not MulticoreParam
# (fork-based): MulticoreParam previously hit a BiocParallel bug ("wrong
# args for environment subassignment"), and fork() doesn't exist on
# Windows anyway -- SnowParam is what actually works cross-platform.
#
# IPO2::optimXCMS()'s internal findChromPeaks() call, and some xcms methods
# in this environment (e.g. adjustRtime()'s ObiwarpParam method), don't take
# a BPPARAM argument, so they're controlled by the session default backend
# instead -- registered at the bottom of this file for that reason.

#' Number of parallel workers to use for peak-picking/alignment steps (and
#' any xcms/IPO2 call that falls back to the session default backend
#' registered at the bottom of this file).
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
#' default registered below.
#'
#' @param workers Worker count. Defaults to `default_worker_count()`;
#'   callers processing the full (not subsampled) file set against
#'   memory-heavy steps -- see `align_and_correspond()` -- may want fewer.
bp_workers <- function(workers = default_worker_count()) {
  BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
}

#' Call an xcms function with `BPPARAM = bpparam`, falling back to no
#' BPPARAM if that method doesn't accept one -- some xcms methods in this
#' Bioconductor devel snapshot don't take BPPARAM despite the docs (e.g.
#' adjustRtime's ObiwarpParam method). Adapts at call time rather than
#' hardcoding which ones do.
#'
#' The fallback path does NOT run single-threaded, despite what it used to
#' say here -- a method that rejects BPPARAM entirely still parallelizes,
#' just via whatever's registered as the *session-wide default* backend
#' (see R/parallel.R's own registration at the bottom of this file) rather
#' than `bpparam`. A caller that needs to control worker count even for
#' that fallback case (e.g. `align_and_correspond()`, for memory reasons)
#' has to temporarily re-register the default around the call instead of
#' relying on this function's `bpparam` argument alone.
#'
#' @param fn The xcms function to call (e.g. xcms::adjustRtime).
#' @param ... Arguments to pass (not including BPPARAM).
#' @param bpparam BiocParallel backend to request. Defaults to
#'   `bp_workers()` (the pipeline's normal worker count).
with_bp_workers <- function(fn, ..., bpparam = bp_workers()) {
  tryCatch(
    fn(..., BPPARAM = bpparam),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("unused argument", msg, fixed = TRUE) && grepl("BPPARAM", msg, fixed = TRUE)) {
        warning(
          "This xcms call doesn't accept BPPARAM in this environment; falling back to the ",
          "session-registered default backend instead of the worker count requested here.",
          call. = FALSE
        )
        fn(...)
      } else {
        stop(e)
      }
    }
  )
}

# Register SnowParam as the session-wide BiocParallel default too, not just
# the explicit BPPARAM passed by with_bp_workers() above. Several xcms calls
# in this Bioconductor snapshot don't accept BPPARAM at all (e.g.
# adjustRtime()'s ObiwarpParam method -- see with_bp_workers()) and silently
# fall back to whatever's registered as default; same for IPO2::optimXCMS()'s
# internal findChromPeaks() call. Without this, that default is whatever
# BiocParallel auto-selects for the platform -- MulticoreParam on
# macOS/Linux, which is the exact backend that hits "wrong args for
# environment subassignment" (why bp_workers() uses SnowParam explicitly in
# the first place). Registering here, in the file every pipeline script
# sources, means every call path gets the known-working backend, not just
# the ones that happen to source R/peak_picking.R too.
#
# Deliberately positioned LAST, after every function above -- in a sandbox
# that can't bind a SnowParam socket, this is the one line that errors (see
# tests/testthat/helper-source.R's source_tolerant()), and every function
# above it must still be defined regardless.
BiocParallel::register(bp_workers())
