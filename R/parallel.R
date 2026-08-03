# Shared parallel-backend configuration for BiocParallel-based steps.
#
# Uses SnowParam (real separate R processes over local sockets), not
# MulticoreParam (fork-based). Two reasons: (1) MulticoreParam is known to
# break IPO2's own loop in this environment ("wrong args for environment
# subassignment") and (2) fork() doesn't exist on Windows, so MulticoreParam
# silently degrades to single-core there regardless — SnowParam is what
# actually gets real parallelism cross-platform (macOS/Linux/Windows alike).
#
# IPO2's own internal optimization loop is a separate concern: it doesn't
# accept a BPPARAM argument directly, so it's controlled by whatever backend
# is registered as the session default (see R/peak_picking.R) rather than
# by bp_workers() below.

#' Number of parallel workers to use for peak-picking/alignment steps (and
#' IPO2's loop, via the session default registered in R/peak_picking.R).
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

#' Call an xcms function with `BPPARAM = bp_workers()`, falling back to
#' calling it without BPPARAM if that specific method doesn't accept one.
#' Some xcms methods in this environment's Bioconductor devel snapshot
#' don't expose a BPPARAM argument even where the documented xcms API
#' usually does (confirmed for adjustRtime's ObiwarpParam method — "unused
#' argument (BPPARAM = ...)" at runtime, despite docs saying otherwise).
#' Rather than hardcode which ones do, this adapts at call time.
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
