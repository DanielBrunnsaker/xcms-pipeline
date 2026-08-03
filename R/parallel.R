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

#' Number of parallel workers to use for peak-picking/alignment steps,
#' leaving some headroom for the OS/other processes rather than claiming
#' every core.
#'
#' @param headroom Cores to leave free.
default_worker_count <- function(headroom = 2) {
  max(1, parallel::detectCores() - headroom)
}

#' Build the BiocParallel backend for peak-picking/alignment steps —
#' separate from IPO2's own loop, which is controlled by the session-wide
#' default registered in R/peak_picking.R.
bp_workers <- function() {
  BiocParallel::SnowParam(workers = default_worker_count(), progressbar = TRUE)
}
