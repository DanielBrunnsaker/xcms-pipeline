# Install project dependencies for the Docker image. Tries the exact
# versions pinned in renv.lock first; Bioconductor's devel repo doesn't
# retain old package builds the way CRAN archives releases, and it churns
# fast, so a specific devel-era version can become unavailable within days
# of being pinned.
#
# When restore() fails for that reason:
# 1. Patch the specific packages known to hit this (Bioconductor core +
#    IPO2) with whatever devel currently serves.
# 2. Retry restore() — a single failed attempt can abort before reaching
#    later-queued packages, silently leaving them uninstalled without
#    naming them in the error (this bit us once: writexl went missing from
#    a build with no mention of it in the failure output).
# 3. Audit the *entire* lockfile against what's actually installed and
#    install anything still missing, regardless of why it was skipped —
#    a final safety net so a partial restore() never silently ships short.
# Re-snapshots at the end so this image's own lockfile reflects reality.

try_restore <- function() {
  tryCatch(
    {
      renv::restore()
      TRUE
    },
    error = function(e) {
      message("renv::restore() failed: ", conditionMessage(e))
      FALSE
    }
  )
}

restore_ok <- try_restore()

if (!restore_ok) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  core_bioc_pkgs <- c("BiocParallel", "mzR", "MSnbase", "xcms")
  missing_bioc <- setdiff(core_bioc_pkgs, rownames(installed.packages()))
  if (length(missing_bioc) > 0) {
    message(sprintf("Installing current devel versions of: %s", paste(missing_bioc, collapse = ", ")))
    BiocManager::install(missing_bioc, update = TRUE, ask = FALSE)
  }

  if (!"IPO2" %in% rownames(installed.packages())) {
    renv::install("wmoldham/IPO2")
  }

  # Pick up anything else restore() hadn't reached before aborting.
  try_restore()

  renv::snapshot(prompt = FALSE)
}

# Final audit: install anything the lockfile expects that still isn't
# present, regardless of why it was skipped.
locked_pkgs <- names(renv::lockfile_read("renv.lock")$Packages)
missing_pkgs <- setdiff(locked_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message(sprintf("Still missing after restore/fallback: %s", paste(missing_pkgs, collapse = ", ")))
  renv::install(missing_pkgs)
  renv::snapshot(prompt = FALSE)
}
