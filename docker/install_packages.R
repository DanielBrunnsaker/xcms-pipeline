# Install project dependencies for the Docker image. Tries the exact
# versions pinned in renv.lock first; Bioconductor's devel repo doesn't
# retain old package builds the way CRAN archives releases, and it churns
# fast, so a specific devel-era version can become unavailable within days
# of being pinned.
#
# The original lockfile's package list is read ONCE, right at the start,
# before anything below can touch it. Every check against "what should be
# installed" uses that captured list, never a re-read of renv.lock —
# renv::snapshot() rewrites the lockfile to match whatever's *currently*
# installed, so re-reading after a snapshot would silently agree with
# reality instead of catching what's still missing (this already happened
# once: writexl went missing from a build, an intermediate snapshot() quietly
# dropped it from the lockfile, and a later audit against that same
# already-rewritten file saw nothing wrong).
#
# When restore() fails (e.g. an unavailable devel-era version):
# 1. Patch the specific packages known to hit this (Bioconductor core +
#    IPO2, from gitlab.com/CarlBrunius/IPO2) with whatever devel currently
#    serves.
# 2. Retry restore() — a single failed attempt can abort before reaching
#    later-queued packages, silently leaving them uninstalled without
#    naming them in the error.
# 3. Audit the captured original package list against what's actually
#    installed and install anything still missing, regardless of why it
#    was skipped — a final safety net so a partial restore() never
#    silently ships short.
# Snapshots exactly once, at the very end, once everything is confirmed
# installed.

original_locked_pkgs <- names(renv::lockfile_read("renv.lock")$Packages)

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
  missing_bioc <- intersect(core_bioc_pkgs, setdiff(original_locked_pkgs, rownames(installed.packages())))
  if (length(missing_bioc) > 0) {
    message(sprintf("Installing current devel versions of: %s", paste(missing_bioc, collapse = ", ")))
    BiocManager::install(missing_bioc, update = TRUE, ask = FALSE)
  }

  if ("IPO2" %in% original_locked_pkgs && !"IPO2" %in% rownames(installed.packages())) {
    # gitlab.com/CarlBrunius/IPO2 -- the package the reference pipeline
    # (MetaboComp/xcms_pipeline) itself depends on for optimXCMS(). Its
    # DESCRIPTION doesn't declare nloptr despite calling nloptr::nloptr()
    # internally, so install it explicitly too.
    renv::install("nloptr")
    renv::install("gitlab::CarlBrunius/IPO2")
  }

  # Pick up anything else restore() hadn't reached before aborting.
  try_restore()
}

# Final audit against the ORIGINAL package list captured above.
missing_pkgs <- setdiff(original_locked_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message(sprintf("Still missing after restore/fallback: %s", paste(missing_pkgs, collapse = ", ")))
  renv::install(missing_pkgs)
}

# Snapshot exactly once, now that everything should actually be installed.
renv::snapshot(prompt = FALSE)
