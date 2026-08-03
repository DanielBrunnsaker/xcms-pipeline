# Install project dependencies for the Docker image. Tries the exact
# versions pinned in renv.lock first; Bioconductor's devel repo doesn't
# retain old package builds the way CRAN archives releases, and it churns
# fast, so a specific devel-era version can become unavailable within days
# of being pinned. When restore() fails for that reason, fall back to
# installing the current devel version of just the packages that failed,
# then re-snapshot so this image's own lockfile reflects what's actually
# inside it (rather than a lockfile that no longer matches reality).

restore_ok <- tryCatch(
  {
    renv::restore()
    TRUE
  },
  error = function(e) {
    message("renv::restore() failed, falling back for unavailable packages: ", conditionMessage(e))
    FALSE
  }
)

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

  renv::snapshot(prompt = FALSE)
}
