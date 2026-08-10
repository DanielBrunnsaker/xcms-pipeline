# Standalone QC quality check: flags likely-faulty QC injections (missed
# injections, empty vials, degraded runs) via raw peak count and aligned
# feature count vs. group peers. Run after reviewing the sample sheet,
# before run_peak_picking.R, so a bad QC never becomes IPO's
# "representative" QC.
#
# Usage:
#   Rscript scripts/check_qc_quality.R <folder> [include_ltqc]
#
# Checks sQC and ltQC per column x polarity group, unless [include_ltqc] is
# "false" (default "true").
# - Raw peak count (per file, before alignment) is the primary, stricter
#   signal -- far less sensitive to concentration than TIC (not used at
#   all here), so it isolates "did this injection clearly fail" from
#   "is the signal magnitude a bit different." Aligned feature count is a
#   secondary, much more lenient check: even a technically-fine injection
#   is a poor normalization/correction anchor if its peaks don't correspond
#   well with its peers.
# - sQC and ltQC are never pooled for either metric's threshold -- different
#   matrices can have systematically different baselines, so each type is
#   judged against its own median instead of a shared one. groupChromPeaks()
#   also groups by sample_type for the same reason. Entirely data-driven
#   (no fixed absolute floor) -- see check_qc_quality() for why.
# - include_ltqc "false": check sQC alone if ltQC isn't trusted.
# - A type with too few good files just doesn't get its own threshold.
#
# Writes `qc_flagged_global`/`qc_flagged_batch`/`qc_flagged` (either --
# what select_ipo_subset() excludes on) and `qc_flag_reason` back into the
# sheet, plus an interactive HTML report
# (<folder>/metadata/qc_quality_report.html).

source("R/sample_sheet.R")
source("R/spectrum_mode.R")
source("R/instrument_params.R")
source("R/parallel.R")
source("R/qc_quality.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/check_qc_quality.R <folder> [include_ltqc: true|false]", call. = FALSE)
}
folder <- args[1]
include_ltqc <- if (length(args) >= 2) as.logical(args[2]) else TRUE
if (is.na(include_ltqc)) {
  stop('include_ltqc must be "true" or "false"', call. = FALSE)
}

sheet_path <- file.path(folder, "metadata", "sample_sheet.xlsx")
if (!file.exists(sheet_path)) {
  stop(sprintf(
    "Sample sheet not found: %s\n(Run generate_sample_sheet.R first, or place a reviewed sheet there by hand.)",
    sheet_path
  ), call. = FALSE)
}

sample_sheet <- as.data.frame(readxl::read_excel(sheet_path))
validate_sample_sheet(sample_sheet)

sample_sheet$qc_flagged_global <- FALSE
sample_sheet$qc_flagged_batch <- FALSE
sample_sheet$qc_flagged <- FALSE
sample_sheet$qc_flag_reason <- NA_character_

# Excluded (include=FALSE) rows are skipped for checking, but sample_sheet
# itself stays whole — it gets written back in full at the end, so
# filtering it directly here would silently drop those rows from the file
# rather than just skip processing them.
included <- get_included(sample_sheet)
if (any(!included)) {
  message(sprintf("Excluding %d file(s) marked include=FALSE from QC checks.", sum(!included)))
}
checkable_sheet <- sample_sheet[included, , drop = FALSE]

message(sprintf("Including ltQC in checks: %s\n", include_ltqc))

groups <- split(checkable_sheet, paste(checkable_sheet$column, checkable_sheet$polarity, sep = "_"))
results_by_group <- list()

for (group_name in names(groups)) {
  group_sheet <- groups[[group_name]]
  qc_sheet <- if (include_ltqc) {
    group_sheet[group_sheet$is_qc, , drop = FALSE]
  } else {
    group_sheet[group_sheet$sample_type == "sQC", , drop = FALSE]
  }

  if (nrow(qc_sheet) < 2) {
    message(sprintf(
      "\n=== Group: %s: fewer than 2 QC files, nothing to compare, skipping ===",
      group_name
    ))
    next
  }

  message(sprintf("\n=== Checking QC quality for group: %s (%d QC files) ===", group_name, nrow(qc_sheet)))
  result <- check_qc_quality(qc_sheet)
  results_by_group[[group_name]] <- result

  match_idx <- match(result$filepath, sample_sheet$filepath)
  sample_sheet$qc_flagged_global[match_idx] <- result$flagged_global
  sample_sheet$qc_flagged_batch[match_idx] <- result$flagged_batch
  sample_sheet$qc_flagged[match_idx] <- result$flagged
  sample_sheet$qc_flag_reason[match_idx] <- result$reason

  print(result[, c(
    "filepath", "raw_peak_count", "aligned_feature_count", "flagged_global", "flagged_batch", "reason"
  )])

  n_flagged <- sum(result$flagged)
  if (n_flagged > 0) {
    message(sprintf("Flagged %d of %d QC file(s) as likely faulty.", n_flagged, nrow(qc_sheet)))
  } else {
    message("No QC files flagged.")
  }
}

if (length(results_by_group) > 0) {
  report_path <- file.path(folder, "metadata", "qc_quality_report.html")
  generate_qc_report(results_by_group, report_path)
}

backup_file(sheet_path)
writexl::write_xlsx(sample_sheet, sheet_path)
message(sprintf("\nUpdated sample sheet written to: %s", sheet_path))
