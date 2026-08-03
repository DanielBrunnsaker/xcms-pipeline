# Standalone QC quality check: flags likely-faulty QC injections (missed
# injections, empty vials, degraded runs) via a quick default-parameters
# peak-picking pass, using TIC and aligned feature count vs. their group's
# peers. Run this after reviewing the sample sheet and before
# run_peak_picking.R, so a bad QC never gets picked as the "representative"
# QC for IPO2 optimization.
#
# Usage:
#   Rscript scripts/check_qc_quality.R <folder> [include_ltqc]
#
# Reads <folder>/metadata/sample_sheet.xlsx and checks sQC + ltQC files
# together (pooled) within each column x polarity group, unless
# [include_ltqc] is "false" (defaults to "true"). Pooling is deliberate: if
# the majority of one QC type is actually faulty (seen in practice), a
# median/MAD computed from that type alone is corrupted — both have a 50%
# breakdown point, so once more than half a population is bad, "normal"
# stops meaning anything. Pooling with the other QC type gives the
# threshold a better chance of being computed from a majority-good
# population. groupChromPeaks() still distinguishes sQC from ltQC internally
# via `sampleGroups = qc_sheet$sample_type`, so correspondence isn't
# affected by the pooling, only the outlier-flagging population is. If
# neither QC type has enough good files, there's nothing left to check
# against — no method can rescue a batch with no reliable reference at all.
# Set include_ltqc to "false" if ltQC shouldn't be trusted/considered at all
# for this project (e.g. an unreliable external reference matrix) — sQC
# will then be checked alone.
#
# The TIC check prefers an absolute, instrument-specific floor
# (`int_threshold` in R/instrument_params.R) over a threshold computed from
# this batch's own distribution, when available — see check_qc_quality()
# for why (median/MAD break down once a majority of a batch is actually bad).
#
# Writes the results back into the sheet as two new columns: `qc_flagged`
# (logical) and `qc_flag_reason`. Also writes a PDF report
# (<folder>/metadata/qc_quality_report.pdf) with TIC and
# aligned-feature-count barplots per group — blue bars passed, red bars
# flagged, with the pass/fail threshold marked.

source("R/sample_sheet.R")
source("R/spectrum_mode.R")
source("R/instrument_params.R")
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

sample_sheet$qc_flagged <- FALSE
sample_sheet$qc_flag_reason <- NA_character_

message(sprintf("Including ltQC in checks: %s\n", include_ltqc))

groups <- split(sample_sheet, paste(sample_sheet$column, sample_sheet$polarity, sep = "_"))
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
  sample_sheet$qc_flagged[match_idx] <- result$flagged
  sample_sheet$qc_flag_reason[match_idx] <- result$reason

  print(result[, c("filepath", "tic", "aligned_feature_count", "flagged", "reason")])

  n_flagged <- sum(result$flagged)
  if (n_flagged > 0) {
    message(sprintf("Flagged %d of %d QC file(s) as likely faulty.", n_flagged, nrow(qc_sheet)))
  } else {
    message("No QC files flagged.")
  }
}

if (length(results_by_group) > 0) {
  report_path <- file.path(folder, "metadata", "qc_quality_report.pdf")
  generate_qc_report(results_by_group, report_path)
}

backup_file(sheet_path)
writexl::write_xlsx(sample_sheet, sheet_path)
message(sprintf("\nUpdated sample sheet written to: %s", sheet_path))
