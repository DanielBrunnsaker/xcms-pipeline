# The pipeline: IPO-optimized peak picking, retention-time alignment,
# correspondence, and gap filling -> one aligned feature table per column x
# polarity group (RP_POS, RP_NEG, HILIC_POS never mixed). Modeled after
# github.com/MetaboComp/xcms_pipeline, stopping at the aligned feature
# table (blank filtering/batch correction/clustering are a later stage).
# Retention-time/correspondence params are also optimized once per group
# (legacy IPO::optimizeRetGroup(), see R/retgroup_optimization.R), cached
# like the centWave search.
#
# Usage:
#   Rscript scripts/run_peak_picking.R <folder> [ipo_scope] [ipo_subset_size] [ipo_fresh] [retgroup_qc_type]
#
# <folder>: same one passed to generate_sample_sheet.R (or any folder with a
# hand-curated metadata/sample_sheet.xlsx meeting validate_sample_sheet()'s
# requirements).
#
# [ipo_scope] (default "global"):
# - global: one IPO optimization per group, across all its files.
# - batch: separate optimization per batch (captures batch-to-batch drift,
#   costs one search per batch). Params saved per batch
#   (<folder>/output/<column>_<polarity>/<batch>/) and collected side by
#   side in .../batch_centwave_params.csv. Peak-picked results combine via
#   c() before alignment -- always one feature table per group either way.
#
# [ipo_subset_size] (default 4): files IPO2 optimizes against
# (select_ipo_subset()). In "global" scope with multiple batches, this is
# the max batches to draw one representative file from each. Larger =
# better batch coverage, higher per-trial cost.
#
# [ipo_fresh] (default "false"): "true" ignores any cached result/
# checkpoint for every group/batch and re-optimizes from scratch, and also
# discards any cached picked-peaks result (picked_peaks.rds -- see below).
# Default resumes an interrupted run from its checkpoint/cache instead of
# restarting.
#
# [retgroup_qc_type] (default "auto"): force which QC type
# (sQC|ltQC, case-insensitive) feeds retention-time/correspondence
# optimization (R/retgroup_optimization.R) for every group, instead of the
# automatic sQC-first/ltQC-fallback pick. "auto" keeps that automatic
# behavior. Still errors per group if the forced type doesn't have enough
# non-flagged files -- see the "QC batch coverage" printout at the start of
# each group for whether sQC/ltQC actually cover every batch.
#
# Output per group, under <folder>/output/<column>_<polarity>/:
#   - picked_peaks.rds       cached pick_peaks() result (per batch, under
#                            <batch>/, in batch scope) -- an interrupted run
#                            resumes without re-running findChromPeaks() on
#                            an already-picked batch/group. Auto-invalidated
#                            if files/params changed since it was written.
#   - peaks/xdata.rds        full aligned XCMSnExp object
#   - peaks/peak_table.csv   flat per-peak table (includes gap-filled peaks)
#   - feature_table.csv      aligned feature table: one row per feature, one
#                            column per sample

source("R/sample_sheet.R")
source("R/spectrum_mode.R")
source("R/instrument_params.R")
source("R/parallel.R")
source("R/peak_picking.R")
source("R/retgroup_optimization.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    paste(
      "Usage: Rscript scripts/run_peak_picking.R <folder> [ipo_scope: global|batch]",
      "[ipo_subset_size] [ipo_fresh: true|false] [retgroup_qc_type: auto|sQC|ltQC]"
    ),
    call. = FALSE
  )
}
folder <- args[1]
ipo_scope <- if (length(args) >= 2) args[2] else "global"
ipo_subset_size <- if (length(args) >= 3) as.integer(args[3]) else 4
ipo_fresh <- if (length(args) >= 4) as.logical(args[4]) else FALSE
retgroup_qc_type <- if (length(args) >= 5) args[5] else "auto"

if (!ipo_scope %in% c("global", "batch")) {
  stop('ipo_scope must be "global" or "batch"', call. = FALSE)
}
if (is.na(ipo_subset_size) || ipo_subset_size < 1) {
  stop("ipo_subset_size must be a positive integer", call. = FALSE)
}
if (is.na(ipo_fresh)) {
  stop('ipo_fresh must be "true" or "false"', call. = FALSE)
}
if (!tolower(retgroup_qc_type) %in% c("auto", "sqc", "ltqc")) {
  stop('retgroup_qc_type must be "auto", "sQC", or "ltQC"', call. = FALSE)
}

sheet_path <- file.path(folder, "metadata", "sample_sheet.xlsx")
output_root <- file.path(folder, "output")

if (!file.exists(sheet_path)) {
  stop(sprintf(
    "Sample sheet not found: %s\n(Run generate_sample_sheet.R first, or place a reviewed sheet there by hand.)",
    sheet_path
  ), call. = FALSE)
}

sample_sheet <- as.data.frame(readxl::read_excel(sheet_path))
validate_sample_sheet(sample_sheet)

# Snapshot the sheet exactly as used for this run -- so a feature_table.csv
# from months ago can be traced back to precisely what it was built from,
# even if the live sample_sheet.xlsx gets edited and rerun later. One
# snapshot per run (not per group): every group in this run shares the same
# sheet state.
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
sheet_snapshot_path <- file.path(
  output_root, sprintf("sample_sheet_snapshot_%s.xlsx", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
file.copy(sheet_path, sheet_snapshot_path)
message(sprintf("Sample sheet snapshot for this run: %s", sheet_snapshot_path))

included <- get_included(sample_sheet)
if (any(!included)) {
  message(sprintf("Excluding %d file(s) marked include=FALSE.", sum(!included)))
}
sample_sheet <- sample_sheet[included, , drop = FALSE]

message(sprintf(
  "IPO optimization scope: %s, subset size: %d, fresh: %s\n", ipo_scope, ipo_subset_size, ipo_fresh
))

groups <- split(sample_sheet, paste(sample_sheet$column, sample_sheet$polarity, sep = "_"))

for (group_name in names(groups)) {
  group_sheet <- groups[[group_name]]
  out_dir <- file.path(output_root, group_name)

  message(sprintf("\n=== QC batch coverage for group: %s ===", group_name))
  report_qc_batch_coverage(group_sheet)

  # Same snapshot as the root-level one (copied, not regenerated, so every
  # copy is byte-identical) -- keeps this group's folder self-contained and
  # traceable even if it's later copied/archived on its own.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(sheet_snapshot_path, file.path(out_dir, basename(sheet_snapshot_path)))

  if (ipo_scope == "batch") {
    batch_sheets <- split(group_sheet, group_sheet$batch)
    picked_list <- list()
    batch_params_rows <- list()

    for (batch_name in names(batch_sheets)) {
      batch_sheet <- batch_sheets[[batch_name]]
      message(sprintf(
        "\n=== Processing group: %s, batch: %s (%d files) ===",
        group_name, batch_name, nrow(batch_sheet)
      ))

      batch_out_dir <- file.path(out_dir, batch_name)
      centwave_param <- run_ipo_optimization(batch_sheet, batch_out_dir, ipo_subset_size, fresh = ipo_fresh)
      picked_list[[batch_name]] <- pick_peaks_cached(
        batch_sheet$filepath, centwave_param, get_spectrum_modes(batch_sheet),
        cache_path = file.path(batch_out_dir, "picked_peaks.rds"), fresh = ipo_fresh
      )
      batch_params_rows[[batch_name]] <- cbind(
        batch = batch_name, centwave_param_to_row(centwave_param)
      )
    }

    # Side-by-side view of every batch's actual centWave params (each
    # batch's own ipo_params.rds/ipo_history.csv already has this
    # individually) -- same idea as the reference's opt_params.csv.
    write.csv(
      do.call(rbind, batch_params_rows),
      file.path(out_dir, "batch_centwave_params.csv"),
      row.names = FALSE
    )

    # Plain c(), not xcms::c() -- xcms registers a combine method on the
    # shared "c" generic (dispatch finds it regardless of package) but
    # doesn't export a symbol named "c" from its own namespace, so
    # `xcms::c` itself fails even though plain c(...) works.
    xdata <- do.call(c, picked_list)
    ordered_sheet <- do.call(rbind, batch_sheets)
  } else {
    message(sprintf(
      "\n=== Processing group: %s (%d files) ===", group_name, nrow(group_sheet)
    ))

    centwave_param <- run_ipo_optimization(group_sheet, out_dir, ipo_subset_size, fresh = ipo_fresh)
    xdata <- pick_peaks_cached(
      group_sheet$filepath, centwave_param, get_spectrum_modes(group_sheet),
      cache_path = file.path(out_dir, "picked_peaks.rds"), fresh = ipo_fresh
    )
    ordered_sheet <- group_sheet
  }

  message(sprintf("\n--- Optimizing retention-time/correspondence params for group: %s ---", group_name))
  retgroup_params <- run_retgroup_optimization(
    xdata, out_dir, ordered_sheet$sample_type,
    batches = ordered_sheet$batch, injection_order = ordered_sheet$injection_order,
    qc_flagged = ordered_sheet$qc_flagged,
    qc_type = if (tolower(retgroup_qc_type) == "auto") NULL else retgroup_qc_type,
    fresh = ipo_fresh
  )

  message(sprintf("\n--- Aligning and grouping peaks for group: %s ---", group_name))
  xdata <- align_and_correspond(xdata, ordered_sheet$sample_type, retgroup_params)
  feature_table <- build_feature_table(xdata, ordered_sheet$sample_label)
  save_peak_picking_outputs(xdata, feature_table, out_dir, ordered_sheet$sample_label)
}

message("\nAll groups processed.")
