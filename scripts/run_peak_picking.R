# The pipeline: IPO-optimized peak picking, retention-time alignment,
# correspondence, and gap filling, producing one aligned feature table per
# column x polarity group (e.g. RP_POS, RP_NEG, HILIC_POS never mixed).
# Modeled after github.com/MetaboComp/xcms_pipeline (pure_xcms_pipeline.R),
# stopping at the aligned feature table — no blank filtering, batch
# correction, or clustering (that's a separate later stage).
#
# Usage:
#   Rscript scripts/run_peak_picking.R <folder> [ipo_scope] [ipo_subset_size]
#
# <folder>: same one passed to generate_sample_sheet.R (or any folder with a
# hand-curated metadata/sample_sheet.xlsx meeting validate_sample_sheet()'s
# requirements).
#
# [ipo_scope] (default "global"):
# - global: one IPO optimization per column x polarity group, across all
#   its files regardless of batch.
# - batch: a separate optimization per batch (captures batch-to-batch
#   drift, at the cost of running IPO's search once per batch). Each
#   batch's params are still saved separately
#   (<folder>/output/<column>_<polarity>/<batch>/), but peak-picked results
#   are combined via xcms::c() before alignment — always exactly one
#   aligned feature table per group, regardless of scope.
#
# [ipo_subset_size] (default 4): how many files IPO2 optimizes against
# (select_ipo_subset() in R/peak_picking.R). In "global" scope spanning
# multiple batches, this becomes the max number of batches to draw one
# representative file from each. Raising it improves batch coverage but
# multiplies IPO2's per-trial cost.
#
# Output per group, under <folder>/output/<column>_<polarity>/:
#   - peaks/xdata.rds        full aligned XCMSnExp object
#   - peaks/peak_table.csv   flat per-peak table (includes gap-filled peaks)
#   - feature_table.csv      aligned feature table: one row per feature, one
#                            column per sample

source("R/sample_sheet.R")
source("R/spectrum_mode.R")
source("R/instrument_params.R")
source("R/parallel.R")
source("R/peak_picking.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript scripts/run_peak_picking.R <folder> [ipo_scope: global|batch] [ipo_subset_size]",
    call. = FALSE
  )
}
folder <- args[1]
ipo_scope <- if (length(args) >= 2) args[2] else "global"
ipo_subset_size <- if (length(args) >= 3) as.integer(args[3]) else 4

if (!ipo_scope %in% c("global", "batch")) {
  stop('ipo_scope must be "global" or "batch"', call. = FALSE)
}
if (is.na(ipo_subset_size) || ipo_subset_size < 1) {
  stop("ipo_subset_size must be a positive integer", call. = FALSE)
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

message(sprintf("IPO optimization scope: %s, subset size: %d\n", ipo_scope, ipo_subset_size))

groups <- split(sample_sheet, paste(sample_sheet$column, sample_sheet$polarity, sep = "_"))

for (group_name in names(groups)) {
  group_sheet <- groups[[group_name]]
  out_dir <- file.path(output_root, group_name)

  if (ipo_scope == "batch") {
    batch_sheets <- split(group_sheet, group_sheet$batch)
    picked_list <- list()

    for (batch_name in names(batch_sheets)) {
      batch_sheet <- batch_sheets[[batch_name]]
      message(sprintf(
        "\n=== Processing group: %s, batch: %s (%d files) ===",
        group_name, batch_name, nrow(batch_sheet)
      ))

      batch_out_dir <- file.path(out_dir, batch_name)
      centwave_param <- run_ipo_optimization(batch_sheet, batch_out_dir, ipo_subset_size)
      picked_list[[batch_name]] <- pick_peaks(
        batch_sheet$filepath, centwave_param, get_spectrum_modes(batch_sheet)
      )
    }

    xdata <- do.call(xcms::c, picked_list)
    ordered_sheet <- do.call(rbind, batch_sheets)
  } else {
    message(sprintf(
      "\n=== Processing group: %s (%d files) ===", group_name, nrow(group_sheet)
    ))

    centwave_param <- run_ipo_optimization(group_sheet, out_dir, ipo_subset_size)
    xdata <- pick_peaks(group_sheet$filepath, centwave_param, get_spectrum_modes(group_sheet))
    ordered_sheet <- group_sheet
  }

  message(sprintf("\n--- Aligning and grouping peaks for group: %s ---", group_name))
  xdata <- align_and_correspond(xdata, ordered_sheet$sample_type)
  feature_table <- build_feature_table(xdata, ordered_sheet$sample_label)
  save_peak_picking_outputs(xdata, feature_table, out_dir)
}

message("\nAll groups processed.")
