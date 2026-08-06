# The pipeline: IPO-optimized peak picking, retention-time alignment,
# correspondence, and gap filling, producing one aligned feature table per
# column x polarity group (e.g. RP_POS, RP_NEG, HILIC_POS never mixed).
# Modeled after github.com/MetaboComp/xcms_pipeline (pure_xcms_pipeline.R),
# stopping at the aligned feature table — no blank filtering, batch
# correction, or clustering (that's a separate later stage).
#
# Usage:
#   Rscript scripts/run_peak_picking.R <folder> [ipo_scope] [ipo_subset_size] [ipo_fresh]
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
#   (<folder>/output/<column>_<polarity>/<batch>/), and also collected
#   side by side in <folder>/output/<column>_<polarity>/
#   batch_centwave_params.csv for a quick per-batch comparison. Peak-picked
#   results are combined via c() before alignment — always exactly
#   one aligned feature table per group, regardless of scope.
#
# [ipo_subset_size] (default 4): how many files IPO2 optimizes against
# (select_ipo_subset() in R/peak_picking.R). In "global" scope spanning
# multiple batches, this becomes the max number of batches to draw one
# representative file from each. Raising it improves batch coverage but
# multiplies IPO2's per-trial cost.
#
# [ipo_fresh] (default "false"): if "true", ignore any cached IPO2 result
# and any mid-search checkpoint left by an interrupted prior run for every
# group/batch, and re-optimize from scratch. Without this, an interrupted
# run resumes each unfinished group/batch's search from its checkpoint
# (see patched_optimPP()/run_ipo_optimization() in R/peak_picking.R) rather
# than starting over.
#
# Retention-time alignment and correspondence parameters are themselves
# optimized once per group (legacy IPO::optimizeRetGroup(), see
# R/retgroup_optimization.R), cached the same way as the centWave search.
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
source("R/retgroup_optimization.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript scripts/run_peak_picking.R <folder> [ipo_scope: global|batch] [ipo_subset_size] [ipo_fresh: true|false]",
    call. = FALSE
  )
}
folder <- args[1]
ipo_scope <- if (length(args) >= 2) args[2] else "global"
ipo_subset_size <- if (length(args) >= 3) as.integer(args[3]) else 4
ipo_fresh <- if (length(args) >= 4) as.logical(args[4]) else FALSE

if (!ipo_scope %in% c("global", "batch")) {
  stop('ipo_scope must be "global" or "batch"', call. = FALSE)
}
if (is.na(ipo_subset_size) || ipo_subset_size < 1) {
  stop("ipo_subset_size must be a positive integer", call. = FALSE)
}
if (is.na(ipo_fresh)) {
  stop('ipo_fresh must be "true" or "false"', call. = FALSE)
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
      picked_list[[batch_name]] <- pick_peaks(
        batch_sheet$filepath, centwave_param, get_spectrum_modes(batch_sheet)
      )
      batch_params_rows[[batch_name]] <- cbind(
        batch = batch_name, centwave_param_to_row(centwave_param)
      )
    }

    # Consolidated, human-readable view of which centWave params were
    # actually used per batch -- each batch's own ipo_params.rds/
    # ipo_history.csv already has this individually, nested under its own
    # subfolder; this collects them side by side in one file (same idea as
    # the reference pipeline's own opt_params.csv).
    write.csv(
      do.call(rbind, batch_params_rows),
      file.path(out_dir, "batch_centwave_params.csv"),
      row.names = FALSE
    )

    # Plain c(), not xcms::c() -- xcms registers a combine method for
    # XCMSnExp on the shared "c" generic (S4 method dispatch finds it
    # regardless of which package registered it), but doesn't export a
    # symbol literally named "c" from its own namespace, so `xcms::c`
    # itself fails ("'c' is not an exported object from 'namespace:xcms'")
    # even though plain c(...) resolves correctly via dispatch.
    xdata <- do.call(c, picked_list)
    ordered_sheet <- do.call(rbind, batch_sheets)
  } else {
    message(sprintf(
      "\n=== Processing group: %s (%d files) ===", group_name, nrow(group_sheet)
    ))

    centwave_param <- run_ipo_optimization(group_sheet, out_dir, ipo_subset_size, fresh = ipo_fresh)
    xdata <- pick_peaks(group_sheet$filepath, centwave_param, get_spectrum_modes(group_sheet))
    ordered_sheet <- group_sheet
  }

  message(sprintf("\n--- Optimizing retention-time/correspondence params for group: %s ---", group_name))
  retgroup_params <- run_retgroup_optimization(xdata, out_dir, ordered_sheet$sample_type)

  message(sprintf("\n--- Aligning and grouping peaks for group: %s ---", group_name))
  xdata <- align_and_correspond(xdata, ordered_sheet$sample_type, retgroup_params)
  feature_table <- build_feature_table(xdata, ordered_sheet$sample_label)
  save_peak_picking_outputs(xdata, feature_table, out_dir)
}

message("\nAll groups processed.")
