# Stage 2: IPO-optimized peak picking for a project, run independently per
# column x polarity group (e.g. RP_POS, RP_NEG, HILIC_POS never mixed).
#
# Usage:
#   Rscript scripts/02_run_peak_picking.R <project_name>
#
# Assumes data/<project>/metadata/sample_sheet.xlsx has already been reviewed
# (see 01_generate_sample_sheet.R).

source("R/peak_picking.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/02_run_peak_picking.R <project_name>", call. = FALSE)
}
project <- args[1]

sheet_path <- file.path("data", project, "metadata", "sample_sheet.xlsx")
output_root <- file.path("data", project, "output")

sample_sheet <- as.data.frame(readxl::read_excel(sheet_path))

groups <- split(sample_sheet, paste(sample_sheet$column, sample_sheet$polarity, sep = "_"))

for (group_name in names(groups)) {
  group_sheet <- groups[[group_name]]
  message(sprintf(
    "\n=== Processing group: %s (%d files) ===", group_name, nrow(group_sheet)
  ))

  out_dir <- file.path(output_root, group_name)
  centwave_param <- run_ipo_optimization(group_sheet, out_dir)
  run_peak_picking_group(group_sheet, centwave_param, out_dir)
}

message("\nAll groups processed.")
