# xcms-pipeline

Untargeted LC-MS pipeline: derives sample metadata from mzML filenames, checks QC
injection quality, then runs IPO-optimized `xcms` peak picking, alignment,
correspondence, and gap filling into one aligned feature table per column x
polarity method (e.g. RP_POS). Stops there — blank filtering, batch correction,
and clustering are a later stage.

Modeled after [github.com/MetaboComp/xcms_pipeline](https://github.com/MetaboComp/xcms_pipeline).

## Requirements

Docker (WSL2 backend on Windows). Nothing else needed locally — R, Bioconductor,
and all packages are built into the image.

## Setup

```
git clone <repo-url> xcms-pipeline
cd xcms-pipeline
docker build -t xcms-pipeline .
```

Builds the full Bioconductor stack (`xcms`, `mzR`, `MSnbase`, `IPO2`, `IPO`) —
takes a while.

## Usage

Every command takes a folder of `.mzML` files via a bind mount. Replace
`/path/to/data` below with wherever your data lives (Windows: `C:\path\to\data`).

For long runs, add `2>&1 | tee /path/to/data/logs/$(date +%Y%m%d_%H%M%S).log`
to keep a log on disk (make the `logs` folder first).

**1. Generate the sample sheet**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/generate_sample_sheet.R /data
```
Writes `/path/to/data/metadata/sample_sheet.xlsx`. Review it: fill in
`sample_group`, `instrument` (must match `R/instrument_params.R`, e.g. `MRT`),
and `notes`. Batch/plate/QC-type/injection-order are parsed from filenames
automatically. Set `include` to `FALSE` to exclude a row without deleting it
(default `TRUE`).

A file whose name doesn't match the expected pattern doesn't block the rest
of the sheet — it's still added, with just its filename and `needs_review =
TRUE`/`parse_error` explaining why, `include` defaulted to `FALSE`. Fill in
its `batch`/`column`/`polarity`/`sample_name`/etc. by hand and flip `include`
back to `TRUE` once it's ready.

**2. Check QC quality**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/check_qc_quality.R /data
```
Flags likely-faulty QC injections (missed injections, empty vials) via raw
peak count (primary) and aligned feature count (secondary, more lenient).
Writes `qc_flagged`/`qc_flag_reason` back into the sheet and an interactive
report to `metadata/qc_quality_report.html`. Pass `false` to check `sQC`
alone instead of also checking `ltQC`:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/check_qc_quality.R /data false
```

**3. Run peak picking**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data
```
Optional args: `[ipo_scope]` (`global` default, or `batch` for a separate IPO
optimization per batch), `[ipo_subset_size]` (default `4`), `[ipo_fresh]`
(default `false`), `[retgroup_qc_type]` (default `auto`), `[center_sample]`
(default `qc`):
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data batch 8
```
An interrupted run resumes from its last checkpoint by default. Pass `true` for
`[ipo_fresh]` to ignore any cache/checkpoint and start over:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data global 4 true
```
Each group's retention-time/correspondence optimization normally picks sQC
over ltQC automatically (falling back only if sQC is too thin) -- a
"QC batch coverage" printout at the start of each group shows how many
batches each type actually covers. It searches against up to 5 files per
batch (spread evenly across `injection_order`), not every QC file in the
group -- every search evaluation re-aligns whichever files are passed in
from scratch, so this keeps a large study's search cost from scaling with
its total QC count while still representing every batch. Pass `sQC` or
`ltQC` as `[retgroup_qc_type]` to force one instead of relying on the
automatic pick:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data batch 5 false ltQC
```
The obiwarp alignment step's reference file ("center sample" -- the one every
other file gets warped against) is picked the same way by default (`qc`): a
non-flagged sQC (falling back to ltQC) closest to the group's median
`injection_order`, logged so it's traceable. Pass `middle` as `[center_sample]`
to use xcms's own default instead -- whichever file sits at the positional
median index, with no regard for quality:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data batch 5 false auto middle
```

By default, parallel steps use all cores minus 2. Override with
`XCMS_PIPELINE_CORES`:
```
docker run --rm -e XCMS_PIPELINE_CORES=8 -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data
```
The final alignment/correspondence/gap-filling step (`align_and_correspond()`)
runs against the *full* group (every sample/blank/QC, not a small subsample),
so it always uses a fixed 8 workers for `groupChromPeaks()`/`fillChromPeaks()`
regardless of that setting -- independent of how many cores the rest of the
pipeline uses. `adjustRtime()` within that same step always runs
single-threaded, not parallel at all -- a memory-safety fix, see
`align_and_correspond()`'s comments.

## Output

Lands back on the host under `/path/to/data/` (never baked into the image):

```
metadata/
  sample_sheet.xlsx           reviewed sample metadata
  qc_quality_report.html      interactive QC charts
output/
  sample_sheet_snapshot_<timestamp>.xlsx   the sheet exactly as used for
                                            this run -- edit and rerun the
                                            live sheet freely, this stays
                                            put so old output stays traceable
                                            (also copied into each group's
                                            own folder below, so a group's
                                            output is traceable on its own
                                            even if copied elsewhere)
output/<column>_<polarity>/
  sample_sheet_snapshot_<timestamp>.xlsx   same file, copied here too
  ipo_params.rds               optimized centWave params (per-batch under
                                <batch>/ in batch scope)
  ipo_checkpoint.rds            mid-search checkpoint, only present if
                                interrupted
  ipo_history.rds/.csv          optimizer result summary
  batch_centwave_params.csv     batch scope only: every batch's params,
                                side by side
  picked_peaks.rds              cached findChromPeaks() result (per-batch
                                under <batch>/ in batch scope) -- resume
                                skips re-picking an already-done batch/group;
                                auto-invalidated if files/params changed
  retgroup_params.rds           optimized alignment/correspondence params
  retgroup_history.rds          full optimizeRetGroup() result
  peaks/xdata.rds               full aligned XCMSnExp object
  peaks/peak_table.csv          per-peak table, with sample_name, is_filled,
                                ms_level, feature columns added
  feature_table.csv             aligned feature table: one row per feature,
                                one column per sample
```

`ipo_fresh` also governs `retgroup_params.rds`.

## Notes

- Every script is re-runnable: sample sheet/QC steps back up the sheet before
  overwriting it; peak picking caches IPO results, checkpoints mid-search, and
  caches each batch/group's picked peaks so an interrupted run doesn't
  re-`findChromPeaks()` from scratch.
- `R/instrument_params.R` holds per-instrument/column/polarity CentWave
  settings and search bounds — edit it to add or tune instruments.
- Rebuild the image after any code change; a container won't pick up new
  source on its own.

See `DEVELOPMENT.md` for running tests.
