# xcms-pipeline

Untargeted LC-MS pipeline: derives sample metadata from mzML filenames, checks QC
injection quality, then runs IPO2-optimized `xcms` peak picking, retention-time
alignment, correspondence, and gap filling to produce one aligned feature table per
column x polarity method (e.g. RP_POS). Stops at the aligned feature table — blank
filtering, batch correction, and clustering are a separate later stage.

Modeled after [github.com/MetaboComp/xcms_pipeline](https://github.com/MetaboComp/xcms_pipeline).

## Requirements

Docker (with the WSL2 backend on Windows). Nothing else needs to be installed
locally — R, Bioconductor, and all packages are built into the image.

## Setup on a new machine

```
git clone <repo-url> xcms-pipeline
cd xcms-pipeline
docker build -t xcms-pipeline .
```

The build compiles the full Bioconductor stack (`xcms`, `mzR`, `MSnbase`, `IPO2`
from [gitlab.com/CarlBrunius/IPO2](https://gitlab.com/CarlBrunius/IPO2) for
centWave optimization — the same package the reference pipeline itself depends
on — and the legacy Bioconductor `IPO` package for retention-time/correspondence
optimization) — expect it to take a while.

## Usage

Point every command at a folder of `.mzML` files via a bind mount. Replace the path
below with wherever your data actually lives (Windows: `C:\path\to\data`, macOS/Linux:
`/path/to/data`).

Runs (especially peak picking) can take a long time and aren't always watched live —
append `2>&1 | tee /path/to/data/logs/$(date +%Y%m%d_%H%M%S).log` to any command below
to keep a full record on disk as well as on screen (create the `logs` folder first).

**1. Generate the sample sheet**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/generate_sample_sheet.R /data
```
Writes `/path/to/data/metadata/sample_sheet.xlsx`. Open it and review: fill in
`sample_group`, `instrument` (must match an entry in `R/instrument_params.R`, e.g.
`MRT`), and `notes`. Batch/plate/QC-type/injection-order are all derived from the
filenames automatically. Set `include` to `FALSE` for any row you want excluded
from QC checking and peak picking without deleting it from the sheet — defaults
to `TRUE`.

**2. Check QC quality**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/check_qc_quality.R /data
```
Flags likely-faulty QC injections (missed injections, empty vials) via TIC and
aligned feature count, writing `qc_flagged`/`qc_flag_reason` back into the sheet and
an interactive report to `metadata/qc_quality_report.html` (colored by batch, hover
for exact values, faded = flagged). Optional second argument `false` excludes `ltQC`
entirely instead of pooling it with `sQC`:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/check_qc_quality.R /data false
```

**3. Run peak picking**
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data
```
Optional arguments: `[ipo_scope]` (`global` default, or `batch` for a separate IPO
optimization per batch), `[ipo_subset_size]` (default `4` — how many files/batches
IPO2 optimizes against), and `[ipo_fresh]` (default `false`):
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data batch 8
```
If interrupted mid-search (killed container, crash, power loss), each group/batch's
IPO2 optimization checkpoints its best result so far to `ipo_checkpoint.rds` —
re-running resumes from that checkpoint rather than starting the ~100-evaluation
search over. Pass `true` for `[ipo_fresh]` to ignore any cached result or checkpoint
and re-optimize every group/batch from scratch instead:
```
docker run --rm -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data global 4 true
```

By default, parallel steps use all detected cores minus 2. Override with the
`XCMS_PIPELINE_CORES` env var (applies to peak picking and the QC quality check
alike, since both run parallel steps):
```
docker run --rm -e XCMS_PIPELINE_CORES=8 -v /path/to/data:/data xcms-pipeline Rscript scripts/run_peak_picking.R /data
```

## Output

Everything lands back on the host under `/path/to/data/` (bind-mounted, never baked
into the image):

```
metadata/
  sample_sheet.xlsx           reviewed sample metadata
  qc_quality_report.html      interactive QC charts
output/<column>_<polarity>/
  ipo_params.rds               optimized centWave parameters (global scope
                                only -- under <batch>/ipo_params.rds per
                                batch instead in batch scope)
  ipo_checkpoint.rds            mid-search checkpoint (deleted on success;
                                only present if interrupted partway through)
  ipo_history.rds/.csv          IPO2 optimization result summary (final
                                solution, score, nloptr status/iterations --
                                not a full per-iteration trace)
  batch_centwave_params.csv     batch scope only -- every batch's actual
                                centWave params side by side, one row each
  retgroup_params.rds           optimized obiwarp + correspondence
                                bandwidth/bin-size (legacy
                                IPO::optimizeRetGroup())
  retgroup_history.rds          full optimizeRetGroup() result object
  peaks/xdata.rds               full aligned XCMSnExp object
  peaks/peak_table.csv          flat per-peak table (includes gap-filled
                                peaks), enriched with sample_name, is_filled,
                                ms_level, and feature (which aligned feature
                                this raw peak belongs to, if any)
  feature_table.csv             aligned feature table: one row per feature,
                                one column per sample
```

`retgroup_params.rds` is cached and participates in `[ipo_fresh]` the same way
as `ipo_params.rds` — delete it, or pass `true` for `[ipo_fresh]`, to force that
group's retention-time/correspondence search to redo. It doesn't checkpoint
mid-search yet, unlike the centWave search — `IPO::optimizeRetGroup()` is a
first integration, untested outside this project so far.

## Notes

- Every script is re-runnable: the sample sheet and QC steps back up the sheet
  before overwriting it, and peak picking caches IPO2 results (and checkpoints
  mid-search progress) so a re-run doesn't re-optimize from scratch.
- `R/instrument_params.R` holds per-instrument/column/polarity CentWave starting
  values and IPO2 search bounds — edit it to add or tune instruments.
- Rebuild the Docker image after pulling any code change; a container doesn't
  pick up new source automatically.
