# Development

Notes for working on this repo — not needed just to run the pipeline (see
README.md for that).

## Tests

Unit tests cover the pure/deterministic logic (filename parsing, sample sheet
handling, instrument param lookup, IPO subset selection) — not the parts that
need real mzML files or a live Bioconductor run. No Docker required, just
`testthat`:
```
install.packages("testthat")   # one-time
Rscript tests/testthat.R
```
Also runs inside the image if you'd rather not install anything locally:
```
docker run --rm xcms-pipeline Rscript tests/testthat.R
```

Run them before any commit that touches `R/` — they're fast, and they exist
specifically to catch a logic slip before it costs a real (slow) Docker run.
A green test suite says nothing about whether the actual xcms/IPO2/IPO
integration works, though — only a real pipeline run confirms that.
