# xcms-pipeline, containerized to freeze the exact R/Bioconductor
# environment recorded in renv.lock (R 4.6.0, Bioconductor 3.23 devel,
# including IPO2 from gitlab.com/CarlBrunius/IPO2 for centWave optimization,
# and the legacy Bioconductor IPO package for retention-time/correspondence
# optimization) — avoids re-fighting the version-specific issues solved
# during initial setup (IPO's centWave optimizer hitting a BiocParallel
# bpstopOnError incompatibility, a separate BiocParallel MulticoreParam
# worker bug, mzR/Rcpp ABI mismatches).
#
# Build:
#   docker build -t xcms-pipeline .
#
# Run (mount your data folder, then invoke whichever script you need):
#   docker run --rm -v /path/to/mzml/folder:/data xcms-pipeline \
#     Rscript scripts/generate_sample_sheet.R /data
#   docker run --rm -v /path/to/mzml/folder:/data xcms-pipeline \
#     Rscript scripts/check_qc_quality.R /data
#   docker run --rm -v /path/to/mzml/folder:/data xcms-pipeline \
#     Rscript scripts/run_peak_picking.R /data

FROM bioconductor/bioconductor_docker:devel

# netCDF headers for mzR, in case not already present in the base image
RUN apt-get update && apt-get install -y --no-install-recommends \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy renv scaffolding first (before application code) so package
# installation is cached by Docker across rebuilds that only touch R/ or
# scripts/
COPY .Rprofile renv.lock ./
COPY renv/activate.R renv/settings.json ./renv/
COPY docker/install_packages.R ./docker/install_packages.R

# Sandboxing protects a shared host library from a project; redundant inside
# an already-isolated container, and avoids permission-write noise.
ENV RENV_CONFIG_SANDBOX_ENABLED=FALSE

# The base image's own Rprofile pre-loads BiocManager before renv activates
# this project, which renv's namespace check otherwise flags on every
# status()/snapshot() call. Harmless and structural to the base image (not
# fixable from here), so silenced rather than left as recurring noise.
ENV RENV_CONFIG_NAMESPACES_CHECK=FALSE
RUN Rscript -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
RUN Rscript docker/install_packages.R

COPY R/ ./R/
COPY scripts/ ./scripts/
