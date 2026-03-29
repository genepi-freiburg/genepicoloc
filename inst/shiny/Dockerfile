# genepicoloc Colocalization Viewer
#
# Build:
#   docker build -t genepicoloc-shiny .
#
# Run with data mounted:
#   docker run -d -p 3838:3838 \
#     -v /path/to/data:/app/data:ro \
#     --name genepicoloc-shiny genepicoloc-shiny
#
# Data directory should contain:
#   data/
#   ├── <trait>_annot_filt.RDS         # Coloc results (flat layout)
#   ├── gene_annotation_hg38_full.tsv  # Gene annotations (optional, 3.9MB)
#   └── regional_plots/                # Regional plot data (optional)
#       └── <trait>/
#           └── regions/
#               └── <region>/
#                   ├── sumstats_1.RDS
#                   └── traits/
#                       └── <study>__<trait>.RDS

FROM rocker/shiny:latest

# Install system dependencies for plotly
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      libssl-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Install R packages (ordered by change frequency for layer caching)
RUN R -e "install.packages(c('data.table', 'visNetwork', 'plotly'), repos='https://cloud.r-project.org/')"

# Create app directory
WORKDIR /app

# Copy app files
COPY R/ /app/R/
COPY www/ /app/www/
COPY app.R /app/

# Data is mounted at runtime via -v flag
# Create empty data dir as fallback
RUN mkdir -p /app/data

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
