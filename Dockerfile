FROM rocker/r-ver:4.4.0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Install system dependencies including those needed for devtools
RUN apt-get update && apt-get install -y \
    git \
    curl \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    make \
    gcc \
    g++ \
    pandoc \
    zlib1g-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libudunits2-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# First, install basic packages with detailed output
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); \
    install.packages("remotes"); \
    install.packages("BiocManager"); \
    install.packages("usethis"); \
    print(installed.packages()[,"Package"])'

# Install devtools with verbose output to see any errors
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); \
    install.packages("devtools", dependencies=TRUE, verbose=TRUE); \
    print(installed.packages()[,"Package"])'

# Install additional required packages for ichorCNA
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); \
    install.packages(c("optparse", "getopt", "tidyverse"))'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(version = "3.20", ask = FALSE); \
    BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "biomaRt", "HMMcopy"), ask = FALSE)'

# Verify that devtools is installed and functioning
RUN Rscript -e 'library(devtools); sessionInfo()'

# Clone and install ichorCNA using remotes instead of devtools
RUN git clone https://github.com/broadinstitute/ichorCNA.git /opt/ichorCNA && \
    Rscript -e 'remotes::install_local("/opt/ichorCNA", dependencies = TRUE)'

# Install HMMcopy utils
RUN apt-get update && apt-get install -y cmake && \
    git clone https://github.com/shahcompbio/hmmcopy_utils.git /opt/hmmcopy_utils && \
    cd /opt/hmmcopy_utils && \
    cmake . && \
    make && \
    find bin -type f -executable -exec cp {} /usr/local/bin/ \; && \
    find util -type f -executable -exec cp {} /usr/local/bin/ \;

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
