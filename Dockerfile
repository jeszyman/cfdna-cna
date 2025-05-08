FROM ubuntu:22.04

ENV PATH=/opt/conda/bin:$PATH
ENV R_LIBS_USER=/opt/conda/envs/ichor/lib/R/library

RUN apt-get update && apt-get install -y \
    curl wget git build-essential libcurl4-openssl-dev \
    libssl-dev libxml2-dev libgit2-dev libfontconfig1-dev \
    pandoc

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/mamba.sh && \
    bash /tmp/mamba.sh -b -p /opt/conda && \
    rm /tmp/mamba.sh

COPY ./config/ichor-conda-env.yaml /opt/ichor-conda-env.yaml
RUN /opt/conda/bin/mamba env create -f /opt/ichor-conda-env.yaml

RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate ichor && \
    Rscript -e 'BiocManager::install(version = \"3.20\")' && \
    git clone https://github.com/broadinstitute/ichorCNA.git /opt/ichorCNA && \
    Rscript -e 'devtools::install_local(\"/opt/ichorCNA\", dependencies = TRUE)'"

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
