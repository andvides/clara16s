# ============================================================
# CLARA16s - Emu
# Species-level taxonomic abundance for full-length 16S reads
# ============================================================

FROM ubuntu:22.04

LABEL maintainer="CLARA16s"
LABEL description="Emu for species-level taxonomic profiling of 16S reads"
LABEL org.opencontainers.image.title="Emu - CLARA16s"

ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------
# 1. System dependencies
# ------------------------------------------------------------

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    git \
    bc \
    procps \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    r-base \
    r-cran-optparse \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# 2. Install Miniforge (AMD64 / ARM64)
# ------------------------------------------------------------

ARG TARGETARCH

RUN case "${TARGETARCH}" in \
        amd64) CONDA_ARCH="x86_64" ;; \
        arm64) CONDA_ARCH="aarch64" ;; \
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;; \
    esac && \
    wget -q \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${CONDA_ARCH}.sh" \
        -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:${PATH}"

# ------------------------------------------------------------
# 3. Install Emu
# ------------------------------------------------------------

WORKDIR /opt

ARG EMU_REF=master

RUN git clone https://github.com/treangenlab/emu.git && \
    cd /opt/emu && \
    git checkout "${EMU_REF}"

# Create an isolated Conda environment
RUN conda config --set channel_priority strict && \
    conda env create \
        -p /opt/conda/envs/emu \
        -f /opt/emu/environment.yml && \
    conda clean -afy

# Activate Emu environment without requiring conda activate
ENV PATH="/opt/conda/envs/emu/bin:/opt/emu:${PATH}"

# ------------------------------------------------------------
# 4. Install R dependencies
# ------------------------------------------------------------

RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(c('Biostrings','dada2','phyloseq'), ask=FALSE, update=FALSE)"

# ------------------------------------------------------------
# 5. Runtime directories and environment
# ------------------------------------------------------------

ENV EMU_DATABASE_DIR=/db

RUN mkdir -p /db /data /results /.config /.cache && \
    chmod -R a+rw /.config /.cache

WORKDIR /app

# ------------------------------------------------------------
# 6. CLARA16s auxiliary scripts
# ------------------------------------------------------------

COPY scripts/ /usr/local/bin/

RUN chmod -R a+rX /usr/local/bin && \
    find /usr/local/bin -type f -exec chmod a+x {} \;

# ------------------------------------------------------------
# 7. Default command
# ------------------------------------------------------------

WORKDIR /data

CMD ["emu"]
