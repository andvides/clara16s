FROM ubuntu:22.04

LABEL org.opencontainers.image.title="CLARA16s Kraken2-Bracken"
LABEL org.opencontainers.image.description="Kraken 2 and Bracken tools for CLARA16s/LEMMI16s"
LABEL org.opencontainers.image.source="https://github.com/DerrickWood/kraken2"

ARG DEBIAN_FRONTEND=noninteractive
ARG KRAKEN2_VERSION=2.17.2
ARG BRACKEN_VERSION=3.1

ENV KRAKEN2_DIR=/opt/kraken2 \
    BRACKEN_DIR=/opt/bracken \
    DB_DIR=/databases \
    PATH="/opt/venv/bin:/opt/kraken2:/opt/bracken:${PATH}"

# System dependencies for Kraken 2, Bracken, R, Python and existing pipeline scripts.
RUN apt-get update && apt-get install -y --no-install-recommends \
    bc \
    build-essential \
    ca-certificates \
    git \
    gpg \
    liblist-moreutils-perl \
    libparallel-forkmanager-perl \
    ncbi-blast+ \
    perl \
    python3 \
    python3-pip \
    python3-venv \
    rsync \
    wget \
    xz-utils \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R and Biostrings, retaining the dependencies used by the existing pipeline.
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
      > /etc/apt/sources.list.d/cran.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base \
    && R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('Biostrings', ask=FALSE); install.packages('optparse', repos='https://cloud.r-project.org')" \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies in a virtual environment to avoid modifying system Python.
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir pandas ete3 legacy-cgi tqdm \
    && ln -sf /usr/bin/python3 /usr/local/bin/python

# Install Kraken 2 from its official source repository.
RUN git clone --depth 1 --branch "${KRAKEN2_VERSION}" \
    https://github.com/DerrickWood/kraken2.git /tmp/kraken2-src \
    && bash /tmp/kraken2-src/install_kraken2.sh "${KRAKEN2_DIR}" \
    && rm -rf /tmp/kraken2-src

# Install Bracken from its official source repository.
RUN git clone --depth 1 --branch "v${BRACKEN_VERSION}" \
      https://github.com/jenniferlu717/Bracken.git /tmp/bracken-src \
    && cd /tmp/bracken-src && bash install_bracken.sh \
    && mkdir -p "${BRACKEN_DIR}" \
    && cp -a /tmp/bracken-src/. "${BRACKEN_DIR}/" \
    && rm -rf /tmp/bracken-src

# External database mount point. Database downloads/builds are intentionally
# excluded from the image to keep it smaller and independently versionable.
RUN mkdir -p "${DB_DIR}" /work

# Copy CLARA16s/LEMMI16s helper scripts when building from the repository root.
COPY scripts/*.r /usr/bin/
COPY scripts/*.sh /usr/bin/
COPY scripts/*.py /usr/bin/

WORKDIR /work

# Fail the build if the main executables are not available.
RUN kraken2 --version \
    && bracken -h >/dev/null 2>&1 \
    && bracken-build -h >/dev/null 2>&1 \
    && Rscript -e "library(Biostrings); library(optparse)" \
    && python -c "import pandas, ete3, tqdm"

CMD ["bash"]

