FROM ubuntu:22.04 AS build
ENV DEBIAN_FRONTEND=noninteractive

####
#
# 1) install the common required stuff here
#

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    git \
    ca-certificates \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    r-base \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade pip setuptools wheel

RUN apt-get update && apt-get install bc -y
RUN apt-get update && apt-get install wget -y

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#

WORKDIR /opt

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    chmod +x /tmp/miniconda.sh && \
    /tmp/miniconda.sh -b -p /miniconda && \
    rm /tmp/miniconda.sh
ENV PATH="/miniconda/bin:$PATH"

RUN git clone https://github.com/treangenlab/emu.git

ENV CONDA_ACCEPT_CHANNELS=true
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create --file emu/environment.yml && \
    conda clean --all -y

RUN conda install conda-forge::r-optparse
RUN mkdir -p /.config /.cache && chmod -R a+rw /.config /.cache

ENV PATH="/opt/emu:${PATH}"

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('Biostrings','dada2','phyloseq'), ask=FALSE, update=FALSE)"


####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#

WORKDIR /app
#ADD scripts/*.sh /usr/bin/
ADD scripts/*.* /usr/bin/


# note: a LEMMI container needs the commands wget, ps, and bc
