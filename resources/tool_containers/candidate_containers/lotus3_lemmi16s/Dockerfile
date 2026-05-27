FROM ubuntu:22.04

####
#
# 1) install the common required stuff here
#

ENV PACKAGES="make gcc wget libc6-dev zlib1g-dev ca-certificates xz-utils"
RUN apt-get update --fix-missing -y && apt-get install -y --no-install-recommends ${PACKAGES}
RUN apt-get install --fix-missing python3 git build-essential -y
#RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  List::MoreUtils!)"
#RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  Parallel::ForkManager!)"

RUN apt-get update && apt-get install bc -y
RUN apt-get update && apt-get install -y wget bzip2 git  python3-pip

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#

# Install conda, then install LotuS2
#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
#    chmod +x /tmp/miniconda.sh && \
#    /tmp/miniconda.sh -b -p /miniconda && \
#    rm /tmp/miniconda.sh
#ENV PATH="/miniconda/bin:$PATH"

#ENV CONDA_ACCEPT_CHANNELS=true

#RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
#    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

#RUN conda create -y -n lotus2_env -c bioconda -c conda-forge lotus3 && \
#   conda clean -afy
#RUN conda create -c conda-forge -c bioconda --strict-channel-priority -n LotuS3 LotuS3
#RUN conda clean -afy
#RUN rm -rf /miniconda/envs/LotuS3/share/lotus3-3.03-1/DB/*
WORKDIR /app
RUN git clone https://github.com/hildebra/LotuS3.git && \
    cd LotuS3 && \
    perl helpers/autoInstall.pl && \
    cp lotus3 /usr/bin

# Activate environment by default
SHELL ["bash", "-c"]
#RUN echo "conda activate LotuS3" >> ~/.bashrc
#ENV CONDA_DEFAULT_ENV=LotuS3

RUN pip3 install pandas
RUN pip install biom-format


####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#
WORKDIR /data
ADD scripts/sdm_PacBio_LSSU.txt .
ADD scripts/*.sh /usr/bin/
ADD scripts/*.py /usr/bin/
# note: a LEMMI container needs the commands wget, ps, and bc
