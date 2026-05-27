FROM ubuntu:22.04 AS build
ENV DEBIAN_FRONTEND=noninteractive

####
#
# 1) install the common required stuff here
#

ENV PACKAGES="\
    make \
    gcc \
    g++ \
    wget \
    libc6-dev \
    zlib1g-dev \
    ca-certificates \
    xz-utils \
    rsync \
    perl \
    git \
    ncbi-blast+ \
"

RUN apt-get update --fix-missing -y && apt-get install -y --no-install-recommends ${PACKAGES}
RUN apt-get install --fix-missing python3 git build-essential -y
RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  List::MoreUtils!)"
RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  Parallel::ForkManager!)"

RUN apt-get update && apt-get install bc -y

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#

# Add CRAN repository for R (Ubuntu 22.04)
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
         > /etc/apt/sources.list.d/cran.list

RUN apt-get update && apt-get install -y r-base && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
          BiocManager::install('Biostrings', ask = FALSE)"
RUN R -e "install.packages('optparse', repos='https://cloud.r-project.org')"

ENV KRAKEN2_DIR=/kraken2

RUN git clone https://github.com/DerrickWood/kraken2.git $KRAKEN2_DIR 

WORKDIR $KRAKEN2_DIR
#RUN git checkout v2.1.2
RUN git checkout v2.1.6
RUN ls
RUN ./install_kraken2.sh $KRAKEN2_DIR
ENV PATH="$KRAKEN2_DIR:$PATH"
#RUN cp /kraken2/src/* /kraken2/scripts  
#ENV PATH="${PATH}:$KRAKEN2_DIR/scripts:$KRAKEN2_DIR:$KRAKEN2_DIR/bin:$KRAKEN2_DIR/src:/usr/bin"
##ENV PATH="/usr/bin:${PATH}"

RUN apt-get install ncbi-blast+ -y


ENV BRACKEN_DIR=/bracken2
RUN git clone https://github.com/jenniferlu717/Bracken.git $BRACKEN_DIR
WORKDIR $BRACKEN_DIR
RUN git checkout v3.1
RUN sh install_bracken.sh

####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#
RUN apt-get update && apt-get install python3-pip -y
RUN pip3 install pandas
RUN cp /usr/bin/python3 /usr/bin/python

#SILVA database
ENV SILVA_DIR=/silvaDB
WORKDIR $SILVA_DIR
#RUN kraken2-build --download-taxonomy --db silva138_kraken
#RUN kraken2-build --db . --special silva


RUN  pip install ete3
RUN  pip install legacy-cgi
RUN  pip install tqdm



WORKDIR /app

ADD scripts/*.r /usr/bin/
ADD scripts/*.sh /usr/bin/
ADD scripts/*.py /usr/bin/
# note: a LEMMI container needs the commands wget, ps, and bc
