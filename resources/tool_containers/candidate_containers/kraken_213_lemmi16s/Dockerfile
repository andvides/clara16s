FROM ubuntu:latest

####
#
# 1) install the common required stuff here
#

ENV PACKAGES make gcc wget libc6-dev zlib1g-dev ca-certificates xz-utils
RUN apt-get update --fix-missing -y && apt-get install -y --no-install-recommends ${PACKAGES}
RUN apt-get install --fix-missing python3 git build-essential -y
RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  List::MoreUtils!)"
RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  Parallel::ForkManager!)"

RUN apt-get update && apt-get install bc -y

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#

ENV KRAKEN2_DIR /kraken2

RUN git clone https://github.com/DerrickWood/kraken2.git $KRAKEN2_DIR 

WORKDIR $KRAKEN2_DIR

#RUN git checkout v2.1.2
RUN git checkout v2.1.3 

RUN ./install_kraken2.sh $KRAKEN2_DIR

RUN apt-get install ncbi-blast+ -y

RUN cp /kraken2/src/* /kraken2/scripts  

ENV PATH="${PATH}:$KRAKEN2_DIR/scripts:$KRAKEN2_DIR:$KRAKEN2_DIR/bin:$KRAKEN2_DIR/src:/usr/bin"

#ENV PATH="/usr/bin:${PATH}"

ENV BRACKEN_DIR /bracken2

RUN git clone https://github.com/jenniferlu717/Bracken.git $BRACKEN_DIR

WORKDIR $BRACKEN_DIR

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
ENV SILVA_DIR /silvaDB
WORKDIR $SILVA_DIR
RUN kraken2-build --db . --special silva

WORKDIR /app

ADD scripts/*.sh /usr/bin/
ADD scripts/*.py /usr/bin/
# note: a LEMMI container needs the commands wget, ps, and bc
