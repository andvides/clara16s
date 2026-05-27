FROM ubuntu:latest

####
#
# 1) install the common required stuff here
#
ENV PACKAGES make gcc wget libc6-dev zlib1g-dev ca-certificates xz-utils
RUN apt-get update --fix-missing -y && apt-get install -y --no-install-recommends ${PACKAGES}
RUN apt-get install --fix-missing python3 git build-essential -y
#RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  List::MoreUtils!)"
#RUN yes | perl -MCPAN -e "CPAN::Shell->notest(qw!install  Parallel::ForkManager!)"

RUN apt-get update && apt-get install bc -y
RUN apt-get install unzip
####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#
#RUN apt install -y snakemake

RUN yes | apt-get install libboost-all-dev
RUN yes | apt-get install libreadline-dev
RUN yes | apt-get install libhdf5-dev
RUN apt-get update
RUN yes | apt-get install libgsl-dev

ENV MOTHUR /app
WORKDIR $MOTHUR
RUN wget https://github.com/mothur/mothur/releases/download/v1.48.0/Mothur.Ubuntu_20.zip
RUN unzip Mothur.Ubuntu_20.zip
RUN ls
RUN cp mothur/* /usr/bin/

RUN wget https://github.com/dstreett/FLASH2/archive/refs/tags/2.2.00.zip
RUN unzip 2.2.00.zip 
RUN make -C FLASH2-2.2.00/
RUN cp FLASH2-2.2.00/flash2 /usr/bin/

####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#

RUN apt-get update && apt-get install python3-pip -y
RUN pip3 install pandas
RUN cp /usr/bin/python3 /usr/bin/python

ADD /scripts/*.sh /usr/bin/
ADD /scripts/*.py /usr/bin/

