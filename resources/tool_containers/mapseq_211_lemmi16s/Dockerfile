FROM ubuntu:20.04 AS build
ARG DEBIAN_FRONTEND=noninteractive
####
#
# 1) install the common required stuff here
#

ENV PACKAGES make gcc wget libc6-dev zlib1g-dev ca-certificates xz-utils
RUN apt-get update --fix-missing -y && apt-get install -y --no-install-recommends ${PACKAGES}
RUN apt-get install --fix-missing python3 git build-essential -y

RUN apt-get update && apt-get install bc -y

RUN apt-get update && apt-get install python3-pip -y
RUN cp /usr/bin/python3 /usr/bin/python

RUN pip3 install --no-cache-dir \
    pandas 

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#

ENV MAPSEQ_DIR /MAPseq
WORKDIR $MAPSEQ_DIR
RUN wget https://github.com/jfmrod/MAPseq/releases/download/v2.1.1/mapseq-2.1.1-linux.tar.gz
RUN tar xvzf mapseq-2.1.1-linux.tar.gz
ENV PATH ${PATH}:$MAPSEQ_DIR/mapseq-2.1.1-linux/

ENV SW_NAME=Trimmomatic
ENV VERSION=0.38
ENV DEST=/usr/src/$SW_NAME
ENV ADAPTERPATH=$DEST/$VERSION/$SW_NAME-$VERSION/adapters

RUN apt-get -y update --fix-missing && \
	apt-get install -y build-essential zip wget default-jre 
WORKDIR $DEST/$VERSION
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip -O $SW_NAME-$VERSION.zip && \
	unzip $SW_NAME-$VERSION.zip && \
	echo "alias trimmomatic='java -jar $DEST/$VERSION/$SW_NAME-$VERSION/trimmomatic-0.38.jar'" >> ~/.bash_aliases && \
        echo "source ~/.bash_aliases" >> ~/.bashrc 


RUN apt install seqprep -y
####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#

ADD /scripts/*.sh /usr/bin/
ADD /scripts/*.py /usr/bin/
RUN cd ~/ 
# note: a LEMMI16s container needs the commands wget, ps, and bc
