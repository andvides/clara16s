FROM quay.io/qiime2/core:2023.2

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
RUN apt-get install unzip

RUN apt install nano less

####
#
# ***2) HERE INSTALL YOUR SOFTWARE.
#
#RUN apt install -y snakemake
RUN apt-get update && apt-get install python3-pip -y

RUN pip install xmltodict
RUN pip install ncbi-datasets-pylib

#ENV TOOLS_PATH /app
#WORKDIR $TOOLS_PATH
RUN pip install git+https://github.com/qiime2/q2-types
RUN pip install git+https://github.com/bokulich-lab/q2-types-genomics.git
RUN pip install git+https://github.com/bokulich-lab/RESCRIPt.git
RUN qiime dev refresh-cache
####
#
# ***3) ADD the additional files, all the sh and python scripts, including the 4 LEMMI_ scripts in /usr/bin so they can be executed 
# don't forget to make the LEMMI_* scripts executable (chmod +x)
#

#RUN apt-get update && apt-get install python3-pip -y
#RUN pip3 install pandas
#RUN cp /usr/bin/python3 /usr/bin/python

#ENV TOOLS_PATH /app

ADD /scripts/*.sh /usr/bin/
ADD /scripts/*.py /usr/bin/
