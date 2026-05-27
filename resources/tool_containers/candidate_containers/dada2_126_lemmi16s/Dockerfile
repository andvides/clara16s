FROM blekhmanlab/dada2
ARG DEBIAN_FRONTEND=noninteractive
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


RUN R -e 'install.packages("Rcpp", repos="https://rcppcore.github.io/drat")'

ADD /scripts/*.sh /usr/bin/
ADD /scripts/*.py /usr/bin/
ADD /scripts/*.R /usr/bin/
