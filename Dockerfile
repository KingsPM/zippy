FROM debian:jessie
MAINTAINER dbrawand@nhs.net

RUN apt-get update && apt-get -y upgrade

# prepare genome
RUN apt-get -y install less make wget
ADD Makefile /zippy/Makefile
RUN cd zippy && make genome-download
ADD package-requirements.txt /zippy/package-requirements.txt
RUN cd zippy && make install
RUN cd zippy && make genome-index

# get annotation
RUN cd zippy && make variation-download
RUN cd zippy && make refgene-download

# add some convenience utils
RUN apt-get -y install curl less vim

# install zippy
ADD . /zippy
RUN cd /zippy && make webservice

EXPOSE 80

RUN echo "ServerName localhost" >> /etc/apache2/apache2.conf
CMD /bin/bash /zippy/zippyd.sh

