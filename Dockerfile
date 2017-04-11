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

# install zippy
#ADD . /zippy
#RUN cd zippy && make webservice

EXPOSE 80

CMD /bin/bash
