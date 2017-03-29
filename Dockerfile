FROM debian:jessie
MAINTAINER dbrawand@nhs.net

RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install git-core less make

RUN mkdir ~/.ssh && ssh-keyscan -t rsa github.com >> ~/.ssh/known_hosts
RUN git clone --branch develop http://github.com/viapath/zippy.git

RUN cd zippy && make install
#RUN cd zippy && make resources
#RUN cd zippy && make webservice

EXPOSE 80

CMD /bin/bash
