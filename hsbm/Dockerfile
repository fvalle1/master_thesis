FROM ubuntu:latest
LABEL maintainer='Filippo Valle <filippo.valle@edu.unito.it>'

RUN apt-get update
RUN apt-get install sudo apt-utils --yes
RUN apt-get install libboost-all-dev python3 python3-scipy python-numpy \
              libcgal-dev libcairo2-dev python-cairo libsparsehash-dev   \
	      autoconf-archive pkg-config --yes

ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Europe/Rome /etc/localtime
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get install ipython libgtk-3-dev libjs-mathjax pandoc --yes
RUN pip3 install gobject matplotlib pandas jupyter cairocffi


RUN echo "deb http://downloads.skewed.de/apt/ bionic main" >> /etc/apt/sources.list
RUN echo "deb-src http://downloads.skewed.de/apt/ bionic main" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keys.openpgp.org --recv-key 612DEFB798507F25
RUN apt-get update && apt-get install python3-graph-tool --yes


RUN useradd -m -g users filippo

ENV PYTHONIOENCODING=utf8

WORKDIR /home/filippo/
#COPY hSBM.ipynb /home/filippo/
COPY sbmtm.py /home/filippo/
ADD run.py /home/filippo/

RUN chmod a+x /home/filippo/run.py

RUN mkdir /home/filippo/files

USER filippo
ENTRYPOINT ["python", "/home/filippo/run.py"]
CMD ["2"]
