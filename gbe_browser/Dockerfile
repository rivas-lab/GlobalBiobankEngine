FROM conda-ubuntu
#FROM continuumio/anaconda
#FROM python:2.7
#RUN mkdir /opt/biobankengine && apt-key adv --keyserver keys.gnupg.net --recv-key 6212B7B7931C4BB16280BA1306F90DE5381BA480 && echo "deb http://cran.cnr.berkeley.edu/bin/linux/debian jessie-cran3/" >> /etc/apt/sources.list && conda remove --force readline -y && apt-get update && apt-get install r-base r-base-dev -y
RUN mkdir /opt/biobankengine && /root/anaconda2/bin/conda remove --force readline -y && apt-get update && apt-get install git-core python-pip r-base r-base-dev -y && pip install pip --upgrade
ADD requirements.txt /opt/biobankengine 
ADD requirements-scidb.txt /opt/biobankengine
WORKDIR /opt/biobankengine
RUN pip install -r requirements.txt
RUN pip install -r requirements-scidb.txt
RUN pip install dash==0.21.0  # The core dash backend
RUN pip install dash-renderer==0.11.3  # The dash front-end
RUN pip install dash-html-components==0.8.0  # HTML components
RUN pip install dash-core-components==0.18.1  # Supercharged components
EXPOSE 8000
ENV SCIDB_URL='http://172.18.0.1:8080'
