FROM ubuntu:latest

RUN set -xe && apt-get update && apt-get install -y python3 python3-pip 
RUN pip install --upgrade pip

RUN apt-get update
RUN apt-get install -y build-essential wget curl rsync apt-utils python3-setuptools

RUN apt-get clean
#RUN apt-get install python3-setuptools
#RUN pip install --upgrade setuptools

RUN rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
RUN pip install pandas
RUN pip install numpy
RUN pip install psycopg2-binary
RUN apt-get update
RUN apt-get update \
    && apt-get install -y postgresql-server-dev-all gcc python3-dev musl-dev
RUN apt-get install libpq-dev
RUN apt-get install -y git ncbi-entrez-direct tabix samtools bioperl 

RUN apt-get install default-jre python3-pip python3-dev libpq-dev postgresql postgresql-contrib

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH $CONDA_DIR/bin:$PATH

## install kallisto version
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
RUN tar -xzf kallisto_linux-v0.43.1.tar.gz
ENV PATH $PATH:kallisto_linux-v0.43.1/kallisto

RUN rm kallisto_linux-v0.43.1.tar.gz
RUN rm -rf kallisto_linux-v0.43.1
RUN rm ._kallisto_linux-v0.43.1

##
WORKDIR /app
COPY . /app

RUN chmod +x environment_docker.sh
RUN ./environment_docker.sh 

## django requirements
RUN python -m pip install -r requirements_django.txt

RUN python main.py --docker --envs --partial
RUN python main.py --docker --seqdl --partial
RUN python main.py --docker --soft --partial