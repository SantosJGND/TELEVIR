FROM ubuntu:22.04

ARG APP_USER=televir_user
RUN useradd -ms /bin/bash ${APP_USER}
ENV DEBIAN_FRONTEND=noninteractive
ENV TERM=xterm
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes

RUN set -x \
    && groupadd -r --gid=990 slurm \
    && useradd -r -g slurm --uid=990 slurm

RUN apt-get update && apt-get upgrade -y

RUN apt-get install -y postgresql postgresql-contrib libpq-dev

RUN apt-get install -y python3.11 python3.11-venv python3.11-dev python3-venv python3-pip python3.10-venv
RUN python3 -m pip install --upgrade pip
RUN python3 -m venv /opt/venv

RUN /opt/venv/bin/pip install wheel
RUN /opt/venv/bin/pip install setuptools numpy==1.24.3 pandas==2.0.1 
RUN /opt/venv/bin/pip install --upgrade pip

RUN /opt/venv/bin/pip install psycopg2 sqlalchemy>=2.0 python-decouple==3.8 danio PyYAML xopen==1.7.0 fastq_filter==0.3.0
RUN /opt/venv/bin/pip install --default-timeout=100 future

RUN apt-get update
RUN apt-get install -y build-essential wget curl rsync apt-utils python3-setuptools python3-tk

RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN apt-get update \
    && apt-get install -y gcc musl-dev 

RUN apt-get install -y git ncbi-entrez-direct tabix samtools bioperl 

RUN apt-get -y install default-jre 

ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH 
RUN conda install --name base conda-anaconda-tos
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 

RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
RUN tar -xzf kallisto_linux-v0.43.1.tar.gz
ENV PATH=$PATH:kallisto_linux-v0.43.1

RUN apt-get install zlib1g-dev make g++

WORKDIR /opt/televir

ADD https://api.github.com/repos/SantosJGND/TELEVIR/git/refs/heads/install-stdl version.json
RUN mkdir -p /opt/televir-repo && \
    wget --quiet https://github.com/SantosJGND/TELEVIR/archive/install-stdl.zip -O televir.zip && \
    unzip -q televir.zip && \
    cp -r TELEVIR-install-stdl/. /opt/televir-repo/ && \    
    rm -rf TELEVIR-install-stdl && \
    rm televir.zip && \
    chown -R ${APP_USER}:slurm /opt/televir-repo

ARG REQUEST_SEQ_FILE=""
RUN if [ -n "$REQUEST_SEQ_FILE" ] && [ -f "$REQUEST_SEQ_FILE" ]; then \
    cp $REQUEST_SEQ_FILE /opt/request_sequences.fa.gz; \
    fi

COPY televir.env /opt/televir/televir.env

WORKDIR /
COPY entrypoint.sh entrypoint_original.sh

RUN sed "s/APP_USER/${APP_USER}/g" entrypoint_original.sh > entrypoint.sh && rm entrypoint_original.sh
RUN chmod a+x entrypoint.sh

EXPOSE 8080

ENTRYPOINT [ "/entrypoint.sh" ]
