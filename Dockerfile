FROM ubuntu:18.04

RUN apt-get update && apt-get install python3 python3-pip -y

RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
    python pip install --upgrade pip
    python pip install pandas
    python pip install numpy


# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

COPY . /usr/local/bin/
WORKDIR /usr/local/bin/

CMD ["python", "main.py", "--docker"]

