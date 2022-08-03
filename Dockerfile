FROM ubuntu:latest

RUN set -xe && apt-get update && apt-get install -y python3 python3-pip 
RUN pip install --upgrade pip

RUN apt-get update
RUN apt-get install -y build-essential 
RUN apt-get install -y wget
RUN apt-get clean
RUN apt-get install -y apt-utils
RUN apt-get install python3-setuptools
RUN pip install --upgrade setuptools
RUN rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
RUN pip install pandas
RUN pip install numpy
RUN pip install psycopg2-binary

RUN apt-get update \
    && apt-get install -y postgresql-server-dev-all gcc python3-dev musl-dev
RUN apt-get install libpq-dev
RUN apt-get install -y git

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH $CONDA_DIR/bin:$PATH

WORKDIR /app
COPY . /app

RUN chmod +x environment_docker.sh
RUN ./environment_docker.sh

CMD ["python", "main.py", "--docker"]
