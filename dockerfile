FROM ubuntu:latest

WORKDIR /app
ADD . /app

RUN set -xe \
    && apt-get -y update \
    && apt-get -y install python3-pip

RUN pip install --upgrade pip
