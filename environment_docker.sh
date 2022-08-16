#!/bin/bash

apt-get update
apt-get install default-jre
apt-get install -y bioperl
apt-get install python3-pip python3-dev libpq-dev postgresql postgresql-contrib tabix

python3 -m venv .venv

source .venv/bin/activate

python -m pip install --upgrade pip

## software requirements

python -m pip install -r requirements.txt

## install kallisto version
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
tar -xzf kallisto_linux-v0.43.1.tar.gz
mv kallisto_linux-v0.43.1/kallisto .venv/bin

rm kallisto_linux-v0.43.1.tar.gz
rm -rf kallisto_linux-v0.43.1
rm ._kallisto_linux-v0.43.1
