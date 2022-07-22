#!/bin/bash

sudo apt install default-jre
sudo apt-get install -y bioperl

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

## django requirements
python -m pip install -r requirements_django.txt