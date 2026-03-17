#/bin/bash

#CONDA_DIR=/opt/conda
PYTHON_PATH=.venv

sudo apt-get update
sudo apt install -y python3.10-venv python3-setuptools
sudo apt-get -y install libpq-dev default-jre install zlib1g-dev make g++
sudo apt install -y tabix ncbi-blast+ samtools ncbi-entrez-direct bioperl
sudo apt-get install -y build-essential wget curl rsync apt-utils unzip


python -m venv $PYTHON_PATH

$PYTHON_PATH/bin/pip install numpy==1.23.0 psycopg2-binary pandas==1.4.3 sqlalchemy==1.4.41 python-decouple
$PYTHON_PATH/bin/pip install --default-timeout=100 future

$PYTHON_PATH/bin/python main.py --envs --setup_conda --partial
$PYTHON_PATH/bin/python main.py --envs --partial

$PYTHON_PATH/bin/python main.py --seqdl --partial
$PYTHON_PATH/bin/python main.py --soft --partial
$PYTHON_PATH/bin/python main.py --deploy --partial

