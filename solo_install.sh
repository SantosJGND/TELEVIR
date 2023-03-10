#/bin/bash

#CONDA_DIR=/opt/conda
PYTHON_PATH=.venv

sudo apt-get update
sudo apt install -y python3.10-venv python3-setuptools
sudo apt-get -y install libpq-dev default-jre install zlib1g-dev make g++
sudo apt install -y tabix ncbi-blast+ samtools ncbi-entrez-direct bioperl
sudo apt-get install -y build-essential wget curl rsync apt-utils unzip



#
#wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#/bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
#rm ~/miniconda.sh
#echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
#echo "conda activate base" >> ~/.bashrc
#echo "export PATH=$CONDA_DIR/bin:$PATH" >> ~/.bashrc
#source ~/.bashrc


python -m venv $PYTHON_PATH

$PYTHON_PATH/bin/pip install numpy==1.23.0 psycopg2-binary pandas==1.4.3 sqlalchemy==1.4.41 python-decouple
$PYTHON_PATH/bin/pip install --default-timeout=100 future

#echo "---> Create televir dbs in /televir/mngs_benchmark/ ..."
#
#if [ ! -d "/televir/mngs_benchmark" ]; then
#    mkdir -p /televir/mngs_benchmark
#fi
#
#cd insaflu_web/TELEVIR

$PYTHON_PATH/bin/python main.py --envs --setup_conda --partial
$PYTHON_PATH/bin/python main.py --envs --partial

$PYTHON_PATH/bin/python main.py --seqdl --partial
$PYTHON_PATH/bin/python main.py --soft --partial
$PYTHON_PATH/bin/python main.py --deploy --partial


#cp utility.db /televir/mngs_benchmark/utility.db
#cp install_scripts/config.py /televir/mngs_benchmark/config.py

#echo "---> Finshed creating televir dbs in /televir/mngs_benchmark/ ... done"