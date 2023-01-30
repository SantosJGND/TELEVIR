
#!/bin/bash

## install requirements
sudo apt-get install gdal-bin
sudo apt-get install postgis postgresql-devel postgresql-client
python -m venv .venv
python -m pip install wheel
python -m pip install -r requirements.txt

## change np.int to np.int64:
# l. 346 : .venv/lib/python3.9/site-packages/networkx/readwrite/graphml.py
# l.223: .venv/lib/python3.9/site-packages/networkx/readwrite/gexf.py

## kallisto needs to be in path
SOFTWARE_DIR=/software
mkdir $SOFTWARE_DIR
cd $SOFTWARE_DIR

wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
tar -xzf kallisto_linux-v0.43.1.tar.gz
PATH $PATH:/software/kallisto_linux-v0.43.1/kallisto