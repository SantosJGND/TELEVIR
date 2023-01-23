

sudo apt-get install gdal-bin
sudo apt-get install postgis postgresql-devel postgresql-client
python -m venv .venv
python -m pip install wheel
python -m pip install -r requirements.txt

## change np.int to np.int64:
# l. 346 : .venv/lib/python3.9/site-packages/networkx/readwrite/graphml.py
# l.223: .venv/lib/python3.9/site-packages/networkx/readwrite/gexf.py
