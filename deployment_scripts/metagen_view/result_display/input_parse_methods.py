import os
import zipfile

import pandas as pd
import requests
from bs4 import BeautifulSoup


def get_config(path: str):
    """get config data for path.
    returns a dict with the config data.

    :param path:
    :return: config_data
    """
    with open(path) as f:
        conf = f.readlines()
        conf = [x.strip() for x in conf if "=" in x]
        conf = {x.split("=")[0]: x.split("=")[1].replace('"', "") for x in conf}

    return conf


def scrape_description(accid: str, existing_description: str = None) -> str:
    """
    Scrape the description for the relevant information.
    """
    if accid.count("_") > 1:
        accid = accid.split("_")[:-1]
        accid = "_".join(accid)
    url = f"https://www.ncbi.nlm.nih.gov/nuccore/{accid}"
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET",
        "Access-Control-Allow-Headers": "Content-Type",
        "Access-Control-Max-Age": "3600",
        "User-Agent": "Mozilla/4.0 (X11; Ubuntu; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0",
    }

    req = requests.get(url, headers)
    soup = BeautifulSoup(req.content, "html.parser")
    try:
        title = soup.find_all("div", class_="rprtheader")[0].h1.text
    except:
        title = existing_description

    if title == "nan":
        return str(existing_description)
    else:
        return str(title)


def fastqc_parse(fastqc_path: str, stdin_fastqc_name: str = "stdin_fastqc"):
    """parse fastqc.
    returns a dict with the first 10 lines of the fastqc_data.txt file

    :param fastqc_path:
    :return: fastqc_data
    """

    if not os.path.isfile(fastqc_path):
        fqreads = pd.DataFrame(columns=["measure", "value"]).set_index("measure")
        return fqreads

    with zipfile.ZipFile(fastqc_path) as zf:

        fqreads = zf.read(f"{stdin_fastqc_name}/fastqc_data.txt").decode("utf-8")
        fqreads = fqreads.split("\n")[5:10]
        fqreads = [x.split("\t") for x in fqreads]
        fqreads = pd.DataFrame(fqreads, columns=["measure", "value"]).set_index(
            "measure"
        )
        fqreads.rename({"Total Sequences": "Total_Sequences"}, axis=0, inplace=True)

    return fqreads
