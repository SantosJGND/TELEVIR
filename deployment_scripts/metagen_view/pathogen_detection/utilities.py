import os
import zipfile

import pandas as pd
import requests
from bs4 import BeautifulSoup


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
        "User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0",
    }

    try:
        req = requests.get(url, headers)
        soup = BeautifulSoup(req.content, "html.parser")
        title = soup.find_all("div", class_="rprtheader")[0].h1.text
    except:
        title = existing_description

    if title == "":

        return str(existing_description)
    else:
        return str(title)


def read_sam_coordinates(samfile: str) -> pd.DataFrame:
    """Read the sam file and return a dataframe with the coordinates."""

    try:
        df = pd.read_csv(samfile, sep="\t", header=None).rename(
            columns={0: "read_id", 7: "ax", 8: "ay", 2: "bx", 3: "by"}
        )

        df = df.sort_values("ax", ascending=True)

    except:
        df = pd.DataFrame(columns=["read_id", "ax", "ay", "bx", "by"])

    return df


def plotly_dotplot(df: pd.DataFrame):
    """Plot the dotplot using plotly."""

    if df.empty:
        return None

    import plotly.graph_objects as go
    from plotly.offline import plot

    x_base = df.bx.min()
    total_segdf = []

    for i, row in df.iterrows():
        x_coords = [row.ax, row.ay]
        y_coords = [row.bx + x_base, row.by + x_base]

        x_base += row.by
        segdf = pd.DataFrame(
            dict(x=x_coords, y=y_coords, label=[f"{i}_{row.read_id}"] * 2)
        )

        total_segdf.append(
            go.Scatter(x=x_coords, y=y_coords, mode="lines", name=f"{i}_{row.read_id}")
        )

    layout = {
        "xaxis_title": "Reference Sequence",
        "yaxis_title": "contigs",
        "legend": {"family": "Courier New, monospace", "size": 8},
        "font": {"family": "Courier New, monospace", "size": 8},
        "height": 250,
        "width": 760,
        "template": "none",
        # "colorway": px.colors.qualitative.Bold,
        "margin": dict(l=50, r=15, b=30, t=10, pad=4),
        "legend": dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    }

    plot_div = plot(
        {"data": total_segdf, "layout": layout},
        output_type="div",
    )

    return plot_div
