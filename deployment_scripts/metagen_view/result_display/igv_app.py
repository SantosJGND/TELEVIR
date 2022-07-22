import os
import urllib.request as urlreq

import dash_bio as dashbio
from dash import html
from django_plotly_dash import DjangoDash

assets_path = os.getcwd() + "/django_plotly_dash/app/igv_app"

app = DjangoDash("igv_app")  # replaces dash.Dash
app.assets_folder = assets_path

directory = ""
bamfile = "snps.sorted.bam"
indexBam = "snps.sorted.bam.bai"
# source_data = urlreq.urlopen("https://git.io/pileup-synth4.json").read().decode("utf-8")


# source_data = {}
# with open(os.path.join(directory, bamfile), "rb") as f:
#    print(f.read().decode("utf-8"))
#    source_data = f.read()

app.layout = html.Div(
    dashbio.Pileup(
        id="tracks-pileup",
        range={
            "contig": "NODE_1_length_1228_cov_15.616368",
            "start": 140,
            "stop": 1500,
        },
        reference={
            "label": "nerw",
            "url": app.get_asset_static_url(
                "genome.2bit"
            ),  # os.path.join(directory, "genome.2bit"),
            #            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit",
        },
        tracks=[
            {
                "viz": "coverage",
                "label": "alignments",
                "source": "bam",
                "sourceOptions": {
                    "url": app.get_asset_static_url(
                        "sample.bam"
                    ),  # os.path.join(directory, bamfile),
                    "indexUrl": os.path.join(directory, indexBam),
                },
            },
            {
                "viz": "pileup",
                "label": "alignments",
                "source": "bam",
                "sourceOptions": {
                    "url": app.get_asset_static_url(
                        "sample.bam.bai"
                    ),  # os.path.join(directory, bamfile),
                    "indexUrl": os.path.join(directory, indexBam),
                },
            },
        ],
    )
)
