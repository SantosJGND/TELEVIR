import os
from ftplib import FTP
from pathlib import Path

import pandas as pd


def entrez_ncbi_taxid_command(lines, tempfile, outdir, outfile):

    Path(tempfile).touch()

    with open(tempfile, "w") as ftemp:
        ftemp.write(lines)

    os.system(
        f"cat {tempfile} | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element TaxId,Title >> {outdir}{outfile}"
    )


def get_args():

    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a file with taxid and accession number"
    )

    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="Input file with accession numbers",
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="Output file with taxid and accession numbers",
    )

    parser.add_argument(
        "-t",
        "--tempfile",
        type=str,
        required=True,
        help="Temporary file to store accession numbers",
    )

    parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        required=True,
        help="Output directory to store the output file",
    )

    args = parser.parse_args()

    return args


def main():

    args = get_args()

    new_descriptions = pd.DataFrame(columns=["taxid", "description"])

    input_df = pd.read_csv(args.infile, sep="\t")

    # analyze pandas dataframe by chuncks of 5000 rows
    for i in range(0, len(input_df), 5000):

        # get the lines of the dataframe
        lines = "\n".join(input_df["acc"][i:i+5000])

        entrez_ncbi_taxid_command(
            lines, args.tempfile, args.outdir, args.outfile)

        # read the output file
        new_df = pd.read_csv(
            f"{args.outdir}{args.outfile}", sep="\t", header=None).rename(columns={0: "taxid", 1: "description"})

        # append the new dataframe to the old one
        new_descriptions = pd.concat([new_descriptions, new_df])

        # remove the output file
        os.remove(f"{args.outdir}{args.outfile}")
