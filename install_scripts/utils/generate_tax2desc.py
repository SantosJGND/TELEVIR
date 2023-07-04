import os
from ftplib import FTP
from pathlib import Path

import pandas as pd


def entrez_ncbi_taxid_command(lines, tempfile, outdir, outfile):
    Path(tempfile).touch()

    with open(tempfile, "w") as ftemp:
        ftemp.write(lines)

    os.system(
        f"cat {tempfile} | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element TaxId,Title > {outdir}{outfile}"
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
    chunksize = 1000

    # read file in chunks using pandas
    d = 0

    for chunk in pd.read_csv(args.infile, chunksize=chunksize, header=None, sep="\t"):
        print(f"Processing chunk {d}")
        chunk = chunk[0].tolist()
        chunk = "\n".join(chunk)

        entrez_ncbi_taxid_command(chunk, args.tempfile, args.outdir, args.outfile)

        new_names = pd.read_csv(
            f"{args.outdir}{args.outfile}",
            sep="\t",
            header=None,
            names=["taxid", "description"],
        )
        new_names = new_names.drop_duplicates()

        new_descriptions = pd.concat([new_descriptions, new_names])
        new_descriptions = new_descriptions.drop_duplicates()

    new_descriptions.to_csv(f"{args.outdir}{args.outfile}", sep="\t", index=False)


if __name__ == "__main__":
    main()
