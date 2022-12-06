def get_args():
    """get parse arguments"""

    parser = argparse.ArgumentParser(description="parse arguments")

    parser.add_argument("--suffix", type=str, default=None, help="run suffix")
    parser.add_argument(
        "--taxd", type=str, help="directory of acc2taxid and tax2description metadata"
    )
    parser.add_argument(
        "--query", type=str, help="query string. split using ; for multiple queries"
    )
    parser.add_argument("--output", type=str, default="./", help="output directory")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        logging.info("Input error")
        print(e)
        sys.exit(1)

    return (
        args.suffix,
        args.taxd,
        args.query,
        args.output,
    )


def main():
    """main function"""

    suffix, classd, query, output = get_args()
    try:
        taxd = pd.read_csv(
            "{}taxid2desc.tsv".format(classd), sep="\t"
        )  # taxid2desc.tsv
    except FileNotFoundError:
        logging("taxid2desc.tsv not found. exiting.")
        sys.exit(1)

    try:
        seqid = pd.read_csv(
            "{}acc2taxid.tsv".format(classd), sep="\t"
        )  # taxid2desc.tsv
    except FileNotFoundError:
        logging("acc2taxid.tsv not found. exiting.")
        sys.exit(1)

    query = query.split(",")
    grep_statement = "grep "

    for q in query:
        grep_statement += f"-e {q} "

    tempfile = "temp.txt"
    tempfile = os.path.join(output, tempfile)
    grep_statement = (
        f"{grep_statement} {os.path.join(classd,'taxid2desc.tsv')} > {tempfile}"
    )
    os.system(grep_statement)

    if os.path.getsize(tempfile) == 0:
        logging.info("no hits found")
        os.system(f"rm {tempfile}")
        sys.exit(1)

    query_result = (
        pd.read_csv(tempfile, sep="\t", header=None)
        .rename(columns={0: "taxid", 1: "description"})
        .drop_duplicates(subset="taxid", keep="first")
    )
    os.system(f"rm {tempfile}")

    query_result = pd.merge(seqid, query_result, on="taxid")

    print(query_result)

    logging.info(
        f"{query_result.taxid.nunique()} taxid(s) found across {query_result.file.nunique()} file(s)."
    )

    if suffix == None:
        return
    else:
        output = os.path.join(output, "{}.report.tsv".format(suffix))
        logging.info("writing to file {}".format(output))
        query_result.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    import argparse
    import logging
    import os
    import sys

    import pandas as pd

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    main()
