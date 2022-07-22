import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_dotplot(df: pd.DataFrame, out_file: str, title: str, inv_color: str = "red"):
    """Plot the dotplot from bamfile
    query and reference coordinates column names are ax, ay, bx, by.

    :param df: The dataframe with the dotplot.
    :param out_dir: The output directory.
    :param title: The title of the plot.
    :param inv_color: The color of the inverted regions.
    """
    fig, ax = plt.subplots(figsize=(11, 3))
    x_base = df.bx.min()

    for i, row in df.iterrows():
        x_coords = [row.ax, row.ay]
        y_coords = [row.bx + x_base, row.by + x_base]

        color = "black"
        if row.ay < row.ax:
            color = inv_color

        ax.plot(x_coords, y_coords, color=color)
        x_base += row.by

    ax.set_title(title)
    ax.set_xlabel("Reference")
    ax.set_ylabel("Contigs")

    plt.savefig(out_file, bbox_inches="tight")


class Bedgraph:
    """Class to store and work with bedgraph files


    Methods:
    read_bedgraph: reads bedgraph wiith coverage column.
    get_coverage_array: returns numpy array of coverage column
    get_bins: generates unzipped bins from start & end positions.
    plot_coverage: barplot of coverage by window in bdgraph.
    """

    def __init__(self, bedgraph_file):
        print("Reading {}".format(bedgraph_file))
        self.bedgraph = self.read_bedgraph(bedgraph_file)
        self.coverage = self.get_coverage_array(self.bedgraph)
        self.bins = self.get_bins(self.bedgraph)

    def read_bedgraph(self, coverage_file) -> pd.DataFrame:
        coverage = pd.read_csv(coverage_file, sep="\t", header=None).rename(
            columns={0: "read_id", 1: "start", 2: "end", 3: "coverage"}
        )

        return coverage

    def get_bins(self, bedgraph: pd.DataFrame) -> np.ndarray:
        """
        Get the bins.

        :param coverage_file: The coverage file.
        """
        bins = bedgraph.start.to_list() + bedgraph.end.to_list()
        bins = list(set(bins))
        bins = np.array(sorted(bins))

        return bins

    def get_coverage_array(self, coverage: pd.DataFrame) -> np.ndarray:
        """
        Get the coverage of the remapping.

        :param coverage_file: The coverage file.
        """
        coverage_values = np.array(coverage.coverage.to_list())

        return coverage_values

    def plot_coverage(self, output_file):
        """
        Plot the coverage of the remapping.

        :param coverage_file: The coverage file. bedgraph produced with samtools.
        :param output_file: The output file.
        """

        fig, ax = plt.subplots(figsize=(11, 3))

        # np.histogram(coverage.coverage, bins=[coverage.start, coverage.end])
        # coverage.plot(x="start", y="coverage", kind="scatter")
        if len(self.coverage) <= 1:
            print("No coverage")
            return

        ax.bar(
            x=self.bins[:-1],
            height=self.coverage,
            width=np.diff(self.bins),
            align="edge",
            fc="skyblue",
            ec="none",
        )

        ax.set_xlabel("Reference")
        ax.set_ylabel("Coverage")

        plt.savefig(output_file, bbox_inches="tight")


def get_args():
    """
    Get arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(
            description="Plot the coverage of the remapping."
        )
        parser.add_argument("--coverage_file", help="The coverage file.", required=True)
        parser.add_argument(
            "--output_file",
            help="The output file.",
            required=False,
            default="coverage.png",
        )
        args = parser.parse_args()

        return args

    except argparse.ArgumentError as e:
        print("ArgumentError: {}".format(e))
        print(e)
        sys.exit(1)


def main():
    """
    Main function.
    """

    args = get_args()

    bedgraph = Bedgraph(args.coverage_file)
    bedgraph.plot_coverage(args.output_file)


if __name__ == "__main__":
    main()
