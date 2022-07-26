import logging
import os

import pandas as pd

from pathogen_detection.object_classes import Remap_Target
from pathogen_detection.utilities import scrape_description


def process_class(r2, maxt=6):
    """
    Process classification results.
    """
    r2 = r2.drop_duplicates(subset=["taxid"], keep="first")
    r2 = r2.reset_index(drop=True)
    r2 = r2.sort_values("counts", ascending=False)

    taxids_tokeep = []
    nr2 = []

    if "length" in r2.columns:
        r2["length"] = r2["length"].astype(int)
        r2c = r2.copy().sort_values("length", ascending=False)

        for i in range(r2.shape[0]):
            if i < maxt:
                taxids_tokeep.append(r2c.taxid[i])
                nr2.append(r2c.loc[i])
                if r2.taxid.tolist()[i] not in taxids_tokeep:
                    taxids_tokeep.remove(r2.taxid[i])
                    nr2.append(r2.loc[i])
            else:

                break
    else:

        r2 = r2.head(maxt)

    if len(nr2):
        r2 = pd.concat(nr2, axis=1).T

    return r2


def merge_classes(r1, r2, maxt=6, exclude="phage"):
    """
    merge tables of taxids to columns.
    """
    if "description" in r1.columns:
        r1 = (
            r1[~r1.description.str.contains(exclude)]
            .drop_duplicates(subset=["taxid"], keep="first")
            .sort_values("counts", ascending=False)
        )

    r1 = r1[["taxid", "counts"]]

    r2pres = 1

    if len(r2):
        r2pres = 2
        if "description" in r2.columns:
            r2 = r2[~r2.description.str.contains(exclude)]

        r1.taxid = r1.taxid.astype(str)
        r2.taxid = r2.taxid.astype(str)

        shared = pd.merge(r1, r2, on=["taxid"], how="inner").sort_values(
            "counts_x", ascending=False
        )
        maxt = maxt - shared.shape[0]

        if maxt < 0:
            r1 = shared
        else:
            r2 = (
                pd.merge(r2, shared, indicator=True, how="outer")
                .query('_merge=="left_only"')
                .drop("_merge", axis=1)
            )
            r2 = process_class(r2, maxt=maxt)

            r1 = (
                pd.concat([shared, r2, r1.head(maxt)], axis=0)
                .drop_duplicates(subset=["taxid"], keep="first")
                .reset_index(drop=True)
            )

    return r1.head(maxt * r2pres)


def filter_files_to_map(nset) -> list:
    """return at most two files, give priority to refseq and virosaurus.

    args: pandas data frame of taxid, acc, files.
    """

    ref1 = "refseq_viral.genome.fna.gz"
    ref2 = "virosaurus90_vertebrate-20200330.fas.gz"

    files_to_map = []
    files_count = nset.file.value_counts()

    if ref1 in nset.file.unique():
        files_to_map.append(ref1)
    if ref2 in nset.file.unique():
        files_to_map.append(ref2)

    if len(files_to_map) == 0:
        files_to_map.append(files_count.index[0])
    if len(files_to_map) == 1 and files_count.shape[0] > 1:
        files_count = files_count[files_count.index != files_to_map[0]]
        files_to_map.append(files_count.index[0])

    return files_to_map


class Metadata_handler:
    def __init__(self, metadata_paths, sift_query: str = "phage"):
        """
        Initialize metadata handler.

        Args:
            metadata_paths: dictionary of paths to metadata files.
            sift_query: string to filter sift report.

        """
        self.metadata_paths = metadata_paths
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler())

        self.input_taxonomy_to_descriptor_path = metadata_paths[
            "input_taxonomy_to_descriptor_path"
        ]
        self.input_accession_to_taxid_path = metadata_paths[
            "input_accession_to_taxid_path"
        ]
        self.input_protein_accession_equivalent_path = metadata_paths[
            "input_protein_accession_to_protid_path"
        ]

        self.input_protein_accession_to_taxid_path = metadata_paths[
            "input_protein_accession_to_taxid_path"
        ]

        self.accession_to_taxid: pd.DataFrame
        self.taxonomy_to_description: pd.DataFrame
        self.protein_to_accession: pd.DataFrame
        self.protein_accession_to_csv: pd.DataFrame
        self.sift_query = sift_query
        self.sift_report = pd.DataFrame(
            [[0, 0, 0]], columns=["input", "output", "removed"]
        )

    def results_process(self, df: pd.DataFrame, sift: bool = True) -> pd.DataFrame:
        """
        Process results.
        merge df with metadata to create taxid columns.
        summarize merged dataframe to get counts per taxid.
        if sift is true, filter results to only include self.sift_query.
        """

        df = self.merge_report_to_metadata(df)

        df = self.map_hit_report(df)

        if sift:
            sifted_df = self.sift_report_filter(df, query=self.sift_query)
            self.sift_report = self.sift_summary(df, sifted_df)
            df = sifted_df

        return df

    def get_metadata(self):
        """
        Get metadata from files.
        """

        self.accession_to_taxid = pd.read_csv(
            self.input_accession_to_taxid_path, sep="\t", header=0
        )
        self.taxonomy_to_description = pd.read_csv(
            self.input_taxonomy_to_descriptor_path, sep="\t", header=0
        )
        self.protein_to_accession = pd.read_csv(
            self.input_protein_accession_equivalent_path, sep="\t", header=0
        )

        self.protein_accession_to_taxid = pd.read_csv(
            self.input_protein_accession_to_taxid_path, sep="\t", header=0
        )

        self.logger.info("Finished retrieving metadata")

    def merge_report_to_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """

        Args:
            df: classifier output, possessing at least columns: acc, protid, prot_acc or taxid.

        Returns:
            df: classifier output, possessing original columns plus: description.

        """

        if df.shape[0] == 0:
            return pd.DataFrame(columns=["taxid", "description", "file"])

        if "taxid" not in df.columns:

            if "prot_acc" in df.columns:
                return self.merge_check_column_types(
                    df, self.protein_accession_to_taxid, "prot_acc"
                )

            if "protid" in df.columns:
                df = self.merge_check_column_types(
                    df, self.protein_to_accession, "protid"
                )

            if "acc" in df.columns:
                df = self.merge_check_column_types(
                    df, self.accession_to_taxid, column="acc", column_two="acc_in_file"
                )

            else:
                raise ValueError(
                    "No taxid, accid or protid in the dataframe, unable to retrieve description."
                )

        return self.merge_check_column_types(df, self.taxonomy_to_description, "taxid")

    @staticmethod
    def merge_check_column_types(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        column: str,
        column_two: str = "",
    ):
        """
        Merge df1 with df2 on given column, check that the column is of the same type.
        """

        if len(column_two) == 0:
            column_two = column

        if column not in df1.columns:
            df1[column] = ""
        if column_two not in df2.columns:
            df2[column_two] = ""
        if df1[column].dtype != str:
            df1[column] = df1[column].astype(str)
        if df2[column_two].dtype != str:
            df2[column_two] = df2[column_two].astype(str)

        return pd.merge(df1, df2, left_on=column, right_on=column_two)

    @staticmethod
    def sift_report_filter(df, query: str = "phage"):
        """
        Filter df to only include hits with query in description column.
        """
        if "description" not in df.columns:
            df["description"] = ""

        ntab = df[~df.description.str.contains(query)]
        ntab = ntab.drop_duplicates(subset=["qseqid"])

        return ntab

    def sift_summary(self, merged_report: pd.DataFrame, filtered_reads: pd.DataFrame):
        """
        generate report of the difference in sequence ids between merged_report and filtered_reads.
        """

        n_removed = merged_report.taxid.nunique() - filtered_reads.taxid.nunique()

        logging.info(f"{n_removed} reads sifted from the report")
        sift_report_df = pd.DataFrame(
            [
                [
                    merged_report.taxid.nunique(),
                    filtered_reads.taxid.nunique(),
                    n_removed,
                ]
            ],
            columns=["input", "output", "removed"],
        )

        return sift_report_df

    @staticmethod
    def map_hit_report(merged_table: pd.DataFrame):
        """
        create a report of the number of hits per taxid.
        """

        if merged_table.shape[0] == 0:
            return pd.DataFrame(columns=["taxid", "description", "file", "counts"])

        counts = merged_table.taxid.value_counts()
        counts = pd.DataFrame(counts).reset_index()
        counts.columns = ["taxid", "counts"]

        new_table = pd.merge(
            left=merged_table, right=counts, on="taxid"
        ).drop_duplicates(subset="taxid")

        new_table = new_table.sort_values("counts", ascending=False)

        return new_table

    def merge_reports(
        self,
        rclass: pd.DataFrame,
        aclass: pd.DataFrame,
        max_remap: int = 15,
    ):
        """merge the reports and filter them."""

        targets = merge_classes(rclass, aclass, maxt=max_remap)
        targets.dropna(subset=["taxid"], inplace=True)
        targets["taxid"] = targets["taxid"].astype(int)

        return targets

    def generate_mapping_targets(
        self,
        targets,
        prefix: str,
        taxid_limit: int = 9,
        fasta_main_dir: str = "",
    ):
        """
        check for presence of taxid in targets in self.accession_to_taxid.
        if present, find every accession ID associated, and create a ReferenceMap object
        for each accession ID.

        """

        remap_targets = []
        remap_absent = []
        taxf = self.accession_to_taxid

        for taxid in targets.taxid.unique():

            if len(taxf[taxf.taxid == taxid]) == 0:
                remap_absent.append(taxid)

                nset = pd.DataFrame(columns=["taxid"])
                continue

            nset = (
                taxf[taxf.taxid == taxid]
                .reset_index(drop=True)
                .drop_duplicates(subset=["acc"], keep="first")
            )
            ###
            files_to_map = filter_files_to_map(nset)

            ####

            for fileset in files_to_map:

                nsu = nset[nset.file == fileset]

                if nsu.shape[0] > taxid_limit:
                    nsu = nsu.drop_duplicates(
                        subset=["taxid"], keep="first"
                    ).reset_index()

                    # with open(logd + "taxid_map.log", "a") as f:
                    #    f.write(
                    #        f"{taxid} found in db {fileset} with over 100 accs, using : {nsu.acc[0]}\n"
                    #    )

                for pref in nsu.acc.unique():

                    # with open(logd + "taxid_map.log", "a") as f:
                    #    f.write(
                    #        f"{taxid} with acc {pref} found, mapping to db: {', '.join(files_to_map)}\n"
                    #    )

                    nsnew = nsu[nsu.acc == pref].reset_index(drop=True)
                    pref_simple = (
                        pref.replace(".", "_")
                        .replace(";", "_")
                        .replace(":", "_")
                        .replace("|", "_")
                    )

                    self.taxonomy_to_description.taxid = (
                        self.taxonomy_to_description.taxid.astype(int)
                    )
                    print(self.taxonomy_to_description.head())
                    description = self.taxonomy_to_description[
                        self.taxonomy_to_description.taxid == int(taxid)
                    ].description.unique()

                    if len(description) == 0:
                        description = [""]

                    if len(description) > 1:
                        description = sorted(description, key=len)
                    print("pref: ", pref)

                    description = description[-1]
                    description = scrape_description(pref, description)

                    remap_targets.append(
                        Remap_Target(
                            pref,
                            pref_simple,
                            taxid,
                            os.path.join(fasta_main_dir, fileset),
                            prefix,
                            description,
                            [nsnew.acc_in_file[0]],
                        )
                    )

        return remap_targets, remap_absent
