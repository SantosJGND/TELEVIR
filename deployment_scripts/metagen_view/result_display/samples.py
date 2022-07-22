"""
Ceated on 10/05/2022
@author: joao santos
"""
import os
import shutil
from dataclasses import dataclass
from re import S
from typing import Type
from xmlrpc.client import Boolean

import numpy as np
import pandas as pd
import requests
from metagen_view.settings import STATICFILES_DIRS
from pathogen_detection.utilities import plotly_dotplot, read_sam_coordinates

from result_display.input_parse_methods import (
    fastqc_parse,
    get_config,
    scrape_description,
)
from result_display.tables import ContigTable


class SuprunMetaData:
    """
    Class to handle the run directory.

    Attributes:
        run_dir (str): run directory.
        run_name (str): run name.
        run_type (str): run type.
    """

    run_dir: str
    name: str
    sample_name: str
    sample_type: str
    config_path: str
    logdir: str
    filtered_reads_dir: str
    sift_path: str

    assembly_root: os.PathLike
    read_summary_path: os.PathLike

    preprocess_sift: pd.DataFrame
    read_summary: pd.DataFrame

    conf: dict
    processed_reads_r1: str = ""
    processed_reads_r2: str = ""
    assembly_contigs: str = ""

    def __init__(self, sample_name: str, run_dir: str, run_name: str, sample_type: str):
        """
        Args:
            run_dir (str): run directory.
        """
        self.run_dir = run_dir
        self.name = run_name
        self.sample_name = sample_name
        self.sample_type = sample_type

        self.config_path = os.path.join(self.run_dir, "config.sh")
        self.logdir = os.path.join(self.run_dir, "logs")
        self.filtered_reads_dir = os.path.join(self.run_dir, "reads", "hd_filtered")
        self.sift_path = os.path.join(
            self.filtered_reads_dir,
            os.path.basename(self.run_dir) + ".sift",
        )
        self.assembly_root = os.path.join(self.run_dir, "assembly", "output")
        self.read_summary_path = os.path.join(self.logdir, "reads_latest.stats")

        self.preprocess_sift = self.read_FilteringSift()
        self.read_summary = self.get_read_summary()

    def get_assembly_contigs(self):
        """
        Get the contigs from the assembly.
        """
        assembly_result = os.path.join(self.assembly_root, self.conf["ASSEMBLY_SOFT"])
        assembly_file = [x for x in os.listdir(assembly_result) if x.endswith(".gz")]
        if len(assembly_file) == 0:
            self.assembly_contigs = ""
            return
        else:
            assembly_file = assembly_file[0]

        assembly_file = os.path.join(assembly_result, assembly_file)

        static_store = os.path.join(STATICFILES_DIRS[0], "assemblies")
        if not os.path.exists(static_store):
            os.makedirs(static_store)

        if os.path.exists(assembly_file):
            # output = os.path.join(static_store, os.path.basename(assembly_file))
            output = ".".join([self.sample_name, self.name, "assembly.fasta.gz"])

            shutil.copy(
                assembly_file,
                os.path.join(static_store, output),
            )

            self.assembly_contigs = os.path.basename(output)

        else:
            print(
                f"{assembly_file} does not exist. Please check the path and try again."
            )

    @staticmethod
    def store_reads(reads_path: str, output=""):

        if not output:
            output = os.path.basename(reads_path)

        static_store = os.path.join(STATICFILES_DIRS[0], "depleted_reads")
        if not os.path.exists(static_store):
            os.makedirs(static_store)

        if os.path.exists(reads_path):

            shutil.copy(reads_path, os.path.join(static_store, output))

            return os.path.basename(output)

        else:
            print(f"{reads_path} does not exist. Please check the path and try again.")
            return ""

    def get_processed_reads(self):
        """
        store the processed reads in a file.
        """

        if self.conf["ENRICH"] != "true" and self.conf["ENRICH"] != "true":
            return

        reads_fofn = os.path.join(self.logdir, f"{self.name}_latest.fofn")
        if not os.path.exists(reads_fofn):
            return

        reads_fofn = pd.read_csv(
            reads_fofn,
            sep="\t",
            header=None,
        ).iloc[:, 0]

        output_r1 = ".".join(
            [self.sample_name, self.name, "processed_reads_r1.fastq.gz"]
        )
        processed_reads_r1 = os.path.join(
            self.filtered_reads_dir, os.path.basename(reads_fofn[0])
        )

        self.processed_reads_r1 = self.store_reads(processed_reads_r1, output=output_r1)

        if self.sample_type == "PE":

            output_r2 = ".".join(
                [self.sample_name, self.name, "processed_reads_r2.fastq.gz"]
            )
            processed_reads_r2 = os.path.join(
                self.filtered_reads_dir, os.path.basename(reads_fofn[1])
            )

            self.processed_reads_r2 = self.store_reads(
                processed_reads_r2, output=output_r2
            )

    def find_children(self, confch: str = "confch"):
        """
        Find the children of the run directory.

        :param confchild: confchild (str): config file suffix.
        :return: None

        """
        children = os.listdir(self.run_dir)
        children = [x for x in children if confch in x]
        children = [os.path.splitext("_".join(x.split("_")[1:]))[0] for x in children]

        self.children = children

    def get_read_summary(self, fname: str = "reads_latest.stats") -> pd.DataFrame:
        """get read summary data for path.
        returns a dict with the read summary data. sets counts to 0 if file is missing or empty.

            :param path:
            :return: read_summary_data

        """

        try:
            read_summary = (
                (
                    pd.read_csv(self.read_summary_path, sep="\t", header=None).rename(
                        columns={0: "step", 1: "value"}
                    )
                )
                .drop_duplicates(subset="step")
                .set_index("step")
            )
        except:
            read_summary = pd.DataFrame(
                {"step": ["CLEAN", "ENRICHMENT"], "value": [0, 0]}
            )

        if len(read_summary) == 0:
            read_summary = pd.DataFrame(
                {"step": ["CLEAN", "ENRICHMENT"], "value": [0, 0]}
            )

        return read_summary

    def get_contig_summary(self, fname: str = "filter_fasta.log"):
        """get contig summary data for path.

        :param path:

        """

        path = os.path.join(self.logdir, fname)

        try:
            contig_summary = pd.read_csv(path, sep="\t", header=None).rename(
                columns={0: "contig", 1: "length"}
            )
        except:
            contig_summary = pd.DataFrame(
                columns=["contig", "length"],
            )

        self.assembly_summary = contig_summary

    def read_FilteringSift(self):
        """
        read FilteringSift.log file
        """

        try:

            sift_data = pd.read_csv(self.sift_path, sep=",")

        except:
            sift_data = pd.DataFrame(
                [[0, 0, 0]], columns=["input", "output", "removed"]
            )

        return sift_data


class SampleMetaData:
    """
    Class to handle the sample directory.

    Attributes:
        sample_dir (str): sample directory.
        sample_name (str): sample name.
        sample_type (str): sample type.
        sample_combinations (str): sample combinations.
        sample_input (str): sample input.
        sample_technology (str): sample technology.
        sample_report (str): sample report.
        self.supruns (list): list of super runs.
        self.runs (list): list of runs.
    """

    project_name: str
    sample_dir: str
    sample_name: str
    sample_name_simple: str
    qc_dir: str

    reads_clean_dir: os.PathLike
    input_fastqc_report: os.PathLike
    processed_fastqc_report: os.PathLike
    qc_config: os.PathLike
    input_fofn: os.PathLike
    read_tech_file: os.PathLike
    comb_file: os.PathLike
    processed_fastq_r1: os.PathLike
    processed_fastq_r2: os.PathLike

    qcdata: dict
    qc_soft: str
    sample_type: str
    sample_combinations: int
    sample_input: str
    sample_technology: str

    def __init__(self, sample_dir: str, sample_name: str, project_name: str):
        """
        Args:
            sample_dir (str): sample directory.
        """
        self.project_name = project_name
        self.sample_dir = sample_dir
        self.sample_name = sample_name
        self.sample_name_simple = os.path.splitext(sample_name)[0]
        self.qc_dir = os.path.join(self.sample_dir, "sampleQC")
        self.reads_clean_dir = os.path.join(self.qc_dir, "reads", "clean")

        self.input_fastqc_report = os.path.join(
            self.reads_clean_dir, "input_data.fastqc.html"
        )
        self.processed_fastqc_report = os.path.join(
            self.reads_clean_dir, "processed_data.fastqc.html"
        )
        self.qc_config = os.path.join(self.qc_dir, "config.sh")
        self.input_fofn = os.path.join(self.sample_dir, "input.fofn")
        self.read_tech_file = os.path.join(self.sample_dir, "technology.txt")
        self.comb_file = os.path.join(self.sample_dir, "parameter_combinations.txt")
        self.supruns: list = []
        self.runs: list = []

    def QC_summary(self, QCdir: str):
        """get QC summary
        runs fastqc parse on zip files in QCdir

        :param QCdir:
        :return: qc_summary
        """

        qc_summary = {
            "input": fastqc_parse(os.path.join(QCdir, "input_data.fastqc.zip")),
            "processed": fastqc_parse(os.path.join(QCdir, "processed_data.fastqc.zip")),
        }

        return qc_summary

    def get_qc_data(self, reads_dir: str = "reads/clean/"):
        """get qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        self.qcdata = self.QC_summary(self.reads_clean_dir)
        confqc = get_config(self.qc_config)
        self.qc_soft = confqc["QC"]

    def get_input_data(self):
        """get input data for sample_class. Update sample_class.input_data.

        :param sample_class:
        :return: None
        """

        if os.path.exists(self.input_fofn):
            input_list = pd.read_csv(
                self.input_fofn,
                sep="\t",
                header=None,
            ).iloc[:, 0]

            self.sample_type = ["SE", "PE"][int(len(input_list) > 1)]

            self.processed_fastq_r1 = os.path.join(
                self.reads_clean_dir, input_list[0].replace(".fq.gz", ".fastq.gz")
            )

            if self.sample_type == "PE":
                self.processed_fastq_r2 = os.path.join(
                    self.reads_clean_dir,
                    input_list[1].replace(".fq.gz", ".fastq.gz"),
                )

            self.sample_input = ", ".join([os.path.basename(x) for x in input_list])

        else:
            self.sample_input = "NA"
            self.sample_type = "NA"

    def get_read_technology(self):
        """get read technology for sample_class. Update sample_class.read_technology.

        :param sample_class:
        :return: None
        """

        if os.path.exists(self.read_tech_file):
            with open(self.read_tech_file, "r") as f:
                self.sample_technology = f.readline().strip().split("\t")[1]

        else:
            self.sample_technology = "NA"

    def get_sample_combinations(self):
        """get sample combinations for sample_class. Update sample_class.sample_combinations.

        :param sample_class:
        :return: None
        """

        if os.path.exists(self.comb_file):

            combinations = pd.read_csv(
                self.comb_file,
                sep="\t",
                header=None,
            )

            self.sample_combinations = combinations[1].prod()

        else:
            self.sample_combinations = 0

    def child(self, suprun: Type[SuprunMetaData], tag: str, suffix: str = "run_"):
        """
        create instance of runMetaData class.
        """
        name = suffix + str(len(self.runs))
        child = RunMetaData(self.project_name, self.sample_name, suprun, name, tag)
        child.sample = self.sample_name

        child.processed_reads = int(
            self.qcdata["processed"].to_dict()["value"]["Total_Sequences"]
        )
        child.processed_reads_r1 = suprun.processed_reads_r1
        child.processed_reads_r2 = suprun.processed_reads_r2
        child.assembly_contigs = suprun.assembly_contigs

        self.runs.append(child)

    def get_run(self, run_name: str):
        """
        get run by name.

        :param run_name: run name.
        :return: RunMetaData object.
        """

        for run in self.runs:
            if run.name == run_name:
                return run
        return None


class RunMetaData:
    """
    Class to handle the run directory.

    Attributes:
        tag (str): run tag.
        conf (str): run config.
        sample (str): sample name.
        name (str): run name.
        run_dir (str): run directory.
        logdir (str): log directory.
        processed_reads (int): processed reads.
        software (dict): software used.
        params (dict): run parameters.
        runtime (int): run time.
        finished (bool): finished.
        report (pd.DataFrame): report data.
        aclass_summary (pd.DataFrame): assembly summary data.
        aclass_success (bool): assembly success.
        rclass_summary (pd.DataFrame): read summary data.
        rclass_success (bool): read success.
        remap_plan (pd.DataFrame): remap plan data.
        remap_success (bool): remap success.

    """

    def __init__(
        self,
        project_name: str,
        sample_name: str,
        suprun: Type[SuprunMetaData],
        run_name: str,
        tag: str,
    ):
        """
        Args:

            run_dir (str): run directory.
        """
        self.tag = tag
        self.conf = []
        self.sample = sample_name
        self.project_name = project_name
        self.prefix = f"{self.project_name}_{self.sample}_{self.tag}"

        self.config_path = os.path.join(suprun.run_dir, "confch_" + tag + ".sh")
        self.param_file_path = os.path.join(suprun.run_dir, tag + ".args.tsv")
        self.params_path_exists = os.path.exists(self.param_file_path)

        self.runtime_path = os.path.join(suprun.run_dir, tag + ".runtime")
        self.final_report_path = os.path.join(suprun.run_dir, tag + ".final.report.tsv")
        self.final_report_exists = os.path.exists(self.final_report_path)

        self.read_classification_path = os.path.join(
            suprun.run_dir, "classification", "reads", f"{tag}.report.tsv"
        )
        self.read_classification_exists = os.path.exists(self.read_classification_path)

        self.assembly_classification_path = os.path.join(
            suprun.run_dir, "classification", "assembly", "assembly.report.tsv"
        )
        self.assembly_classification_exists = os.path.exists(
            self.assembly_classification_path
        )

        self.merged_log_path = os.path.join(suprun.run_dir, "logs", tag, "targets.tsv")
        self.merged_log_exists = os.path.exists(self.merged_log_path)

        self.remap_plan_path = os.path.join(
            suprun.run_dir, "logs", tag, "taxid_map.log"
        )
        self.remap_plan_exists = os.path.exists(self.remap_plan_path)

        self.remap_report_path = os.path.join(
            suprun.run_dir, "logs", tag, "remap_acc_plan.log"
        )
        self.remap_report_exists = os.path.exists(self.remap_report_path)

        self.suprun = suprun
        self.name = run_name
        self.processed_reads = 0
        self.params = {}
        self.software = {}
        self.runtime = 0
        self.finished = False
        self.report = ""
        self.aclass_summary = ""
        self.aclass_success = False
        self.rclass_summary = ""
        self.rclass_success = False
        self.remap_plan = ""
        self.remap_success = False
        self.minimum_coverage = 0
        self.coverage_max = 1e5
        self.references = []

        self.pprec_sift = self.read_preprocessing_Sift()
        self.context = {}

    @staticmethod
    def store_classification_file(file_path, output: str = ""):

        classification_store = os.path.join(
            STATICFILES_DIRS[0], "classification_reports"
        )

        if not os.path.exists(classification_store):
            os.mkdir(classification_store)

        if not output:
            output = os.path.basename(file_path)

        # print(f"Storing {output}")
        # print(f"{file_path}")

        if not os.path.isfile(os.path.join(classification_store, output)):
            shutil.copy(file_path, os.path.join(classification_store, output))

        return os.path.join(classification_store, output)

    def store_classification_results(self):
        """
        move classification files to static
        """

        if self.read_classification_exists:
            output = f"{self.prefix}.read_classification.tsv"
            output = self.store_classification_file(
                self.read_classification_path, output
            )
            self.read_classification_path = output

        if self.assembly_classification_exists:

            output = f"{self.prefix}.assembly_classification.tsv"
            self.store_classification_file(self.assembly_classification_path, output)
            self.assembly_classification_path = output

        if self.merged_log_exists:
            output = f"{self.prefix}.targets.tsv"
            output = self.store_classification_file(self.merged_log_path, output)
            self.merged_log_path = output

        if self.remap_plan_exists:
            output = f"{self.prefix}.taxid_map.log"
            output = self.store_classification_file(self.remap_plan_path, output)
            self.remap_plan_path = output

        if self.params_path_exists:
            output = f"{self.prefix}.args.tsv"
            # print(self.param_file_path)
            output = self.store_classification_file(self.param_file_path, output)
            self.params_file_path = output

    def read_preprocessing_Sift(self):
        """
        read ClassSift.log file
        """
        sift_path = os.path.join(
            self.suprun.run_dir, "classification", "reads", self.tag + ".sift"
        )

        try:

            sift_data = pd.read_csv(sift_path, sep="\t")

        except:
            sift_data = pd.DataFrame(
                [[self.processed_reads, self.processed_reads, 0]],
                columns=["input", "output", "removed"],
            )

        return sift_data

    @staticmethod
    def process_id(accid: str):
        """
        Process the accession id.

        :param accid:
        :return: modified accid
        """
        if "|" in accid:
            accid = accid.split("|")[2]

        accid = accid.split(":")[0]

        return accid

    def process_report_get_sucess(self):
        """create success column based on rclass and aclass columns"""

        self.report["success"] = self.report[["rclass", "aclass"]].fillna(0).sum(axis=1)
        self.report["success_verbose"] = self.report["success"].apply(
            lambda x: ["none", "reads or contigs", "both"][x]
        )

    def process_report_id(self):
        """extract ID from contig string, necessary du to db specific formatting (kraken2)"""

        self.report["ID"] = self.report["ID"].apply(self.process_id)

    def process_report_description(self):
        """parse description from ncbi using accession id for each line in the report"""
        new_descriptions = []
        try:
            for index, row in self.report.iterrows():
                new_descriptions.append(
                    scrape_description(row["ID"], row["description"])
                )
            self.report["description"] = new_descriptions
        except requests.exceptions.ConnectionError as e:
            print(e)

    @staticmethod
    def sort_groups_report(df):
        """
        Sort the groups in the report.

        :param df:
        :return:
        """
        group_dict = {}

        for taxid in df.taxid.unique():
            loc = df.loc[df.taxid == taxid]
            loc.sort_values(
                by=[
                    "coverage",
                    "Hdepth",
                ],
                ascending=False,
                inplace=True,
            )
            meanC = loc.coverage.mean()
            group_dict[taxid] = {
                "mean": meanC,
                "data": loc,
            }

        taxid_order = sorted(
            group_dict.keys(), key=lambda x: group_dict[x]["mean"], reverse=True
        )

        new_df = pd.concat(
            [group_dict[x]["data"] for x in taxid_order],
        ).reset_index(drop=True)
        return new_df

    def process_report_match_remap_metadata(self):
        """match report data with remap metadata, create columns containing the remap metadata

        new columns:
        bam_path: path to the bam file.
        bai_path: path to the bai file.
        reference_path: path to the reference file.
        reference_index_name: path toreference fasta indexfile.
        covplot: path to coverage plot.
        refa_dotplot: path to dotplot of reference assembly.
        contig_name: contig name.
        contig_length: contig length.
        contig_depth: coverage depth across contig.
        reference_contig: reference fasta string of contig.
        refa_dotplot_exists: boolean if dotplot exists.
        covplot_exists: boolean if covplot exists.
        """

        new_columns = {
            "covplot": None,
            "refa_dotplot": None,
            "contig_name": None,
            "contig_length": None,
            "reference_contig": None,
            "contig_string": None,
            "refa_dotplot_exists": False,
            "covplot_exists": False,
            "bam_path": "None",
            "bai_path": "None",
            "reference_path": "None",
            "reference_index_path": "None",
            "reference_assembly_paf": "None",
            "mapped_scaffolds_path": "None",
            "mapped_scaffolds_index_path": "None",
            "refa_dotplot_path": "None",
        }

        refa_dotplots = {x.reference: new_columns.copy() for x in self.references}

        for i in self.references:

            refa_dotplots[i.reference]["bam_path"] = i.bam_path
            refa_dotplots[i.reference]["bai_path"] = i.bam_index_path
            refa_dotplots[i.reference]["reference_path"] = i.reference_path
            refa_dotplots[i.reference]["reference_index_path"] = i.reference_index_path
            refa_dotplots[i.reference]["reference_assembly_paf"] = i.ref_to_assembly_paf

            refa_dotplots[i.reference]["mapped_scaffolds_path"] = (
                i.mapped_scaffolds_path + ".gz"
            )
            refa_dotplots[i.reference]["mapped_scaffolds_index_path"] = (
                i.mapped_scaffolds_path + ".fai"
            )

            refa_dotplots[i.reference]["covplot"] = i.covplot_path
            refa_dotplots[i.reference]["covplot_exists"] = i.covplot_exists
            refa_dotplots[i.reference]["refa_dotplot_path"] = i.refa_dotplot_path
            refa_dotplots[i.reference]["refa_dotplot_exists"] = i.refa_dotplot_exists
            refa_dotplots[i.reference]["contig_name"] = i.contig_name
            refa_dotplots[i.reference]["contig_length"] = i.contig_len
            refa_dotplots[i.reference]["contig_string"] = i.reference_contig

            # print(refa_dotplots[i.reference])

        for new_col in new_columns.keys():
            # print(new_col)
            self.report[new_col] = self.report["simple_id"].apply(
                lambda x: refa_dotplots[x][new_col]
            )

        # print(self.report["mapped_scaffolds_path"])

    def process_report_clean(self):
        self.report = self.report[self.report.coverage > 0]

    def process_report(self):
        """
        Process the report file.
        """
        self.process_report_clean()
        self.process_report_get_sucess()
        self.process_report_id()
        self.process_report_description()
        self.process_report_match_remap_metadata()

        if self.report.shape[0] > 0:
            self.report = self.sort_groups_report(self.report)

    def get_assembly_stats(self):
        """Assembly contig summary_stats"""

        if self.suprun.assembly_summary.shape[0] > 0:
            assembly_min = self.suprun.assembly_summary["length"].fillna(0).min()
            assembly_max = self.suprun.assembly_summary["length"].fillna(0).max()
            assembly_mean = self.suprun.assembly_summary["length"].fillna(0).mean()
        else:
            assembly_min = 0
            assembly_max = 0
            assembly_mean = 0

        return assembly_min, assembly_max, assembly_mean

    def get_classification_input(self):
        """
        Get the classification read number input for classification.
        if sift was performed, then final processing is grabbed from the output ofsift.
        otherwise, final_processed_reads == post_processed_reads.

        """
        post_processed_reads = self.suprun.read_summary.loc["CLEAN"][0]
        if "ENRICHMENT" in self.suprun.read_summary.index:

            post_processed_reads = self.suprun.read_summary.loc["ENRICHMENT"][0]

        final_processing_reads = post_processed_reads.copy()
        if self.suprun.preprocess_sift.loc[0]["input"]:

            post_processed_reads = (
                self.suprun.preprocess_sift.loc[0]["output"]
                + self.suprun.preprocess_sift.loc[0]["removed"]
            )

        return post_processed_reads, final_processing_reads

    def get_context(self):
        """create dictionary of details for run."""

        processed_reads = self.processed_reads
        post_processed_reads, final_processing_reads = self.get_classification_input()

        post_percent = (int(post_processed_reads) / processed_reads) * 100
        final_processing_percent = (final_processing_reads / processed_reads) * 100

        assembly_number = self.suprun.assembly_summary.shape[0]

        assembly_min, assembly_max, assembly_mean = self.get_assembly_stats()

        #### ######
        minhit_assembly = self.aclass_summary["counts"].min()
        if not minhit_assembly or not self.aclass_summary.shape[0]:
            minhit_assembly = 0

        minhit_reads = self.rclass_summary["counts"].min()
        if np.isnan(minhit_reads):
            minhit_reads = 0

        max_gaps = self.report.ngaps.max()
        if np.isnan(max_gaps):
            max_gaps = 0

        max_prop = self.report.ref_prop.max()
        if np.isnan(max_prop):
            max_prop = 0

        max_mapped = self.report.mapped.max()
        if np.isnan(max_mapped):
            max_mapped = 0

        self.run_detail_report = Run_detail_report(
            self.report.Hdepth.max(),
            self.report.HdepthR.max(),
            max_gaps,
            max_prop,
            max_mapped,
            f"{processed_reads:,}",
            f"{post_processed_reads:,}",
            f"{post_percent:.2f}",
            self.suprun.conf["PHAGE_DEPL"] == "true",
            self.conf["PHAGE_DEPL"] == "true",
            f"{self.suprun.preprocess_sift.loc[0]['removed']:,}",
            f"{final_processing_reads:,}",
            round(final_processing_percent, 3),
            self.remap_success,
            self.remap_plan.taxid.nunique(),
            ", ".join(self.remap_plan["file"].unique()),
        )

        self.contig_classification_results = Contig_classification_results(
            True,
            self.suprun.conf["ASSEMBLE_CLASS"],
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_success,
        )

        test_db = pd.read_csv(self.read_classification_path, sep="\t")

        self.read_classification_results = Read_classification_results(
            True,
            self.conf["CLASSM"],
            self.rclass_summary.shape[0],
            minhit_reads,
            self.rclass_success,
        )

        self.assembly_report = Assembly_results(
            True,
            self.suprun.conf["ASSEMBLY_SOFT"],
            assembly_number,
            f"{assembly_min:,}",
            f"{int(assembly_mean):,}",
            f"{assembly_max:,}",
            f"{int(self.suprun.conf['ASSEMBLY_LTRIM']):,}",
        )

        self.remap_main = Remap_main(
            True,
            self.conf["REMAP_SOFT"],
            self.remap_plan["status"].value_counts().get("success", 0),
            self.minimum_coverage,
            self.coverage_max,
            self.report.shape[0],
        )

        context = {
            "sample": self.sample.split(".")[0],
            "run_name": self.name,
            "final_data": self.report,
            "report_references": {
                "max_depth": self.report.Hdepth.max(),
                "max_depthR": self.report.HdepthR.max(),
                "max_gaps": max_gaps,
                "max_prop": max_prop,
                "max_mapped": max_mapped,
            },
            "input": f"{processed_reads:,}",
            "post_input": f"{post_processed_reads:,}",
            "post_percent": round(post_percent, 2),
            "processing_final": f"{final_processing_reads:,}",
            "processing_final_percent": round(final_processing_percent, 3),
            "sift": {
                "preproc": self.suprun.conf["PHAGE_DEPL"] == "true",
                "pprc_removed": f"{self.suprun.preprocess_sift.loc[0]['removed']:,}",
                "remap": self.conf["PHAGE_DEPL"] == "true",
            },
            "merged": {
                "performed": self.remap_success,
                "merged_number": self.remap_plan.taxid.nunique(),
                "files": ", ".join(self.remap_plan["file"].unique()),
            },
            "class_contigs": {
                "performed": True,
                "assembly_class_soft": self.suprun.conf["ASSEMBLE_CLASS"],
                "assembly_class_number": self.aclass_summary.shape[0],
                "assembly_min_hit": minhit_assembly,
                "success": self.aclass_success,
            },
            "assembly": {
                "performed": True,
                "assembly_soft": self.suprun.conf["ASSEMBLY_SOFT"],
                "assembly_number": assembly_number,
                "assembly_min": f"{assembly_min:,}",
                "assembly_mean": f"{int(assembly_mean):,}",
                "assembly_max": f"{assembly_max:,}",
                "assembly_trim": f"{int(self.suprun.conf['ASSEMBLY_LTRIM']):,}",
            },
            "class_reads": {
                "performed": True,
                "reads_class_soft": self.conf["CLASSM"],
                "reads_class_number": self.rclass_summary.shape[0],
                "reads_class_min_hit": minhit_reads,
                "success": self.rclass_success,
            },
            "remap": {
                "performed": True,
                "remap_soft": self.conf["REMAP_SOFT"],
                "found_total": self.remap_plan["status"]
                .value_counts()
                .get("success", 0),
                "cov_min": self.minimum_coverage,
                "cov_max": self.coverage_max,
                "success": self.report.shape[0],
                "plotly_dotplots": {
                    x.reference: x.plotly_dotplot for x in self.references
                },
            },
            "enrichment": True,
            "depletion": False,
            "software": self.software.set_index("module").to_dict()["software"],
        }

        return context

    def generate_context(self):
        """Generate context for run."""

        self.context = self.get_context()


@dataclass(frozen=True)
class Run_detail_report:
    max_depth: int
    max_depthR: int
    max_gaps: int
    max_prop: int
    max_mapped: int
    input: str
    processed: str
    processed_percent: str
    sift_preproc: bool
    sift_remap: bool
    sift_removed_pprc: str
    processing_final: str
    processing_final_percent: str
    merged: bool
    merged_number: int
    merged_files: str


@dataclass(frozen=True)
class Contig_classification_results:
    performed: bool
    method: str
    classification_number: int
    classification_minhit: int
    classification_success: bool


@dataclass(frozen=True)
class Read_classification_results:
    performed: bool
    method: str
    classification_number: int
    classification_minhit: int
    success: bool


@dataclass(frozen=True)
class Remap_main:
    performed: bool
    success: bool
    method: str
    found_total: int
    coverage_min: int
    coverage_max: int


@dataclass(frozen=True)
class Assembly_results:
    performed: bool
    assembly_soft: str
    assembly_number: int
    assembly_min: int
    assembly_mean: int
    assembly_max: int
    assembly_trim: int


class Reference_Map:
    """
    Class to handle the reference map.

    Attributes:

        ref_map (str): reference map.
        ref_map_path (str): reference map path.
        ref_map_name (str): reference map name.
        ref_map_type (str): reference map type.
        ref_map_data (dict): reference map data.

    """

    def __init__(
        self,
        project: str,
        sample: str,
        run_name: str,
        reference: str,
        taxid: str,
        refdb: str,
        refmap_dir: str,
        refmap_logdir: dict,
    ):
        """
        Args:

            ref_map (str): reference map.
        """
        project_name = project.split(".")[0]
        self.sample = sample
        self.run_name = run_name
        self.reference = reference
        self.refdb = refdb
        self.taxid = taxid
        self.dir = refmap_dir
        self.logdir = refmap_logdir
        self.prefix = f"{project_name}_{sample}_{run_name}_{reference}_{taxid}"

        self.refA_dotplot_path = os.path.join(self.dir, "minimap2", "dotplot.png")
        self.refA_dotplot_exists = True

        self.mapped_reads = 0
        self.scaffold_remap_stats = pd.DataFrame(columns=["scaffold", "remap_stats"])
        self.reference_contig = "NA"
        self.ref_map_data = self.get_ref_map_data()
        self.bedgraph = self.get_bedgraph()
        self.remap_stats = self.get_remap_stats()
        self.mapped_reads = self.get_mapped_reads_total()
        self.get_refA_coverage()
        self.get_refA_dotplot()
        self.reference_contig_lengths()
        self.get_ref_coverage_plot()
        self.plotly_dotplot = self.get_refA_plotly_map()
        self.reference_contig = self.read_reference_contig()
        self.bam_path = self.get_bam_file()
        self.bam_index_path = self.get_bam_index_file()
        self.reference_path = self.get_reference_file()
        self.reference_index_path = self.get_reference_index_file()
        self.ref_to_assembly_paf = self.ref_assembly_map()
        self.mapped_scaffolds_path = self.get_mapped_scaffolds()

        if self.bedgraph.shape[0] > 0 and self.covplot_exists == False:
            from result_display.plot_coverage import Bedgraph

            bedgraph_path = os.path.join(self.logdir, "snps.ref.sorted.bedgraph")

            output = ".".join(
                [
                    self.sample,
                    self.run_name,
                    self.refdb,
                    self.reference,
                    "coverage",
                    "png",
                ]
            )

            covplot_destination = os.path.join(
                STATICFILES_DIRS[0], "ref_coverage", output
            )

            bedgraph = Bedgraph(bedgraph_path)
            bedgraph.plot_coverage(covplot_destination)
            self.covplot_exists = True
            self.covplot_path = os.path.basename(covplot_destination)

    def get_mapped_scaffolds(self):
        """Get mapped scaffolds."""
        scaffolds_path = os.path.join(self.dir, "assembly.fasta")

        if os.path.isfile(scaffolds_path + ".gz"):

            return scaffolds_path
        else:
            return "None"

    def ref_assembly_map(self):
        """get reference assembly map"""

        map_path = os.path.join(
            self.dir, "minimap2", os.path.basename(self.dir) + ".paf"
        )

        if os.path.isfile(map_path) and os.path.getsize(map_path) > 0:
            output = self.prefix + ".paf"
            shutil.copy(
                map_path,
                os.path.join(
                    STATICFILES_DIRS[0],
                    "igv_files",
                    output,
                ),
            )
            map_path = output

            return map_path
        else:
            return "None"

    def get_reference_file(self):
        """Get reference file."""
        reference_path = os.path.join(self.dir, "ref.fa")
        if os.path.isfile(reference_path) and os.path.getsize(reference_path):
            output = self.prefix + "_ref.fa"
            shutil.copy(
                reference_path,
                os.path.join(
                    STATICFILES_DIRS[0],
                    "igv_files",
                    output,
                ),
            )

            reference_path = output

            return reference_path
        else:
            return "None"

    def get_reference_index_file(self):
        """Get reference bai file."""
        reference_bai_path = os.path.join(self.dir, "ref.fa.fai")
        if os.path.isfile(reference_bai_path) and os.path.getsize(reference_bai_path):

            output = self.prefix + "_ref.fa.fai"
            shutil.copy(
                reference_bai_path,
                os.path.join(
                    STATICFILES_DIRS[0],
                    "igv_files",
                    output,
                ),
            )

            reference_bai_path = output

            return reference_bai_path
        else:
            return "None"

    def get_bam_file(self):
        """Get reference A bam file."""

        refA_bam = os.path.join(self.dir, "snps.sorted.bam")
        # print("refA_bam", refA_bam)

        if os.path.isfile(refA_bam) and os.path.getsize(refA_bam):

            output = self.prefix + "_refA.bam"
            shutil.copy(
                refA_bam,
                os.path.join(
                    STATICFILES_DIRS[0],
                    "igv_files",
                    output,
                ),
            )

            refA_bam = output

            return refA_bam
        else:
            return "None"

    def get_bam_index_file(self):
        """Get reference A bai file."""

        refA_bai = os.path.join(self.dir, "snps.sorted.bam.bai")

        if os.path.isfile(refA_bai) and os.path.getsize(refA_bai):
            output = self.prefix + "_refA.bam.bai"

            shutil.copy(
                refA_bai,
                os.path.join(
                    STATICFILES_DIRS[0],
                    "igv_files",
                    output,
                ),
            )

            refA_bai = output

            return refA_bai
        else:
            return "None"

    def get_mapped_reads_total(self):
        """Get total mapped reads."""

        mapped_reads = 0

        reads_kept_file = os.path.join(self.dir, "snps.reads.keep")
        if os.path.isfile(reads_kept_file):
            if os.path.getsize(reads_kept_file):
                mapped_reads = pd.read_csv(
                    reads_kept_file, sep="\t", header=None
                ).shape[0]

        return mapped_reads

    def read_reference_contig(self):
        contig_string = os.path.join(self.logdir, "ref_contig_string.txt")
        reference_contig = "NA"
        if os.path.isfile(contig_string):
            with open(contig_string, "r") as f:
                reference_contig = f.readline().strip().replace(">", "")

        return reference_contig

    def get_ref_map_data(self):
        """get reference map data."""
        read_map_file = os.path.join(self.logdir, "ref_reads_map.tsv")

        if os.path.isfile(read_map_file):
            if os.path.getsize(read_map_file):
                ref_map_data = pd.read_csv(
                    read_map_file,
                    sep="\t",
                    header=None,
                ).rename(columns={0: "reads", 1: "scaffold", 2: "quality"})
                self.mapped_reads = ref_map_data.reads.nunique()
            else:
                ref_map_data = pd.DataFrame(columns=["reads", "scaffold", "quality"])

        else:
            ref_map_data = pd.DataFrame(columns=["reads", "scaffold", "quality"])

        return ref_map_data

    def get_bedgraph(self):

        try:
            bedgraph = pd.read_csv(
                os.path.join(self.logdir, "snps.ref.sorted.bedgraph"),
                sep="\t",
                header=None,
            ).rename(columns={0: "chrom", 1: "start", 2: "end", 3: "value"})

        except:
            bedgraph = pd.DataFrame(columns=["chrom", "start", "end", "value"])

        return bedgraph

    def get_remap_stats(self):

        try:
            remap_stats = pd.read_csv(
                os.path.join(self.logdir, "scaffold_remap.stats"),
                sep="\t",
            )

            coverage_col = [x for x in remap_stats.columns if "cov" in x][0]
            remap_stats = remap_stats.rename(columns={coverage_col: "coverage"})

        except:
            remap_stats = pd.DataFrame(
                columns=[
                    "ID",
                    "Hdepth",
                    "HdepthR",
                    "tag",
                    "aclass",
                    "rclass",
                    "coverage",
                    "ngaps",
                    "ref_prop",
                    "mapped",
                ]
            )

        return remap_stats

    def reference_contig_lengths(self):
        """read contig lengths"""

        fname = os.path.join(self.dir, "snps.contig.len")
        lenf_exists = os.path.isfile(fname)

        contig_len = 0
        contig_name = "john doe"

        if lenf_exists:
            contig_length = pd.read_csv(fname, sep="\t", header=None)

            contig_len = contig_length.iloc[0, 1]
            contig_name = contig_length.iloc[0, 0]

        self.contig_name = contig_name
        self.contig_len = contig_len

    def get_refA_coverage(self):
        """get coverage plots"""

        coverage_exists = os.path.isfile(
            os.path.join(self.dir, "minimap2", "coverage.png")
        )

        if coverage_exists:
            coverage_plot = os.path.join(self.dir, "minimap2", "coverage.png")
        else:
            coverage_plot = None

        self.coverage_plot_path = coverage_plot

    def get_refA_plotly_map(self):
        """get samfile of mapping coordinates bettween reference and assembly"""

        refA_map_path = os.path.join(
            self.dir, "minimap2", os.path.basename(self.dir) + ".paf"
        )
        if os.path.isfile(refA_map_path):

            map_coordinates = read_sam_coordinates(refA_map_path)
            map_plot = plotly_dotplot(map_coordinates)
        else:
            map_plot = None

        return map_plot

    def get_refA_dotplot(self):
        """get coverage plots"""

        dot_exists = os.path.isfile(self.refA_dotplot_path)
        output = ".".join(
            [self.sample, self.run_name, self.refdb, self.reference, "png"]
        )

        refA_dotplot_static = os.path.join(STATICFILES_DIRS[0], "refa_dotplots", output)

        if dot_exists:

            shutil.copy(
                self.refA_dotplot_path,
                os.path.join(STATICFILES_DIRS[0], "refa_dotplots", output),
            )

            self.refa_dotplot_exists = os.path.isfile(refA_dotplot_static)
            self.refa_dotplot_path = refA_dotplot_static
        else:
            self.refa_dotplot_exists = False

        # print(refA_dotplot)
        self.refa_dotplot_exists = os.path.isfile(refA_dotplot_static)
        self.refa_dotplot_path = os.path.basename(refA_dotplot_static)

    def get_ref_coverage_plot(self):
        """get coverage plots"""

        covplot_path = os.path.join(self.logdir, "ref_reads_coverage.png")
        covplot_exists = os.path.isfile(covplot_path)

        output = ".".join(
            [
                self.sample,
                self.run_name,
                self.refdb,
                self.reference,
                "coverage",
                "png",
            ]
        )

        covplot_destination = os.path.join(STATICFILES_DIRS[0], "ref_coverage", output)

        if covplot_exists:

            if not os.path.isdir(os.path.join(STATICFILES_DIRS[0], "ref_coverage")):
                os.mkdir(os.path.join(STATICFILES_DIRS[0], "ref_coverage"))

            shutil.copy(
                covplot_path,
                covplot_destination,
            ),

        # print(covplot)
        self.covplot_exists = covplot_exists
        self.covplot_path = os.path.basename(covplot_destination)
