import logging
import os
import time
from random import randint
from typing import Type

import numpy as np
import pandas as pd

from pathogen_detection.assembly_class import Assembly_class
from pathogen_detection.classification_class import Classifier
from pathogen_detection.metadata_handler import Metadata_handler
from pathogen_detection.object_classes import (
    Assembly_results,
    Contig_classification_results,
    Read_class,
    Read_classification_results,
    Remap_main,
    Remap_Target,
    Run_detail_report,
    RunCMD,
    Sample_runClass,
    Software_detail,
)
from pathogen_detection.preprocess_class import Preprocess
from pathogen_detection.remap_class import Mapping_Manager, Remapping


def get_bindir_from_binaries(binaries, key, value: str = ""):

    if value == "":
        try:
            return os.path.join(binaries["ROOT"], binaries[key]["default"], "bin")
        except KeyError:
            return ""
    else:
        try:
            return os.path.join(binaries["ROOT"], binaries[key][value], "bin")
        except KeyError:
            return ""


class RunDetail_main:
    suprun: str = None

    threads: int
    config: dict
    prefix: str
    ## input

    type: str
    r1_suffix: str
    r2_suffix: str

    r1: Type[Read_class]
    r2: Type[Read_class]

    sample: Type[Sample_runClass]
    ##  metadata
    metadata_tool: Type[Metadata_handler]
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## actions
    quality_control: bool
    sift: bool
    depletion: bool
    enrichment: bool
    assembly: bool
    classification: bool
    remapping: bool
    house_cleaning: bool

    ## methods
    preprocess_method: Software_detail
    depletion_method: Software_detail
    enrichment_method: Software_detail
    assembly_method: Software_detail
    assembly_classification_method: Software_detail
    read_classification_method: Software_detail
    remapping_method: Software_detail
    remap_manager = Mapping_Manager

    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str
    ## output content
    report: pd.DataFrame

    def __init__(self, config: dict, method_args: pd.DataFrame):

        self.logger_level_main = logging.INFO
        self.logger_level_detail = logging.CRITICAL
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(self.logger_level_main)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.runtime = 0
        self.start_time = time.perf_counter()
        self.exec_time = 0

        self.method_args = method_args
        self.config = config
        self.cmd = RunCMD(get_bindir_from_binaries(config["bin"], "PREPROCESS"))
        self.threads = config["threads"]
        self.prefix = config["prefix"]
        self.type = config["type"]
        self.suprun = self.prefix

        self.r1 = Read_class(
            config["r1"],
            config["directories"]["PREPROCESS"],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(config["bin"], "PREPROCESS"),
        )
        self.r2 = Read_class(
            config["r2"],
            config["directories"]["PREPROCESS"],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(config["bin"], "PREPROCESS"),
        )

        ###

        self.min_scaffold_length = config["assembly_contig_min_length"]
        self.minimum_coverage = int(config["minimum_coverage_threshold"])
        self.maximum_coverage = 100000
        ### metadata
        self.metadata_tool = Metadata_handler(
            self.config, sift_query=config["sift_query"], prefix=self.prefix
        )

        self.max_remap = config["max_output_number"]
        self.taxid_limit = config["taxid_limit"]

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = config["actions"]["DEPLETE"]
        self.enrichment = config["actions"]["ENRICH"]
        self.assembly = config["actions"]["ASSEMBLY"]
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAPPING"]
        self.house_cleaning = config["actions"]["CLEANING"]

        ### methods
        self.preprocess_method = Software_detail(
            "PREPROCESS",
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            "ASSEMBLY",
            method_args,
            config,
            self.prefix,
        )
        self.depletion_method = Software_detail(
            "DEPLETION",
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            "ENRICHMENT",
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            "CONTIG_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            "READ_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            "REMAPPING",
            method_args,
            config,
            self.prefix,
        )

        ### drones
        print("empty drones")
        self.depletion_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )
        self.enrichment_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )

        ### directories.
        self.root = config["directories"]["root"]
        self.filtered_reads_dir = config["directories"]["PREPROCESS"]
        self.log_dir = config["directories"]["log_dir"]

        self.report = pd.DataFrame()
        self.rclass_summary = pd.DataFrame()
        self.aclass_summary = pd.DataFrame()
        self.merged_targets = pd.DataFrame()

        ###

    def Update(self, config: dict, method_args: pd.DataFrame):

        self.method_args = method_args
        # with open(config_json) as json_file:
        #    config = json.load(json_file)

        self.config = config
        self.prefix = config["prefix"]
        self.type = config["type"]
        print("prefix:", self.prefix)
        print("type:", self.type)
        self.start_time = time.perf_counter()

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = config["actions"]["DEPLETE"]
        self.enrichment = config["actions"]["ENRICH"]
        self.assembly = config["actions"]["ASSEMBLY"]
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAPPING"]
        self.house_cleaning = config["actions"]["CLEANING"]

        self.preprocess_method = Software_detail(
            "PREPROCESS",
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            "ASSEMBLY",
            method_args,
            config,
            self.prefix,
        )

        self.depletion_method = Software_detail(
            "DEPLETION",
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            "ENRICHMENT",
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            "CONTIG_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            "READ_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            "REMAPPING",
            method_args,
            config,
            self.prefix,
        )

    def Update_exec_time(self):
        """
        Update the execution time of the pipeline.
        """
        self.exec_time = self.exec_time + time.perf_counter() - self.start_time


class Run_Deployment_Methods(RunDetail_main):
    def __init__(self, config_json: os.PathLike, method_args: pd.DataFrame):
        super().__init__(config_json, method_args)
        self.mapped_instances = []

    def deploy_QC(self):
        self.preprocess_drone = Preprocess(
            self.r1.current,
            self.r2.current,
            self.filtered_reads_dir,
            self.type,
            self.preprocess_method,
            self.r1.clean,
            self.r2.clean,
            self.threads,
            self.subsample,
            logging_level=self.logger_level_detail,
        )

        print("r1 reads: ", self.r1.get_current_fastq_read_number())
        print("r2 reads: ", self.r2.get_current_fastq_read_number())

        self.preprocess_drone.run()

    def deploy_HD(self):
        self.depletion_drone = Classifier(
            self.depletion_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.depletion_drone.run()

    def deploy_EN(self):
        self.enrichment_drone = Classifier(
            self.enrichment_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.enrichment_drone.run()

    def deploy_ASSEMBLY(self):
        self.assembly_drone = Assembly_class(
            self.r1.current,
            self.assembly_method,
            self.type,
            min_scaffold_length=self.min_scaffold_length,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.assembly_drone.run()

    def deploy_CONTIG_CLASSIFICATION(self):

        self.contig_classification_drone = Classifier(
            self.contig_classification_method,
            self.assembly_drone.assembly_file_fasta_gz,
            r2="",
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.contig_classification_drone.run()

    def deploy_READ_CLASSIFICATION(self):

        self.read_classification_drone = Classifier(
            self.read_classification_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,  #
        )

        self.read_classification_drone.run()

    def deploy_REMAPPING(self):

        self.remap_manager = Mapping_Manager(
            self.metadata_tool.remap_targets,
            self.r1,
            self.r2,
            self.remapping_method,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            self.threads,
            self.minimum_coverage,
            get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            self.logger_level_detail,
            self.house_cleaning,
        )

        self.remap_manager.run_mappings()
        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()

    def clean_unique(self):
        if self.type == "SE":
            self.clean_unique_SE()
        else:
            self.clean_unique_PE()

    def clean_unique_SE(self):
        WHERETO = os.path.dirname(self.r1.current)
        unique_reads = os.path.join(WHERETO, "unique_reads.lst")

        cmd = [
            "zgrep",
            "'^@'",
            self.r1.current,
            "|",
            "awk",
            "'{print $1}'",
            "|",
            "sed",
            "'s/^@//'",
            "|",
            "sort",
            "|",
            "uniq",
            ">",
            unique_reads,
        ]

        self.cmd.run_bash(cmd)
        if not os.path.exists(unique_reads) or os.path.getsize(unique_reads) == 0:
            self.logger.error(
                f"No unique reads found in {self.r1.current}, skipping unique read cleaning"
            )
            return

        self.r1.read_filter_inplace(self.r1.current, unique_reads)

    def clean_unique_PE(self):

        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, "common_reads.lst")

        cmd_find_common = [
            f"seqkit common -n -i {self.r1.current} {self.r1.current} | paste - - - - | cut -f1 | uniq | sed 's/^@//g' > {common_reads}"
        ]

        self.cmd.run(cmd_find_common)
        if os.path.getsize(common_reads) == 0:
            self.logger.info("No common reads found")
            return

        self.r1.read_filter_inplace(self.r1.current, common_reads)
        self.r2.read_filter_inplace(self.r2.current, common_reads)

    def trimmomatic_sort(self):
        if self.type == "SE":
            return

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            f"trimmomatic PE -phred33 -threads {self.threads} {self.r1.current} {self.r2.current} -baseout {tempfq}.fastq.gz MINLEN:20"
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq + "_1P.fastq.gz" in os.listdir(tempdir):
            os.remove(self.r1.current)
            os.remove(self.r2.current)
            os.rename(tempfq + "_1P.fastq.gz", self.r1.current)
            os.rename(tempfq + "_2P.fastq.gz", self.r2.current)


class RunMain_class(Run_Deployment_Methods):
    def __init__(self, config_json: os.PathLike, method_args: pd.DataFrame):
        super().__init__(config_json, method_args)

    def Run(self):
        print("quality control: ", self.quality_control)
        print("enrichment: ", self.enrichment)
        print("depletion: ", self.depletion)
        print("remapping: ", self.remapping)
        print("assembly: ", self.assembly)
        print("classification: ", self.classification)
        print("sift: ", self.sift)

        if self.quality_control:
            self.deploy_QC()

            self.r1.is_clean()
            self.r2.is_clean()

        if self.enrichment:
            self.deploy_EN()

            self.r1.enrich(self.enrichment_drone.classified_reads_list)
            self.r2.enrich(self.enrichment_drone.classified_reads_list)

        if self.depletion:
            self.deploy_HD()

            self.r1.deplete(self.depletion_drone.classified_reads_list)
            self.r2.deplete(self.depletion_drone.classified_reads_list)

        if self.enrichment or self.depletion or self.assembly:
            self.clean_unique()
            self.trimmomatic_sort()

        if self.assembly:
            self.deploy_ASSEMBLY()

        if self.classification:
            self.deploy_READ_CLASSIFICATION()
            self.deploy_CONTIG_CLASSIFICATION()
            self.metadata_tool.match_and_select_targets(
                self.read_classification_drone.classification_report,
                self.contig_classification_drone.classification_report,
                self.max_remap,
                self.taxid_limit,
            )
            self.aclass_summary = self.metadata_tool.aclass
            self.rclass_summary = self.metadata_tool.rclass
            self.merged_targets = self.metadata_tool.merged_targets

        if self.remapping:
            self.deploy_REMAPPING()
            self.report = self.remap_manager.report

        self.Update_exec_time()

    #### SUMMARY FUNCTIONS ####

    def Summarize(self):

        self.logger.info(f"prefix: {self.prefix}")
        with open(os.path.join(self.log_dir, self.prefix + "_latest.fofn"), "w") as f:
            f.write(self.r1.current + "\n")
            if self.type == "PE":
                f.write(self.r2.current + "\n")

        with open(os.path.join(self.log_dir, "reads_latest.stats"), "w") as f:
            f.write(f"CLEAN\t{self.r1.read_number_clean}\n")
            f.write(f"ENRICHED\t{self.r1.read_number_enriched}\n")

    def generate_output_data_classes(self):
        ### transfer to sample class
        processed_reads = (
            self.sample.r1.read_number_clean + self.sample.r2.read_number_clean
        )

        post_processed_reads = self.sample.reads_after_processing
        final_processing_reads = (
            self.r1.current_fastq_read_number() + self.r2.current_fastq_read_number()
        )

        post_percent = (int(post_processed_reads) / processed_reads) * 100
        final_processing_percent = (final_processing_reads / processed_reads) * 100

        ### transfer to assembly class / drone.

        minhit_assembly = self.aclass_summary["counts"].min()
        if not minhit_assembly or not self.aclass_summary.shape[0]:
            minhit_assembly = 0

        minhit_reads = self.rclass_summary["counts"].min()
        if np.isnan(minhit_reads):
            minhit_reads = 0

        files = list(set([t["reference"].target.file for t in self.mapped_instances]))

        self.run_detail_report = Run_detail_report(
            self.remap_manager.max_depth,
            self.remap_manager.max_depthR,
            self.remap_manager.max_gaps,
            self.remap_manager.max_prop,
            self.remap_manager.max_mapped,
            f"{processed_reads:,}",
            f"{post_processed_reads:,}",
            f"{post_percent:.2f}",
            False,
            self.sift,
            f"{self.metadata_tool.sift_report.loc[0]['removed']:,}",
            f"{final_processing_reads:,}",
            round(final_processing_percent, 3),
            self.remapping,
            self.merged_targets.taxid.nunique(),
            ", ".join(files),
        )

        self.contig_classification_results = Contig_classification_results(
            True,
            self.contig_classification_drone.classifier_method.name,
            self.contig_classification_drone.classifier_method.args,
            self.contig_classification_drone.classifier_method.db_name,
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_summary.shape[0] > 0,
        )

        self.read_classification_results = Read_classification_results(
            True,
            self.read_classification_drone.classifier_method.name,
            self.read_classification_drone.classifier_method.args,
            self.read_classification_drone.classifier_method.db_name,
            self.rclass_summary.shape[0],
            minhit_reads,
            self.rclass_summary.shape[0] > 0,
        )

        self.assembly_report = Assembly_results(
            True,
            self.assembly_drone.assembly_method.name,
            self.assembly_drone.assembly_method.args,
            self.assembly_drone.assembly_number,
            f"{self.assembly_drone.assembly_min:,}",
            f"{int(self.assembly_drone.assembly_mean):,}",
            f"{self.assembly_drone.assembly_max:,}",
            f"{int(self.min_scaffold_length):,}",
        )

        self.remap_main = Remap_main(
            True,
            self.report.shape[0],
            self.remapping_method.name,
            len(self.mapped_instances),
            self.minimum_coverage,
            self.maximum_coverage,
        )

    def export_reports(self):

        self.full_report = os.path.join(self.root, f"{self.prefix}_full_report.tsv")
        self.assembly_classification_summary = os.path.join(
            self.root, f"{self.prefix}_aclass_summary.tsv"
        )
        self.read_classification_summary = os.path.join(
            self.root, f"{self.prefix}_rclass_summary.tsv"
        )
        self.merged_classification_summary = os.path.join(
            self.root, f"{self.prefix}_mclass_summary.tsv"
        )

        ### main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )

        ### contig classification report
        self.aclass_summary.to_csv(
            self.assembly_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )

        ### read classification report
        self.rclass_summary.to_csv(
            self.read_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )

        ### merged classification report

        self.merged_targets.to_csv(
            self.merged_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )
